"""Generic structural upload analysis. No network, persistence, or design inference.

Scores are explicit heuristic priorities, never affinity or binding predictions.
"""
from __future__ import annotations

import hashlib
import io
import math
from pathlib import PurePath
import warnings as python_warnings

import Bio
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.Polypeptide import is_aa, protein_letters_3to1_extended
import numpy as np
import scipy
from scipy.spatial import cKDTree

MAX_BYTES = 15 * 1024 * 1024
MAX_ATOMS = 30000
RADII = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80,
         "SE": 1.90, "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98}
MAX_ASA = dict(zip("ARNDCQEGHILKMFPSTWYV", [129,274,195,193,167,225,223,104,224,197,201,236,224,240,159,155,172,285,263,174]))
KD = dict(zip("ACDEFGHIKLMNPQRSTVWY", [1.8,2.5,-3.5,-3.5,2.8,-0.4,-3.2,4.5,-3.9,3.8,1.9,-3.5,-1.6,-3.5,-4.5,-0.8,-0.7,4.2,-0.9,-1.3]))
WEIGHTS = {"exposure": .25, "sidechain_area": .15, "normal_coherence": .10,
           "tertiary_support": .25, "low_hydrophobic_risk": .15, "low_glycan_risk": .10}
SIDECHAIN_ATOMS = {
    "A":"CB", "R":"CB CG CD NE CZ NH1 NH2", "N":"CB CG OD1 ND2",
    "D":"CB CG OD1 OD2", "C":"CB SG", "Q":"CB CG CD OE1 NE2",
    "E":"CB CG CD OE1 OE2", "G":"", "H":"CB CG ND1 CD2 CE1 NE2",
    "I":"CB CG1 CG2 CD1", "L":"CB CG CD1 CD2", "K":"CB CG CD CE NZ",
    "M":"CB CG SD CE", "F":"CB CG CD1 CD2 CE1 CE2 CZ", "P":"CB CG CD",
    "S":"CB OG", "T":"CB OG1 CG2", "W":"CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2",
    "Y":"CB CG CD1 CD2 CE1 CE2 CZ OH", "V":"CB CG1 CG2"}


class UploadValidationError(ValueError):
    """Invalid input; web callers should translate to an HTTP 400 response."""


def _unit(v):
    length = float(np.linalg.norm(v))
    return v / length if length > 1e-9 else np.zeros(3)


def _aa(residue):
    return protein_letters_3to1_extended.get(residue.resname.upper(), "X")


def _sequons(sequence):
    return [i for i in range(len(sequence)-2)
            if sequence[i] == "N" and sequence[i+1] != "P" and sequence[i+2] in "ST"]


def _frame(xyz):
    """Intrinsic frame makes the finite sphere quadrature rigid-motion equivariant."""
    centered = xyz - xyz.mean(axis=0)
    x = _unit(centered[np.argmax(np.sum(centered**2, axis=1))])
    perpendicular = centered - np.outer(centered @ x, x)
    y = _unit(perpendicular[np.argmax(np.sum(perpendicular**2, axis=1))])
    if np.linalg.norm(y) < .5:
        raise UploadValidationError("Degenerate/collinear atomic geometry")
    return np.stack((x, y, np.cross(x, y)))


def _surface(xyz, radii, points):
    z = 1 - 2*(np.arange(points)+.5)/points
    phi = np.arange(points) * math.pi * (3-math.sqrt(5))
    sphere = np.stack((np.sqrt(1-z*z)*np.cos(phi), np.sqrt(1-z*z)*np.sin(phi), z), axis=1) @ _frame(xyz)
    expanded = radii + 1.4
    tree = cKDTree(xyz)
    counts = tree.query_ball_point(xyz, expanded + float(expanded.max()), return_length=True)
    if counts.max() > 1024 or counts.sum() > 6000000:
        raise UploadValidationError("Unphysical local atom density; surface resource limit exceeded")
    neighborhoods = tree.query_ball_point(xyz, expanded + float(expanded.max()))
    area = np.zeros(len(xyz))
    vectors = np.zeros_like(xyz)
    for i, neighbors in enumerate(neighborhoods):
        js = np.array([j for j in neighbors if j != i], dtype=int)
        dots = xyz[i] + expanded[i]*sphere
        exposed = np.ones(points, dtype=bool)
        # Chunk neighbors to bound memory for deliberately dense uploaded structures.
        for start in range(0, len(js), 256):
            part = js[start:start+256]
            delta = dots[:, None, :] - xyz[part][None, :, :]
            exposed &= ~np.any(np.einsum("ijk,ijk->ij", delta, delta) < expanded[part]**2 - 1e-9, axis=1)
            if not exposed.any():
                break
        count = int(exposed.sum())
        if count:
            area[i] = 4*math.pi*expanded[i]**2*count/points
            vectors[i] = sphere[exposed].mean(axis=0)
    return area, vectors, tree


def _clearance(centers, envelope_radii, xyz, radii, tree):
    minimum = float("inf")
    collisions = 0
    for center, radius in zip(centers, envelope_radii):
        near = tree.query_ball_point(center, float(radius + radii.max() + 1.0))
        if near:
            gap = np.linalg.norm(xyz[near]-center, axis=1) - radii[near] - radius
            minimum = min(minimum, float(gap.min()))
            collisions += int(np.sum(gap < -.05))
    return collisions == 0, collisions, round(minimum, 4) if math.isfinite(minimum) else None


def _sweep(center, direction, binder, xyz, radii, tree):
    # Overlapping spheres sample a tapered straight approach, including distal body.
    # These envelopes are conservative geometric heuristics, not fitted antibody poses.
    if binder == "VHH":
        t = np.array([4., 8., 12., 17., 23., 30., 38., 46.])
        r = np.array([2., 4., 6., 8., 10., 11., 11., 11.])
    else:
        t = np.array([4., 8., 13., 19., 26., 34., 44., 55., 66.])
        r = np.array([2., 4., 7., 10., 13., 16., 20., 22., 22.])
    if np.linalg.norm(direction) < .5:
        return False, 0, None
    # Continue the broad distal body to bulk solvent; distant chains must not be missed.
    far = float(np.max((xyz-center) @ direction)) + float(r[-1]) + 2.
    if far > t[-1]:
        extra = np.arange(t[-1]+r[-1]*.8, far+r[-1]*.8, r[-1]*.8)
        t = np.concatenate((t,extra))
        r = np.concatenate((r,np.full(len(extra),r[-1])))
    return _clearance(center + t[:, None]*direction, r, xyz, radii, tree)


def _directions(i, normals, neighbors):
    average = _unit(normals[neighbors[i]].sum(axis=0)) if neighbors[i] else normals[i]
    choices = [normals[i], average, _unit(normals[i]+average)]
    unique = []
    for d in choices:
        if np.linalg.norm(d) > .5 and not any(np.dot(d, p) > .999 for p in unique):
            unique.append(d)
    return unique


def _normalize(residues, chain):
    lines, numbering = [], []
    serial = 0
    for i, residue in enumerate(residues, 1):
        numbering.append({"resi": i, "chain_original": chain.id, "resseq": int(residue.id[1]),
                          "icode": residue.id[2].strip(), "aa": _aa(residue),
                          "resname_original": residue.resname})
        for atom in residue:
            if atom.element.upper() in {"H", "D"}:
                continue
            serial += 1
            x, y, z = (float(v) for v in atom.coord)
            if any(v < -999.999 or v > 9999.999 for v in (x,y,z)):
                raise UploadValidationError("Coordinates exceed normalized PDB field limits")
            occupancy = atom.occupancy if atom.occupancy is not None else 1.0
            lines.append(f"ATOM  {serial:5d} {atom.fullname:4s} {residue.resname:>3s} A{i:4d}    "
                         f"{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{atom.bfactor:6.2f}          {atom.element:>2s}  ")
    return "\n".join(lines + ["TER", "END", ""]), numbering, serial


def analyze_bytes(filename: str, raw: bytes, target_chain: str = "auto",
                  environment: str = "detergent", binder_type: str = "VHH",
                  settings: dict | None = None) -> dict:
    """Analyze an uploaded PDB/mmCIF/FASTA entirely in local memory.

    Raises UploadValidationError for invalid/oversized inputs. Settings accepts
    sphere_points (48..96), max_patches (3..8); no settings provide a certified
    membrane orientation or reinterpret B factors as confidence.
    """
    if not isinstance(filename, str) or not filename.strip():
        raise UploadValidationError("filename must be a nonempty string")
    if not isinstance(raw, bytes) or not raw or len(raw) > MAX_BYTES:
        raise UploadValidationError("Upload must be nonempty bytes, maximum 15 MiB")
    if not isinstance(target_chain, str):
        raise UploadValidationError("target_chain must be an author chain ID or auto")
    if environment not in ("detergent", "soluble", "bilayer"):
        raise UploadValidationError("environment must be detergent, soluble, or bilayer")
    if binder_type not in ("VHH", "Fab"):
        raise UploadValidationError("binder_type must be VHH or Fab")
    settings = {} if settings is None else settings
    if not isinstance(settings, dict):
        raise UploadValidationError("settings must be a dictionary")
    points, max_patches = settings.get("sphere_points", 64), settings.get("max_patches", 6)
    if type(points) is not int or not 48 <= points <= 96:
        raise UploadValidationError("sphere_points must be an integer from 48 to 96")
    if type(max_patches) is not int or not 3 <= max_patches <= 8:
        raise UploadValidationError("max_patches must be an integer from 3 to 8")
    suffix = PurePath(filename).suffix.lower()
    fmt = {".pdb":"pdb", ".cif":"mmcif", ".mmcif":"mmcif", ".fa":"fasta", ".fasta":"fasta", ".faa":"fasta"}.get(suffix)
    if fmt is None:
        raise UploadValidationError("Supported extensions: .pdb, .cif, .mmcif, .fa, .fasta, .faa")
    try:
        content = raw.decode("utf-8-sig")
    except UnicodeError as exc:
        raise UploadValidationError("Upload must be UTF-8 text") from exc
    if "\x00" in content:
        raise UploadValidationError("NUL bytes are not allowed")
    warn = ["EXPLORATORY_GEOMETRY_NOT_BINDING_OR_AFFINITY", "GLYCAN_OCCUPANCY_UNKNOWN"]
    if set(settings) - {"sphere_points", "max_patches"}:
        warn.append("UNSUPPORTED_SETTINGS_IGNORED")
    env = {"mode": environment, "label": {"detergent":"DDM/GDN micelle context (proxy only)", "soluble":"Soluble; full uploaded surface", "bilayer":"Bilayer; orientation unregistered"}[environment],
           "inferred_membrane": False, "orientation_status": "UNKNOWN" if environment != "soluble" else "NOT_APPLICABLE",
           "detergent_risk_model": "sequence_window_and_3D_hydrophobic_cluster_proxy_not_micelle_boundary" if environment == "detergent" else "not_applied"}
    if environment != "soluble":
        warn += ["MEMBRANE_ORIENTATION_UNKNOWN", "MEMBRANE_ACCESS_NOT_CERTIFIED"]
    out = {"status":"ok", "source":{"filename":filename, "sha256":hashlib.sha256(raw).hexdigest(), "format":fmt, "parsed_chain_names":[]},
           "environment":env, "residues":[], "patches":[], "warnings":warn,
           "qc":{"engine_version":"1.5.0", "tools":{"biopython":Bio.__version__, "numpy":np.__version__, "scipy":scipy.__version__},
                 "persistent_upload":False, "score_weights":WEIGHTS.copy(), "assumptions":[
                     "Only first model; highest-occupancy selected altloc; author mmCIF IDs.",
                     "All uploaded non-water heavy atoms obstruct surface and approach; no biological assembly expansion.",
                     "B factors are not verified pLDDT; confidence is unavailable and contributes a fixed 0.85 uncertainty multiplier.",
                     "Approximate Shrake-Rupley SASA, 1.4 A probe; SASA is not an antibody-access guarantee.",
                     "Coarse tapered overlapping-sphere straight approach, not docking or binding validation.",
                     "No membrane plane or exact micelle boundary is inferred from coordinates."]}}
    if fmt == "fasta":
        lines = [line.strip() for line in content.splitlines() if line.strip()]
        if not lines or not lines[0].startswith(">") or sum(line.startswith(">") for line in lines) != 1:
            raise UploadValidationError("FASTA must contain exactly one header and sequence")
        seq = "".join(lines[1:]).upper()
        if not 20 <= len(seq) <= 10000 or set(seq)-set("ACDEFGHIKLMNPQRSTVWYXBZUO"):
            raise UploadValidationError("Protein FASTA requires 20..10000 amino acids, no gaps or stop symbols")
        if target_chain != "auto":
            raise UploadValidationError("FASTA has no selectable structural chain")
        out.update(status="needs_structure", reason="structure_prediction_not_run")
        out["target"] = {"name":lines[0][1:] or PurePath(filename).stem, "chain_original":None,
                         "sequence":seq, "residue_count":len(seq), "atom_count":0, "pdb":"", "numbering":[],
                         "glyco_sequons":[{"resi":i+1, "motif":seq[i:i+3], "occupancy":"unknown"} for i in _sequons(seq)]}
        out["qc"]["geometry_computed"] = False
        return out
    try:
        parser = PDBParser(PERMISSIVE=False, QUIET=False) if fmt == "pdb" else MMCIFParser(QUIET=False, auth_chains=True, auth_residues=True)
        with python_warnings.catch_warnings(record=True) as caught:
            python_warnings.simplefilter("always")
            structure = parser.get_structure("upload", io.StringIO(content))
        out["qc"]["parser_warning_count"] = len(caught)
        models = list(structure)
        if not models:
            raise ValueError("no models")
        all_atoms = sum(len(residue.get_unpacked_list()) for model in models for chain in model for residue in chain.get_unpacked_list())
        if all_atoms > MAX_ATOMS:
            raise UploadValidationError("Upload exceeds 30000 atoms across all models/alternate locations")
        model = models[0]
        chains = list(model)
        out["source"]["parsed_chain_names"] = [chain.id for chain in chains]
        candidates = [(chain, [r for r in chain if is_aa(r, standard=False) and "CA" in r]) for chain in chains]
        if target_chain == "auto":
            chain, residues = max(candidates, key=lambda item: len(item[1]))
        else:
            selected = [(c,r) for c,r in candidates if c.id == target_chain]
            if not selected:
                raise UploadValidationError(f"Unknown author chain: {target_chain!r}")
            chain, residues = selected[0]
        if len(residues) < 20 or len(residues) > 9999:
            raise UploadValidationError("Selected protein chain needs 20..9999 residues with CA atoms")
        atoms = [a for c in chains for r in c if r.id[0] != "W" for a in r if a.element.upper() not in {"H", "D"}]
        xyz = np.array([a.coord for a in atoms], dtype=float)
        if not np.isfinite(xyz).all() or not all(np.isfinite(a.bfactor) and (a.occupancy is None or np.isfinite(a.occupancy)) for a in atoms):
            raise UploadValidationError("Non-finite atomic coordinates/metadata")
    except UploadValidationError:
        raise
    except Exception as exc:
        raise UploadValidationError(f"Malformed {fmt} structure: {str(exc)[:180]}") from exc
    if len(models) > 1:
        warn.append("FIRST_MODEL_ONLY")
    if caught:
        warn.append("PARSER_WARNINGS_REVIEW_INPUT")
    missing_ca = sum(is_aa(r, standard=False) and "CA" not in r for r in chain)
    if missing_ca:
        warn.append("PROTEIN_RESIDUES_WITHOUT_CA_OMITTED_FROM_TARGET")
    if any(a.is_disordered() for a in atoms):
        warn.append("ALTERNATE_LOCATIONS_SELECTED_NOT_ENSEMBLE")
    if any(a.element.upper() not in RADII for a in atoms):
        warn.append("UNKNOWN_ELEMENT_RADIUS_DEFAULT_1_7_A")
    radii = np.array([RADII.get(a.element.upper(), 1.7) for a in atoms])
    pdb, numbering, target_atoms = _normalize(residues, chain)
    sequence = "".join(_aa(r) for r in residues)
    glyco = _sequons(sequence)
    out["target"] = {"name":PurePath(filename).stem, "chain_original":chain.id, "sequence":sequence,
                     "residue_count":len(residues), "atom_count":target_atoms, "pdb":pdb, "numbering":numbering,
                     "glyco_sequons":[{"resi":i+1, "motif":sequence[i:i+3], "occupancy":"unknown"} for i in glyco]}
    area, vectors, tree = _surface(xyz, radii, points)
    indices = {id(atom):i for i,atom in enumerate(atoms)}
    n = len(residues)
    ca = np.array([r["CA"].coord for r in residues], dtype=float)
    normals, centers = np.zeros((n,3)), np.zeros((n,3))
    total, sc = np.zeros(n), np.zeros(n)
    completeness = np.ones(n, dtype=bool)
    for i,r in enumerate(residues):
        js = [indices[id(a)] for a in r if id(a) in indices]
        side = [indices[id(a)] for a in r if id(a) in indices and a.name not in {"N","CA","C","O","OXT"}]
        total[i], sc[i] = area[js].sum(), area[side].sum()
        use = side if sc[i] > 1 else js
        normals[i] = _unit((vectors[use]*area[use,None]).sum(axis=0))
        if area[use].sum() > 0:
            centers[i] = np.average(xyz[use] + (radii[use]+1.4)[:,None]*vectors[use], axis=0, weights=area[use])
        else:
            centers[i] = ca[i]
        required = {"N","CA","C","O"} | set(SIDECHAIN_ATOMS.get(sequence[i], "CB").split())
        if r.resname == "MSE":
            required.discard("SD")
            required.add("SE")
        completeness[i] = all(name in r for name in required)
    chain_breaks = [i+1 for i in range(n-1) if "C" not in residues[i] or "N" not in residues[i+1]
                    or np.linalg.norm(residues[i]["C"].coord-residues[i+1]["N"].coord) > 2.2]
    if chain_breaks:
        warn.append("OBSERVED_SEQUENCE_HAS_CHAIN_BREAKS_MISSING_RESIDUES_NOT_REBUILT")
    if not completeness.all():
        warn.append("INCOMPLETE_RESIDUES_EXCLUDED_FROM_HOTSPOTS")
    rsa = np.clip(total / np.array([MAX_ASA.get(a, 200) for a in sequence]), 0, 1)
    neighborhood = cKDTree(ca).query_ball_point(ca, 10.)
    packing_neighborhood = cKDTree(ca).query_ball_point(ca, 12.)
    # Exclude eight sequence neighbors on either side: a lone helix must not
    # masquerade as tertiary packing just because its turns are locally close.
    contacts = np.array([sum(abs(i-j)>8 for j in js) for i,js in enumerate(packing_neighborhood)])
    packing = np.clip(contacts/12.,0,1)
    terminal = np.array([min(i,n-1-i)<15 for i in range(n)])
    terminal_low_support = terminal & (contacts<6)
    hydro = np.array([(KD.get(a,0)+4.5)/9. for a in sequence])
    local_hydro = np.array([np.average(hydro[js], weights=total[js]+1e-6) for js in neighborhood])
    risk = np.clip(.5*hydro + .5*local_hydro, 0, 1)
    # A sequence+3D belt proxy, not an inferred plane or a TM annotation.
    kd_values = np.array([KD.get(a,0) for a in sequence])
    window_hydro = np.array([kd_values[max(0,i-9):min(n,i+10)].mean() for i in range(n)])
    tm_seed = (window_hydro>=1.2) & (np.arange(n)>=9) & (np.arange(n)<n-9)
    micelle_proximity = np.zeros(n)
    for j in np.flatnonzero(tm_seed):
        strength = float(np.clip((window_hydro[j]-.5)/1.5,0,1))
        distance = np.linalg.norm(ca-ca[j],axis=1)
        micelle_proximity = np.maximum(micelle_proximity,strength*np.exp(-np.maximum(distance-6,0)/8))
    if environment != "soluble":
        risk = np.maximum(risk,micelle_proximity)
    glycan_near = np.zeros(n, dtype=bool)
    for i in glyco:
        glycan_near |= np.linalg.norm(ca-ca[i], axis=1) < 10.
    graph = [[] for _ in residues]
    exposed = (rsa >= .12) & (total >= 15) & (np.linalg.norm(normals,axis=1) > .5)
    for i in range(n):
        if not exposed[i]:
            continue
        for j in neighborhood[i]:
            if j <= i or not exposed[j] or np.dot(normals[i],normals[j]) < .15 or np.linalg.norm(centers[i]-centers[j]) > 11:
                continue
            start, end = centers[i]+normals[i]*2, centers[j]+normals[j]*2
            samples = start + np.linspace(0,1,7)[:,None]*(end-start)
            clear, _, _ = _clearance(samples, np.full(7,.4), xyz, radii, tree)
            if clear:
                graph[i].append(j)
                graph[j].append(i)
    access = np.zeros(n, dtype=int)
    for i in range(n):
        if exposed[i]:
            access[i] = sum(_sweep(centers[i],d,binder_type,xyz,radii,tree)[0] for d in _directions(i,normals,graph))
    # Rank only residues that pass a physical approach probe; never use motif Asn.
    raw_score = .25*rsa + .15*np.clip(sc/100,0,1) + .10 + .25*packing + .15*(1-risk if environment != "soluble" else np.ones(n)) + .1*(~glycan_near)
    raw_score *= .85
    raw_score *= .6+.4*packing
    raw_score[terminal_low_support] *= .6
    if environment != "soluble":
        raw_score *= 1-.55*micelle_proximity
    eligible = exposed & (access > 0) & completeness
    eligible[glyco] = False
    scores = np.where(eligible, raw_score, 0.)
    for i in range(n):
        rw = ["PLDDT_UNVERIFIED"]
        if i in glyco:
            rw.append("GLYCO_SEQUON_ASN_OCCUPANCY_UNKNOWN_EXCLUDED_HOTSPOT")
        elif glycan_near[i]:
            rw.append("NEAR_POTENTIAL_GLYCAN")
        if not completeness[i]:
            rw.append("INCOMPLETE_BACKBONE_OR_SIDECHAIN")
        if contacts[i] < 6:
            rw.append("LOW_NONLOCAL_PACKING_SUPPORT_NOT_DISORDER_PROOF")
        if terminal_low_support[i]:
            rw.append("TERMINAL_LOW_STRUCTURAL_SUPPORT")
        if environment != "soluble" and micelle_proximity[i] >= .55:
            rw.append("TM_LIKE_SEQUENCE_OR_3D_HYDROPHOBIC_CLUSTER_PROXIMITY")
        if environment != "soluble" and risk[i] >= .60:
            rw.append("HYDROPHOBIC_CONTEXT_RISK_UNCERTIFIED")
        out["residues"].append({"resi":i+1, "aa":sequence[i], "ca":ca[i].tolist(),
             "rsa_proxy":round(float(rsa[i]),5), "sidechain_sas_A2":round(float(sc[i]),4),
             "plddt":None, "plddt_verified":False, "glyco_sequon":i in glyco,
             "hydrophobic_risk":round(float(risk[i]),5), "accessible_probe_count":int(access[i]),
             "nonlocal_contact_count":int(contacts[i]), "structural_support":round(float(packing[i]),5),
             "terminal":bool(terminal[i]), "terminal_low_support":bool(terminal_low_support[i]),
             "sequence_hydropathy_19":round(float(window_hydro[i]),5),
             "micelle_proximity_proxy":round(float(micelle_proximity[i]),5),
             "access_status":"EXPOSED_APPROACH_CLEAR" if access[i] else ("APPROACH_BLOCKED" if exposed[i] else "LOW_EXPOSURE"),
             "color_score":round(float(scores[i]),5), "warnings":rw})
    candidates_patches = []
    for anchor in sorted(range(n), key=lambda i:(-scores[i], i)):
        if not eligible[anchor]:
            continue
        # Grow solely over solvent-connected graph edges, within a compact surface radius.
        members, frontier = {anchor}, [anchor]
        while frontier and len(members) < 12:
            node = frontier.pop(0)
            for j in sorted(graph[node], key=lambda j:(-scores[j],j)):
                if j not in members and np.linalg.norm(centers[j]-centers[anchor]) <= 12 and np.dot(normals[j],normals[anchor]) > .25:
                    members.add(j)
                    frontier.append(j)
                    if len(members) >= 12:
                        break
        if len(members) < 3:
            continue
        hot = sorted((i for i in members if eligible[i]), key=lambda i:(-scores[i],i))[:4]
        if len(hot) < 2:
            continue
        center = np.average(centers[hot], axis=0, weights=total[hot])
        directions = [_unit(normals[hot].sum(axis=0))] + _directions(anchor,normals,graph)
        passes = []
        min_collisions = None
        for direction in directions:
            clear, collisions, gap = _sweep(center,direction,binder_type,xyz,radii,tree)
            # Every hotspot must share this one outward direction and clear tip corridor.
            coherent = all(np.dot(normals[i],direction) >= .25 for i in hot)
            tips = [centers[i]+direction*t for i in hot for t in (3.,5.,8.)]
            tips_ok, tip_collisions, _ = _clearance(tips, np.full(len(tips),1.4),xyz,radii,tree)
            count = collisions + tip_collisions
            min_collisions = count if min_collisions is None else min(min_collisions,count)
            if clear and coherent and tips_ok:
                passes.append((direction,gap))
        physical_pass = bool(passes)
        coherence = float(np.linalg.norm(normals[hot].mean(axis=0)))
        components = {"exposure":float(rsa[hot].mean()), "sidechain_area":float(np.clip(sc[hot]/100,0,1).mean()),
                      "normal_coherence":coherence, "tertiary_support":float(packing[hot].mean()),
                      "low_hydrophobic_risk":float(1-risk[hot].mean()) if environment != "soluble" else 1.,
                      "low_glycan_risk":float((~glycan_near[hot]).mean())}
        score = .85*sum(WEIGHTS[k]*v for k,v in components.items()) if physical_pass else 0.
        support_multiplier = .6+.4*components['tertiary_support']
        terminal_multiplier = 1-.4*float(terminal_low_support[hot].mean())
        micelle_multiplier = 1-.55*float(micelle_proximity[hot].mean()) if environment != "soluble" else 1.
        score *= support_multiplier*terminal_multiplier*micelle_multiplier
        pw = ["UNVERIFIED_CONFIDENCE", "COARSE_ENVELOPE_NOT_DOCKED_ANTIBODY"]
        hold = not physical_pass or environment == "bilayer"
        if float(contacts[hot].mean()) < 6 or float((contacts[hot]<6).mean()) > .5:
            pw.append("LOW_NONLOCAL_PACKING_SUPPORT_NOT_DISORDER_PROOF")
            hold = True
        if terminal_low_support[hot].any():
            pw.append("TERMINAL_LOW_STRUCTURAL_SUPPORT")
            hold = True
        if not physical_pass:
            pw.append("NO_COHERENT_COLLISION_FREE_APPROACH")
        if any(glycan_near[hot]):
            pw.append("POTENTIAL_GLYCAN_NEAR_HOTSPOTS_OCCUPANCY_UNKNOWN")
            hold = True
        if environment != "soluble":
            pw.append("MEMBRANE_MICELLE_ACCESS_UNCERTIFIED")
            if float(risk[hot].mean()) >= .60:
                hold = True
                pw.append("HIGH_HYDROPHOBIC_CONTEXT_RISK")
            if float(micelle_proximity[hot].mean()) >= .55:
                pw.append("TM_LIKE_SEQUENCE_OR_3D_HYDROPHOBIC_CLUSTER_PROXIMITY")
                hold = True
        candidates_patches.append({"id":"", "anchor":anchor+1, "members":[i+1 for i in sorted(members)],
             "hotspots":[i+1 for i in hot], "center":center.tolist(), "score":round(score,5),
             "status":"HOLD_REVIEW" if hold else "DESIGN_EXPLORATION", "binder_type":binder_type,
             "rationale":["Connected exposed surface graph with outward-normal and solvent-segment checks.",
                          "Shared tapered-envelope approach passed." if physical_pass else "Shared approach failed; score suppressed.",
                          f"Mean nonlocal CA contacts within 12 A (sequence separation >8): {contacts[hot].mean():.2f}; support and terminal penalties applied.",
                          "Human review and upstream design exploration only; no affinity or hit-rate prediction."],
             "warnings":pw, "approach_direction":passes[0][0].tolist() if passes else None,
             "access_checks":{"physical_pass":physical_pass,"coherent_direction_count":len(passes),
                              "minimum_collision_count":min_collisions, "sampled_min_clearance_A":passes[0][1] if passes else None,
                              "envelope_model":"tapered_overlapping_spheres_v1", "surface_graph_connected":True,
                              "score_components":components,"uncertainty_multiplier":.85,
                              "mean_nonlocal_contacts":round(float(contacts[hot].mean()),4),
                              "terminal_low_support_count":int(terminal_low_support[hot].sum()),
                              "mean_micelle_proximity_proxy":round(float(micelle_proximity[hot].mean()),5),
                              "support_multiplier":support_multiplier,"terminal_multiplier":terminal_multiplier,
                              "micelle_multiplier":micelle_multiplier}})
    # Greedy diversity: limited overlap, different spatial anchors; no forced three-patch result.
    for patch in sorted(candidates_patches, key=lambda p:(-p["score"], p["anchor"])):
        if any(len(set(patch["members"]) & set(p["members"]))/min(len(patch["members"]),len(p["members"])) > .5
               or np.linalg.norm(np.array(patch["center"])-p["center"]) < 9 for p in out["patches"]):
            continue
        patch["id"] = f"P{len(out['patches'])+1:02d}"
        out["patches"].append(patch)
        if len(out["patches"]) >= max_patches:
            break
    if len(out["patches"]) < 3:
        warn.append("FEWER_THAN_THREE_SUPPORTED_PATCHES_NO_FORCED_RECOMMENDATIONS")
    out["qc"].update(geometry_computed=True, sphere_points=points, solvent_probe_A=1.4,
         uploaded_atom_count=all_atoms, occluding_heavy_atom_count=len(atoms), model_count=len(models),
         chain_breaks_after_resi=chain_breaks, incomplete_residue_count=int((~completeness).sum()),
         target_missing_ca_count=missing_ca, exposed_residue_count=int(exposed.sum()),
         approach_clear_residue_count=int((access>0).sum()), surface_graph_edge_count=sum(map(len,graph))//2,
         sasa_total_target_A2=round(float(total.sum()),3),
         structural_support_definition="CA neighbors within 12 A excluding sequence separation <=8; min(count/12,1); low support is not proof of disorder",
         micelle_proxy_definition="19-residue Kyte-Doolittle mean >=1.2 seeds; strength clip((mean-0.5)/1.5,0,1); max spatial decay exp(-max(distance-6,0)/8); no orientation inference",
         tm_like_seed_count=int(tm_seed.sum()),
         rsa_definition="per-residue total heavy-atom SASA / Tien maximum ASA, clipped 0..1; approximate proxy",
         normalized_pdb_precision_A=.001,
         surface_graph_edges=[[i+1,j+1] for i in range(n) for j in graph[i] if i<j])
    return out
