"""Public-fixture regression tests for CLI, scientific nonclaims and local API."""
import io
import json
import subprocess
import sys
import unittest
from pathlib import Path

from Bio.PDB import MMCIFIO, PDBParser

from v15_engine import UploadValidationError, analyze_bytes
from v15_engine.core import MAX_BYTES
from v15_web import create_app

ROOT = Path(__file__).resolve().parents[1]
DEMO = ROOT / "examples" / "1UBQ.pdb"


class PublicResearchContracts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.raw = DEMO.read_bytes()

    def test_cli_human_summary_never_claims_affinity(self):
        p = subprocess.run(
            [sys.executable, "-m", "v15_engine", str(DEMO),
             "--environment", "soluble", "--summary"],
            text=True, capture_output=True, check=True, timeout=30
        )
        self.assertIn("NOT validated binding or affinity", p.stdout)
        self.assertIn("P01:", p.stdout)
        self.assertIn("Residues: 76", p.stdout)
        self.assertEqual(p.stderr, "")

    def test_cli_default_json_is_machine_readable(self):
        p = subprocess.run(
            [sys.executable, "-m", "v15_engine", str(DEMO),
             "--environment", "soluble"],
            text=True, capture_output=True, check=True, timeout=30
        )
        obj = json.loads(p.stdout)
        self.assertEqual(obj["status"], "ok")
        self.assertEqual(obj["target"]["residue_count"], 76)
        self.assertGreaterEqual(len(obj["patches"]), 1)

    def test_repeated_geometry_is_deterministic(self):
        a = analyze_bytes("1UBQ.pdb", self.raw, environment="soluble")
        b = analyze_bytes("1UBQ.pdb", self.raw, environment="soluble")
        self.assertEqual(a["source"]["sha256"], b["source"]["sha256"])
        self.assertEqual(a["patches"], b["patches"])
        self.assertTrue(all(not r["plddt_verified"] for r in a["residues"]))

    def test_mmcif_round_trip_on_public_structure(self):
        structure = PDBParser(QUIET=True).get_structure("public", io.StringIO(self.raw.decode()))
        writer = MMCIFIO()
        writer.set_structure(structure)
        text = io.StringIO()
        writer.save(text)
        result = analyze_bytes("1UBQ.cif", text.getvalue().encode("utf-8"), environment="soluble")
        self.assertEqual(result["status"], "ok")
        self.assertEqual(result["target"]["residue_count"], 76)
        self.assertEqual(result["target"]["chain_original"], "A")
        self.assertEqual(len(result["target"]["numbering"]), 76)

    def test_unknown_chain_is_rejected(self):
        with self.assertRaises(UploadValidationError):
            analyze_bytes("1UBQ.pdb", self.raw, target_chain="nonexistent")

    def test_rejects_untrusted_input(self):
        for filename, data in (
            ("toy.exe", b"hello"),
            ("toy.pdb", b"\0bad"),
            ("toy.pdb", b""),
            ("toy.fasta", b">x\nAAA\n"),
            ("toy.fasta", b">x\n" + b"A" * 20 + b"*\n"),
        ):
            with self.subTest(filename=filename, contents=data[:10]):
                with self.assertRaises(UploadValidationError):
                    analyze_bytes(filename, data)
        with self.assertRaises(UploadValidationError):
            analyze_bytes("oversize.pdb", b"x" * (MAX_BYTES + 1))

    def test_sequence_without_structure_never_generates_hotspots(self):
        result = analyze_bytes("toy.fasta", b">toy\nAAAAANATANPSAAAAAAAAAAAAAAA\n")
        self.assertEqual(result["status"], "needs_structure")
        self.assertFalse(result["patches"])
        self.assertFalse(result["residues"])
        self.assertFalse(result["target"]["pdb"])

    def test_bilayer_without_orientation_is_held(self):
        r = analyze_bytes("1UBQ.pdb", self.raw, environment="bilayer")
        self.assertEqual(r["environment"]["orientation_status"], "UNKNOWN")
        self.assertFalse(r["environment"]["inferred_membrane"])
        self.assertGreater(len(r["patches"]), 0)
        self.assertTrue(all(p["status"] == "HOLD_REVIEW" for p in r["patches"]))

    def test_local_server_accepts_public_demo_and_rejects_remote_origin(self):
        app = create_app()
        app.testing = True
        c = app.test_client()
        health = c.get("/api/health")
        self.assertEqual(health.status_code, 200)
        self.assertFalse(health.json["upload_persistence"])
        self.assertEqual(c.get("/api/sample/slc34a2").status_code, 404)
        self.assertEqual(c.get("/api/health", headers={"Host": "evil.example"}).status_code, 403)
        form = {"file": (io.BytesIO(b">toy\n" + b"A" * 30 + b"\n"), "toy.fasta")}
        resp = c.post("/api/analyze", data=form,
                      headers={"Origin": "https://evil.example", "Host": "localhost"})
        self.assertEqual(resp.status_code, 403)

    def test_port_flag_is_explicit_and_loopback_only(self):
        h = subprocess.run([sys.executable, "-m", "v15_web", "--help"],
                           text=True, capture_output=True, check=True, timeout=15)
        self.assertIn("--port", h.stdout)
        invalid = subprocess.run([sys.executable, "-m", "v15_web", "--port", "0"],
                                 text=True, capture_output=True, timeout=15)
        self.assertEqual(invalid.returncode, 2)
        self.assertIn("port must be between", invalid.stderr)


if __name__ == "__main__":
    unittest.main()
