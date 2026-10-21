"""Only public-domain 1UBQ and a synthetic sequence; no private lab data."""
import hashlib
import io
import json
import unittest
from pathlib import Path
from v15_engine import analyze_bytes, UploadValidationError
from v15_web import create_app

ROOT=Path(__file__).resolve().parents[1]
DEMO=ROOT/'examples/1UBQ.pdb'

class PublicReleaseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.raw=DEMO.read_bytes()
        cls.app=create_app()
        cls.app.testing=True

    def test_public_fixture_integrity(self):
        receipt=json.loads((ROOT/'examples/1UBQ.receipt.json').read_text())
        self.assertEqual(hashlib.sha256(self.raw).hexdigest(),receipt['sha256'])
        self.assertEqual(receipt['license'],'CC0-1.0')

    def test_sequence_has_no_invented_structure(self):
        r=analyze_bytes('synthetic.fasta',b'>toy\nAAAAANATANPSAAAAAAAAAAAAAAA\n')
        self.assertEqual(r['status'],'needs_structure')
        self.assertEqual(r['patches'],[])
        self.assertEqual(r['residues'],[])

    def test_public_structure_cpu(self):
        r=analyze_bytes('1UBQ.pdb',self.raw,environment='soluble')
        self.assertEqual(r['status'],'ok')
        self.assertEqual(r['target']['residue_count'],76)
        self.assertGreater(len(r['patches']),0)
        self.assertTrue(all(x['plddt_verified'] is False for x in r['residues']))
        self.assertTrue(all(x['plddt'] is None for x in r['residues']))
        json.dumps(r,allow_nan=False)

    def test_local_web_endpoints(self):
        c=self.app.test_client()
        r=c.get('/api/health')
        self.assertEqual(r.status_code,200)
        self.assertTrue(r.json['samples']['1ubq'])
        self.assertNotIn('slc34a2',r.json['samples'])
        self.assertFalse(r.json['upload_persistence'])
        a=c.get('/api/sample/1ubq');self.assertEqual(a.status_code,200);a.close()
        self.assertEqual(c.get('/api/sample/slc34a2').status_code,404)
        b=c.get('/vendor/3Dmol-min.js');self.assertEqual(b.status_code,200);b.close()

    def test_web_rejects_invalid_and_accepts_fasta(self):
        c=self.app.test_client()
        bad=c.post('/api/analyze',data={'file':(io.BytesIO(b'bad'), 'bad.exe')})
        self.assertEqual(bad.status_code,400)
        ok=c.post('/api/analyze',data={'file':(io.BytesIO(b'>x\nACDEFGHIKLMNPQRSTVWY\n'),'toy.fasta')})
        self.assertEqual(ok.status_code,200)
        self.assertEqual(ok.json['status'],'needs_structure')

if __name__=='__main__':
    unittest.main()
