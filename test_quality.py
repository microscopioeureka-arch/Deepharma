import json
import tempfile
import unittest
from pathlib import Path
from scientific_quality import build_quality_report, nonfinite_paths, write_quality_report


class QualityTests(unittest.TestCase):
    def test_nested_nonfinite_rejected(self):
        self.assertEqual(nonfinite_paths({"a": [0, float("nan"), float("inf")]}),
                         ["resultado.a[1]", "resultado.a[2]"])

    def test_finite_does_not_certify_biology(self):
        r = build_quality_report({}, {"score": 100.0}, Path(__file__).parent)
        self.assertEqual(r["integridade_numerica"], "aprovada")
        self.assertEqual(r["validade_biologica"], "nao_demonstrada")
        self.assertEqual(len(r["sha256_codigo"]["simulacao.py"]), 64)

    def test_invalid_result_report_is_strict_json(self):
        r = build_quality_report({}, {"energia": float("nan")}, Path(__file__).parent)
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp) / "audit.json"
            write_quality_report(r, target)
            loaded = json.loads(target.read_text(encoding="utf-8"))
        self.assertEqual(loaded["integridade_numerica"], "reprovada")


if __name__ == "__main__":
    unittest.main()
