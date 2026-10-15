"""HTTP contract tests with simulation stubbed. No scientific benchmark is run."""
import importlib
import os
import sys
import tempfile
import types
import unittest
from unittest.mock import patch
from fastapi.testclient import TestClient


class ServerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.env = patch.dict(os.environ, APP_PASSWORD="sênha-teste", SESSION_SECRET="test-only-secret")
        cls.env.start()
        fake = types.ModuleType("simulacao")
        fake.AMINOACIDOS = {"ALA": {}}
        fake.ENVIRONMENTS = {"citoplasma": {}}
        fake.COMPOSTOS_REAIS = {"Paracetamol": {}, "Cafeina": {}, "Cloroquina": {}}
        cls.modules = patch.dict(sys.modules, simulacao=fake)
        cls.modules.start()
        cls.server = importlib.import_module("server")

    @classmethod
    def tearDownClass(cls):
        sys.modules.pop("server", None)
        cls.modules.stop()
        cls.env.stop()

    def setUp(self):
        self.server.JOBS.clear()
        self.client = TestClient(self.server.app)

    def tearDown(self):
        self.client.close()

    def login(self):
        response = self.client.post("/login", data={"senha": "sênha-teste"}, follow_redirects=False)
        self.assertEqual(response.status_code, 303)

    def test_unauthenticated_api(self):
        self.assertEqual(self.client.post("/api/simular", json={}).status_code, 401)

    def test_unicode_login(self):
        self.login()
        self.assertEqual(self.client.get("/simulador").status_code, 200)

    def test_limits_reject_before_execution(self):
        self.login()
        for payload in ({"n_proteinas": 1000}, {"comprimento_proteina": -1},
                        {"composicao_fixa": {"ALA": 1000000}},
                        {"composicao_fixa": {"XXX": 5}},
                        {"virus_customizado": {"n_espiculas": 100000}},
                        {"ambiente": "unknown"}):
            with self.subTest(payload=payload):
                self.assertEqual(self.client.post("/api/simular", json=payload).status_code, 422)

    def test_concurrent_job_rejected(self):
        self.login()
        self.server.JOBS["busy"] = {"estado": "executando"}
        self.assertEqual(self.client.post("/api/simular", json={}).status_code, 429)

    def test_audit_attached_to_completed_job(self):
        self.server.JOBS["test"] = {"estado": "executando", "log": []}
        with tempfile.TemporaryDirectory() as tmp, patch.object(self.server, "OUTPUT_DIR", tmp), \
                patch.object(self.server.simulacao, "executar_simulacao", create=True,
                             return_value={"arquivos": {}, "ranking": []}):
            self.server._rodar_job("test", {})
            self.assertEqual(self.server.JOBS["test"]["estado"], "concluido")
            self.assertEqual(self.server.JOBS["test"]["resultado"]["arquivos"]["qualidade"], "qualidade.json")

    def test_exception_detail_not_exposed(self):
        self.server.JOBS["test"] = {"estado": "executando", "log": []}
        with tempfile.TemporaryDirectory() as tmp, patch.object(self.server, "OUTPUT_DIR", tmp), \
                patch.object(self.server.simulacao, "executar_simulacao", create=True,
                             side_effect=RuntimeError("private path")), self.assertLogs(level="ERROR"):
            self.server._rodar_job("test", {})
        self.assertNotIn("private path", self.server.JOBS["test"]["erro"])
        self.assertEqual(self.server.JOBS["test"]["estado"], "erro")

    def test_nonfinite_result_not_published(self):
        self.server.JOBS["test"] = {"estado": "executando", "log": []}
        with tempfile.TemporaryDirectory() as tmp, patch.object(self.server, "OUTPUT_DIR", tmp), \
                patch.object(self.server.simulacao, "executar_simulacao", create=True,
                             return_value={"ranking": [{"energia": float("nan")}] }), \
                self.assertLogs(level="ERROR"):
            self.server._rodar_job("test", {})
        self.assertEqual(self.server.JOBS["test"]["estado"], "erro")
        self.assertNotIn("resultado", self.server.JOBS["test"])

    def test_thread_start_failure_releases_slot(self):
        self.login()
        with patch.object(self.server, "threading") as threads, self.assertLogs(level="ERROR"):
            threads.Thread.return_value.start.side_effect = RuntimeError("busy")
            self.assertEqual(self.client.post("/api/simular", json={}).status_code, 503)
        self.assertEqual(self.server.JOBS, {})


if __name__ == "__main__":
    unittest.main()
