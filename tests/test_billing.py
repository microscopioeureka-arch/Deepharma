import hashlib
import hmac
import os
import tempfile
import time
import unittest
from unittest import mock

import billing


_ENV_NAMES = (
    "DEEPHARMA_DB_PATH",
    "DEEPHARMA_DATABASE_URL",
    "DATABASE_URL",
    "DEEPHARMA_STORAGE_DURABLE",
    "ALLOW_LIVE_BILLING",
    "ENABLE_BILLING",
    "STRIPE_SECRET_KEY",
    "STRIPE_WEBHOOK_SECRET",
    "STRIPE_PRICE_DEEPLATINUM",
    "STRIPE_PRICE_DEEPHARMAX",
    "DEEPHARMA_PUBLIC_URL",
    "DEEPHARMA_FREE_MONTHLY_LIMIT",
    "DEEPHARMA_PLATINUM_LIMIT_5H",
    "DEEPHARMA_PLATINUM_WINDOW_SECONDS",
)


class BillingTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.original_env = {name: os.environ.get(name) for name in _ENV_NAMES}
        for name in _ENV_NAMES:
            os.environ.pop(name, None)
        os.environ["DEEPHARMA_DB_PATH"] = os.path.join(self.temp_dir.name, "test.sqlite3")
        billing.init_db()

    def tearDown(self):
        for name, value in self.original_env.items():
            if value is None:
                os.environ.pop(name, None)
            else:
                os.environ[name] = value
        self.temp_dir.cleanup()

    def test_passwords_are_verified_without_storing_plaintext(self):
        user = billing.create_user("  Researcher@Example.com ", "safe-test-pass-123")
        self.assertEqual(user["email"], "researcher@example.com")
        self.assertIsNotNone(billing.verify_user("researcher@example.com", "safe-test-pass-123"))
        self.assertIsNone(billing.verify_user("researcher@example.com", "incorrect-password"))
        with billing._connect() as conn:
            row = conn.execute(
                "SELECT password_salt, password_hash FROM users WHERE id = ?", (user["id"],)
            ).fetchone()
        self.assertNotEqual(row["password_hash"], "safe-test-pass-123")
        self.assertNotEqual(row["password_salt"], "safe-test-pass-123")

    def test_free_plan_enforces_five_started_simulations_per_month(self):
        user = billing.create_user("free@example.com", "safe-test-pass-123")
        for index in range(5):
            billing.reserve_simulation(user["id"], f"free-job-{index}")
        with self.assertRaises(billing.SimulationLimitReached):
            billing.reserve_simulation(user["id"], "free-job-over-limit")
        summary = billing.usage_summary(user["id"])
        self.assertEqual(summary["plan"], "free")
        self.assertEqual(summary["used"], 5)
        self.assertEqual(summary["remaining"], 0)

    def test_deeplatinum_uses_a_rolling_window_limit(self):
        user = billing.create_user("platinum@example.com", "safe-test-pass-123")
        with billing._connect() as conn:
            conn.execute(
                "UPDATE users SET plan = 'deeplatinum', subscription_status = 'active' WHERE id = ?",
                (user["id"],),
            )
        os.environ["DEEPHARMA_PLATINUM_LIMIT_5H"] = "2"
        os.environ["DEEPHARMA_PLATINUM_WINDOW_SECONDS"] = "18000"
        billing.reserve_simulation(user["id"], "platinum-job-1")
        billing.reserve_simulation(user["id"], "platinum-job-2")
        with self.assertRaises(billing.SimulationLimitReached):
            billing.reserve_simulation(user["id"], "platinum-job-over-limit")

    def test_deepharmax_has_no_quantity_limit(self):
        user = billing.create_user("max@example.com", "safe-test-pass-123")
        with billing._connect() as conn:
            conn.execute(
                "UPDATE users SET plan = 'deepharmax', subscription_status = 'active' WHERE id = ?",
                (user["id"],),
            )
        for index in range(8):
            billing.reserve_simulation(user["id"], f"max-job-{index}")
        summary = billing.usage_summary(user["id"])
        self.assertIsNone(summary["limit"])
        self.assertEqual(summary["used"], 8)

    def test_billing_is_off_when_required_configuration_is_missing(self):
        self.assertFalse(billing.billing_is_configured())
        self.assertFalse(billing.public_billing_status()["enabled"])

    def test_public_status_returns_validated_recurring_prices(self):
        os.environ.update({
            "ENABLE_BILLING": "true",
            "STRIPE_SECRET_KEY": "sk_test_placeholder",
            "STRIPE_WEBHOOK_SECRET": "whsec_test_placeholder",
            "STRIPE_PRICE_DEEPLATINUM": "price_platinum",
            "STRIPE_PRICE_DEEPHARMAX": "price_max",
            "DEEPHARMA_PUBLIC_URL": "https://example.test",
            "DEEPHARMA_DB_PATH": os.environ["DEEPHARMA_DB_PATH"],
            "DEEPHARMA_STORAGE_DURABLE": "true",
        })
        price = {
            "id": "price_platinum",
            "active": True,
            "currency": "usd",
            "unit_amount": 1200,
            "recurring": {"interval": "month", "interval_count": 1},
        }
        with mock.patch.object(billing, "_stripe_get", side_effect=[price, {**price, "id": "price_max", "unit_amount": 2500}]):
            result = billing.public_billing_status()
        self.assertTrue(result["enabled"])
        self.assertEqual(result["plans"]["deeplatinum"]["amount"], 1200)
        self.assertEqual(result["plans"]["deepharmax"]["amount"], 2500)

    def test_public_status_disables_a_non_monthly_price(self):
        os.environ.update({
            "ENABLE_BILLING": "true",
            "STRIPE_SECRET_KEY": "sk_test_placeholder",
            "STRIPE_WEBHOOK_SECRET": "whsec_test_placeholder",
            "STRIPE_PRICE_DEEPLATINUM": "price_platinum",
            "STRIPE_PRICE_DEEPHARMAX": "price_max",
            "DEEPHARMA_PUBLIC_URL": "https://example.test",
            "DEEPHARMA_DB_PATH": os.environ["DEEPHARMA_DB_PATH"],
        })
        monthly = {
            "id": "price_platinum",
            "active": True, "currency": "usd", "unit_amount": 1200,
            "recurring": {"interval": "month", "interval_count": 1},
        }
        yearly = {
            "id": "price_max",
            "active": True, "currency": "usd", "unit_amount": 2500,
            "recurring": {"interval": "year", "interval_count": 1},
        }
        with mock.patch.object(billing, "_stripe_get", side_effect=[monthly, yearly]):
            result = billing.public_billing_status()
        self.assertFalse(result["enabled"])

    def test_stripe_signature_validation_rejects_forgery_and_old_timestamps(self):
        secret = "whsec_unit_test"
        os.environ["STRIPE_WEBHOOK_SECRET"] = secret
        timestamp = int(time.time())
        payload = b'{"id":"evt_test","type":"example.event","data":{"object":{}}}'
        signed = str(timestamp).encode("ascii") + b"." + payload
        signature = hmac.new(secret.encode("utf-8"), signed, hashlib.sha256).hexdigest()
        header = f"t={timestamp},v1={signature}"

        billing.verify_stripe_signature(payload, header, now=timestamp)
        with self.assertRaises(billing.WebhookSignatureError):
            billing.verify_stripe_signature(payload, "t=" + str(timestamp) + ",v1=forged", now=timestamp)
        with self.assertRaises(billing.WebhookSignatureError):
            billing.verify_stripe_signature(payload, header, now=timestamp + 301)


    def test_postgres_url_can_enable_configured_test_prices_without_sqlite_path(self):
        os.environ.update({
            "ENABLE_BILLING": "true",
            "STRIPE_SECRET_KEY": "sk_test_placeholder",
            "STRIPE_WEBHOOK_SECRET": "whsec_test_placeholder",
            "STRIPE_PRICE_DEEPLATINUM": "price_platinum",
            "STRIPE_PRICE_DEEPHARMAX": "price_max",
            "DEEPHARMA_PUBLIC_URL": "https://example.test",
            "DATABASE_URL": "postgresql://user:pass@db.example.test/deepharma",
            "DEEPHARMA_DB_PATH": "",
            "DEEPHARMA_STORAGE_DURABLE": "true",
        })
        price = {
            "id": "price_platinum",
            "active": True,
            "currency": "brl",
            "unit_amount": 3900,
            "recurring": {"interval": "month", "interval_count": 1},
        }
        with mock.patch.object(
            billing, "_stripe_get",
            side_effect=[price, {**price, "id": "price_max", "unit_amount": 7900}],
        ):
            result = billing.public_billing_status()
        self.assertTrue(result["enabled"])
        self.assertEqual(result["plans"]["deeplatinum"]["amount"], 3900)

    def test_postgres_adapter_converts_parameter_placeholders(self):
        class FakeConnection:
            def __init__(self):
                self.query = None
                self.params = None

            def execute(self, query, params):
                self.query = query
                self.params = params
                return "cursor"

        fake = FakeConnection()
        adapter = billing._PostgresConnectionAdapter(fake)
        result = adapter.execute("SELECT ? WHERE user_id = ?", (1, 2))
        self.assertEqual(result, "cursor")
        self.assertEqual(fake.query, "SELECT %s WHERE user_id = %s")
        self.assertEqual(fake.params, (1, 2))

    def test_postgres_schema_uses_portable_identity_columns(self):
        statements = []

        class FakeCursor:
            rowcount = 1

        class FakeConnection:
            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, traceback):
                return False

            def execute(self, query, params=()):
                statements.append(query)
                return FakeCursor()

        with mock.patch.object(billing, "_using_postgres", return_value=True), \
             mock.patch.object(billing, "_connect", return_value=FakeConnection()):
            billing.init_db()

        combined = "\\n".join(statements)
        self.assertIn("BIGSERIAL PRIMARY KEY", combined)
        self.assertIn("user_id BIGINT NOT NULL", combined)
        self.assertNotIn("AUTOINCREMENT", combined)
        self.assertIn("ON CONFLICT (email) DO NOTHING", combined)

if __name__ == "__main__":
    unittest.main()
