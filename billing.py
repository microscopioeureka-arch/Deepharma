"""Contas, cotas de simulação e integração segura com Stripe Checkout.

A cobrança fica DESATIVADA por padrão. O plano só é liberado após webhook
assinado do provedor; voltar da página de Checkout, sozinho, nunca concede acesso.
"""
from __future__ import annotations

import hashlib
import hmac
import json
import os
import re
import secrets
import sqlite3
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from typing import Any

LEGACY_EMAIL = "legacy-access@invalid.local"
ACTIVE_SUBSCRIPTION_STATUSES = {"active", "trialing"}
PASSWORD_ITERATIONS = 310_000
PLANS = {
    "free": {"label": "Gratuito", "limit": 5, "period": "mês"},
    # Limite sugerido para a versão de trabalho; confirmar antes de habilitar billing live.
    "deeplatinum": {"label": "Deeplatinum", "limit": 40, "period": "5 horas"},
    "deepharmax": {"label": "DeepharmaX", "limit": None, "period": "sem limite"},
}


class EmailAlreadyExists(Exception):
    pass


class InvalidCredentials(Exception):
    pass


class SimulationLimitReached(Exception):
    def __init__(self, summary: dict[str, Any]):
        self.summary = summary
        super().__init__("O limite de simulações deste plano foi atingido.")


class BillingConfigurationError(Exception):
    pass


class WebhookSignatureError(Exception):
    pass


def _db_path() -> str:
    return os.environ.get(
        "DEEPHARMA_DB_PATH",
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "deepharma.sqlite3"),
    )


def _connect() -> sqlite3.Connection:
    conn = sqlite3.connect(_db_path(), timeout=10)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON")
    conn.execute("PRAGMA busy_timeout = 10000")
    return conn


def init_db() -> None:
    """Cria tabelas locais sem substituir nem apagar dados existentes."""
    with _connect() as conn:
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS users (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                email TEXT NOT NULL UNIQUE,
                password_salt TEXT NOT NULL,
                password_hash TEXT NOT NULL,
                plan TEXT NOT NULL DEFAULT 'free',
                subscription_status TEXT NOT NULL DEFAULT 'inactive',
                stripe_customer_id TEXT UNIQUE,
                stripe_subscription_id TEXT UNIQUE,
                created_at INTEGER NOT NULL
            )
            """
        )
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS simulation_usage (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                user_id INTEGER NOT NULL REFERENCES users(id),
                job_id TEXT NOT NULL UNIQUE,
                created_at INTEGER NOT NULL,
                status TEXT NOT NULL DEFAULT 'started'
            )
            """
        )
        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_sim_usage_user_time "
            "ON simulation_usage(user_id, created_at, status)"
        )
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS stripe_webhook_events (
                event_id TEXT PRIMARY KEY,
                processed_at INTEGER NOT NULL
            )
            """
        )
        # Compatibilidade temporária com o formulário de senha compartilhada.
        # Essa conta não pode comprar plano nem representar uma assinatura individual.
        conn.execute(
            """
            INSERT OR IGNORE INTO users
                (email, password_salt, password_hash, plan, subscription_status, created_at)
            VALUES (?, 'legacy', 'legacy', 'free', 'inactive', ?)
            """,
            (LEGACY_EMAIL, int(time.time())),
        )


def _public_user(row: sqlite3.Row | dict | None) -> dict[str, Any] | None:
    if row is None:
        return None
    allowed = (
        "id", "email", "plan", "subscription_status",
        "stripe_customer_id", "stripe_subscription_id", "created_at",
    )
    return {key: row[key] for key in allowed if key in row.keys()}


def get_legacy_user() -> dict[str, Any]:
    init_db()
    with _connect() as conn:
        row = conn.execute("SELECT * FROM users WHERE email = ?", (LEGACY_EMAIL,)).fetchone()
    user = _public_user(row)
    if not user:
        raise RuntimeError("Não foi possível inicializar o acesso legado.")
    return user


def get_user(user_id: int | str) -> dict[str, Any] | None:
    try:
        normalized_id = int(user_id)
    except (TypeError, ValueError):
        return None
    with _connect() as conn:
        row = conn.execute("SELECT * FROM users WHERE id = ?", (normalized_id,)).fetchone()
    return _public_user(row)


def get_user_by_email(email: str) -> dict[str, Any] | None:
    normalized = (email or "").strip().lower()
    with _connect() as conn:
        row = conn.execute("SELECT * FROM users WHERE email = ?", (normalized,)).fetchone()
    return _public_user(row)


def _password_digest(password: str, salt_hex: str) -> str:
    return hashlib.pbkdf2_hmac(
        "sha256", password.encode("utf-8"), bytes.fromhex(salt_hex), PASSWORD_ITERATIONS
    ).hex()


def create_user(email: str, password: str) -> dict[str, Any]:
    normalized = (email or "").strip().lower()
    if not re.fullmatch(r"[^\s@]+@[^\s@]+\.[^\s@]+", normalized) or len(normalized) > 254:
        raise ValueError("Informe um e-mail válido.")
    if normalized == LEGACY_EMAIL:
        raise ValueError("Esse e-mail não pode ser usado.")
    if len(password) < 10:
        raise ValueError("A senha precisa ter pelo menos 10 caracteres.")
    if len(password) > 128:
        raise ValueError("A senha excede o tamanho permitido.")
    salt = secrets.token_hex(16)
    digest = _password_digest(password, salt)
    try:
        with _connect() as conn:
            cursor = conn.execute(
                """
                INSERT INTO users (email, password_salt, password_hash, plan, subscription_status, created_at)
                VALUES (?, ?, ?, 'free', 'inactive', ?)
                """,
                (normalized, salt, digest, int(time.time())),
            )
            user_id = cursor.lastrowid
    except sqlite3.IntegrityError as exc:
        raise EmailAlreadyExists("Já existe uma conta com esse e-mail.") from exc
    user = get_user(user_id)
    if not user:
        raise RuntimeError("Não foi possível recuperar a conta recém-criada.")
    return user


def verify_user(email: str, password: str) -> dict[str, Any] | None:
    normalized = (email or "").strip().lower()
    with _connect() as conn:
        row = conn.execute("SELECT * FROM users WHERE email = ?", (normalized,)).fetchone()
    if row is None or normalized == LEGACY_EMAIL:
        return None
    try:
        candidate = _password_digest(password, row["password_salt"])
    except (ValueError, TypeError):
        return None
    if not hmac.compare_digest(candidate, row["password_hash"]):
        return None
    return _public_user(row)


def _effective_plan(row: sqlite3.Row | dict) -> str:
    plan = row["plan"] if row["plan"] in PLANS else "free"
    if plan in {"deeplatinum", "deepharmax"}:
        if row["subscription_status"] not in ACTIVE_SUBSCRIPTION_STATUSES:
            return "free"
    return plan


def _month_start_epoch(now: int) -> int:
    dt = datetime.fromtimestamp(now, tz=timezone.utc)
    start = datetime(dt.year, dt.month, 1, tzinfo=timezone.utc)
    return int(start.timestamp())


def _window_for(plan: str, now: int) -> tuple[int, int | None, str]:
    if plan == "deeplatinum":
        window = max(60, int(os.environ.get("DEEPHARMA_PLATINUM_WINDOW_SECONDS", "18000")))
        limit = max(1, int(os.environ.get("DEEPHARMA_PLATINUM_LIMIT_5H", "40")))
        return now - window, limit, "últimas 5 horas"
    if plan == "deepharmax":
        return _month_start_epoch(now), None, "mês atual (sem limite)"
    limit = max(1, int(os.environ.get("DEEPHARMA_FREE_MONTHLY_LIMIT", "5")))
    return _month_start_epoch(now), limit, "mês atual"


def _usage_count(conn: sqlite3.Connection, user_id: int, since: int) -> int:
    row = conn.execute(
        """
        SELECT COUNT(*) AS total FROM simulation_usage
        WHERE user_id = ? AND created_at >= ? AND status IN ('started', 'succeeded')
        """,
        (user_id, since),
    ).fetchone()
    return int(row["total"])


def usage_summary(user_id: int | str) -> dict[str, Any]:
    user = get_user(user_id)
    if not user:
        raise ValueError("Conta não encontrada.")
    now = int(time.time())
    plan = _effective_plan(user)
    start, limit, period = _window_for(plan, now)
    with _connect() as conn:
        used = _usage_count(conn, int(user["id"]), start)
    return {
        "plan": plan,
        "plan_label": PLANS[plan]["label"],
        "subscription_status": user["subscription_status"],
        "used": used,
        "limit": limit,
        "remaining": None if limit is None else max(0, limit - used),
        "period": period,
        "billing_enabled": billing_is_configured(),
    }


def reserve_simulation(user_id: int | str, job_id: str) -> dict[str, Any]:
    """Reserva cota atomicamente para evitar duas requisições ultrapassarem o limite."""
    try:
        normalized_id = int(user_id)
    except (TypeError, ValueError) as exc:
        raise ValueError("É necessário entrar em uma conta para simular.") from exc
    now = int(time.time())
    with _connect() as conn:
        conn.execute("BEGIN IMMEDIATE")
        user = conn.execute("SELECT * FROM users WHERE id = ?", (normalized_id,)).fetchone()
        if user is None:
            raise ValueError("Conta não encontrada.")
        plan = _effective_plan(user)
        start, limit, period = _window_for(plan, now)
        used = _usage_count(conn, normalized_id, start)
        if limit is not None and used >= limit:
            raise SimulationLimitReached({
                "plan": plan, "plan_label": PLANS[plan]["label"],
                "used": used, "limit": limit, "remaining": 0, "period": period,
            })
        conn.execute(
            """
            INSERT INTO simulation_usage (user_id, job_id, created_at, status)
            VALUES (?, ?, ?, 'started')
            """,
            (normalized_id, job_id, now),
        )
        used += 1
    return {
        "plan": plan, "plan_label": PLANS[plan]["label"],
        "used": used, "limit": limit,
        "remaining": None if limit is None else max(0, limit - used),
        "period": period,
    }


def finish_simulation(job_id: str, status: str) -> None:
    if status not in {"succeeded", "failed"}:
        raise ValueError("Status de simulação inválido.")
    with _connect() as conn:
        conn.execute(
            "UPDATE simulation_usage SET status = ? WHERE job_id = ?",
            (status, job_id),
        )


def billing_is_configured() -> bool:
    enabled = os.environ.get("ENABLE_BILLING", "").strip().lower() in {"1", "true", "yes"}
    required = (
        "STRIPE_SECRET_KEY",
        "STRIPE_WEBHOOK_SECRET",
        "STRIPE_PRICE_DEEPLATINUM",
        "STRIPE_PRICE_DEEPHARMAX",
        "DEEPHARMA_PUBLIC_URL",
    )
    return enabled and all(os.environ.get(name, "").strip() for name in required)


def _stripe_request(path: str, fields: dict[str, str]) -> dict[str, Any]:
    secret = os.environ.get("STRIPE_SECRET_KEY", "").strip()
    if not secret:
        raise BillingConfigurationError("A integração de pagamentos ainda não foi configurada.")
    data = urllib.parse.urlencode(fields).encode("utf-8")
    req = urllib.request.Request(
        "https://api.stripe.com/v1/" + path.lstrip("/"),
        data=data,
        headers={
            "Authorization": "Bearer " + secret,
            "Content-Type": "application/x-www-form-urlencoded",
        },
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=15) as response:
            parsed = json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        details = exc.read().decode("utf-8", errors="replace")[:400]
        raise BillingConfigurationError(
            f"O provedor recusou a solicitação (HTTP {exc.code}). Confira a configuração."
        ) from exc
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
        raise BillingConfigurationError(
            "Não foi possível alcançar o provedor de pagamentos. Nenhuma assinatura foi ativada."
        ) from exc
    if not isinstance(parsed, dict) or parsed.get("error"):
        raise BillingConfigurationError("O provedor não retornou uma sessão válida.")
    return parsed


def create_checkout_session(user: dict[str, Any], plan: str) -> str:
    if not billing_is_configured():
        raise BillingConfigurationError("Os pagamentos ainda não estão habilitados.")
    if plan not in {"deeplatinum", "deepharmax"}:
        raise BillingConfigurationError("Plano inválido.")
    if user.get("email") == LEGACY_EMAIL:
        raise BillingConfigurationError("Crie uma conta individual antes de assinar.")
    if user.get("stripe_subscription_id") and user.get("subscription_status") in ACTIVE_SUBSCRIPTION_STATUSES:
        raise BillingConfigurationError(
            "Você já tem uma assinatura ativa. Use o portal para gerenciá-la."
        )
    price_id = os.environ.get("STRIPE_PRICE_" + plan.upper(), "").strip()
    if not price_id:
        raise BillingConfigurationError("O preço deste plano ainda não está configurado.")
    base = os.environ["DEEPHARMA_PUBLIC_URL"].rstrip("/")
    user_id = str(user["id"])
    fields = {
        "mode": "subscription",
        "success_url": base + "/billing/success?session_id={CHECKOUT_SESSION_ID}",
        "cancel_url": base + "/billing/cancel",
        "client_reference_id": user_id,
        "customer_email": user["email"],
        "metadata[user_id]": user_id,
        "metadata[plan]": plan,
        "subscription_data[metadata][user_id]": user_id,
        "subscription_data[metadata][plan]": plan,
        "line_items[0][price]": price_id,
        "line_items[0][quantity]": "1",
    }
    result = _stripe_request("checkout/sessions", fields)
    url = result.get("url")
    if not isinstance(url, str) or not url.startswith("https://checkout.stripe.com/"):
        raise BillingConfigurationError("O provedor não retornou um link seguro de Checkout.")
    return url


def create_customer_portal_session(user: dict[str, Any]) -> str:
    if not billing_is_configured():
        raise BillingConfigurationError("Os pagamentos ainda não estão habilitados.")
    customer_id = user.get("stripe_customer_id")
    if not customer_id:
        raise BillingConfigurationError("Não há assinatura para gerenciar nesta conta.")
    base = os.environ["DEEPHARMA_PUBLIC_URL"].rstrip("/")
    result = _stripe_request("billing_portal/sessions", {
        "customer": str(customer_id),
        "return_url": base + "/simulador",
    })
    url = result.get("url")
    if not isinstance(url, str) or not url.startswith("https://billing.stripe.com/"):
        raise BillingConfigurationError("O provedor não retornou um portal seguro.")
    return url


def verify_stripe_signature(payload: bytes, signature_header: str, now: int | None = None) -> None:
    secret = os.environ.get("STRIPE_WEBHOOK_SECRET", "").strip()
    if not secret:
        raise BillingConfigurationError("O webhook de pagamentos ainda não está configurado.")
    if not signature_header:
        raise WebhookSignatureError("Assinatura do webhook ausente.")
    parts: dict[str, list[str]] = {}
    for part in signature_header.split(","):
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        parts.setdefault(key.strip(), []).append(value.strip())
    timestamps = parts.get("t", [])
    signatures = parts.get("v1", [])
    if not timestamps or not signatures:
        raise WebhookSignatureError("Cabeçalho de assinatura inválido.")
    try:
        timestamp = int(timestamps[0])
    except ValueError as exc:
        raise WebhookSignatureError("Data da assinatura inválida.") from exc
    now = int(time.time()) if now is None else int(now)
    if abs(now - timestamp) > 300:
        raise WebhookSignatureError("Assinatura de webhook expirada.")
    signed_payload = str(timestamp).encode("ascii") + b"." + payload
    expected = hmac.new(secret.encode("utf-8"), signed_payload, hashlib.sha256).hexdigest()
    if not any(hmac.compare_digest(expected, candidate) for candidate in signatures):
        raise WebhookSignatureError("Assinatura de webhook inválida.")


def _price_to_plan(price_id: str | None) -> str | None:
    if not price_id:
        return None
    for plan in ("deeplatinum", "deepharmax"):
        configured = os.environ.get("STRIPE_PRICE_" + plan.upper(), "").strip()
        if configured and hmac.compare_digest(price_id, configured):
            return plan
    return None


def _extract_price_id(subscription: dict[str, Any]) -> str | None:
    try:
        price = subscription["items"]["data"][0]["price"]
    except (KeyError, IndexError, TypeError):
        return None
    return price.get("id") if isinstance(price, dict) else price if isinstance(price, str) else None


def process_stripe_webhook(payload: bytes, signature_header: str) -> bool:
    """Verifica assinatura e processa evento de forma idempotente. Retorna False se repetido."""
    verify_stripe_signature(payload, signature_header)
    try:
        event = json.loads(payload.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise WebhookSignatureError("O corpo do webhook não é JSON válido.") from exc
    event_id = event.get("id")
    event_type = event.get("type")
    obj = event.get("data", {}).get("object", {})
    if not isinstance(event_id, str) or not isinstance(event_type, str) or not isinstance(obj, dict):
        raise WebhookSignatureError("Formato de evento inválido.")

    now = int(time.time())
    with _connect() as conn:
        conn.execute("BEGIN IMMEDIATE")
        if conn.execute("SELECT 1 FROM stripe_webhook_events WHERE event_id = ?", (event_id,)).fetchone():
            return False

        if event_type == "checkout.session.completed":
            user_id = obj.get("client_reference_id") or obj.get("metadata", {}).get("user_id")
            customer = obj.get("customer")
            subscription = obj.get("subscription")
            if user_id:
                try:
                    normalized_id = int(user_id)
                except (TypeError, ValueError):
                    normalized_id = -1
                row = conn.execute("SELECT * FROM users WHERE id = ?", (normalized_id,)).fetchone()
                if row:
                    new_customer = customer if isinstance(customer, str) else None
                    new_subscription = subscription if isinstance(subscription, str) else None
                    # Não rebaixa uma assinatura se o evento de subscription chegou primeiro.
                    status = row["subscription_status"]
                    plan = row["plan"]
                    if status not in ACTIVE_SUBSCRIPTION_STATUSES:
                        status, plan = "pending", "free"
                    conn.execute(
                        """
                        UPDATE users SET stripe_customer_id = COALESCE(?, stripe_customer_id),
                            stripe_subscription_id = COALESCE(?, stripe_subscription_id),
                            subscription_status = ?, plan = ? WHERE id = ?
                        """,
                        (new_customer, new_subscription, status, plan, normalized_id),
                    )

        elif event_type in {"customer.subscription.created", "customer.subscription.updated"}:
            metadata = obj.get("metadata") or {}
            user_id = metadata.get("user_id")
            customer = obj.get("customer")
            subscription_id = obj.get("id")
            status = str(obj.get("status") or "inactive")
            price_plan = _price_to_plan(_extract_price_id(obj))
            plan = price_plan if status in ACTIVE_SUBSCRIPTION_STATUSES and price_plan else "free"
            if user_id:
                try:
                    normalized_id = int(user_id)
                except (TypeError, ValueError):
                    normalized_id = -1
                row = conn.execute("SELECT id FROM users WHERE id = ?", (normalized_id,)).fetchone()
            else:
                row = None
            if row is None and isinstance(customer, str):
                row = conn.execute("SELECT id FROM users WHERE stripe_customer_id = ?", (customer,)).fetchone()
            if row:
                normalized_id = int(row["id"])
                conn.execute(
                    """
                    UPDATE users SET plan = ?, subscription_status = ?,
                        stripe_customer_id = COALESCE(?, stripe_customer_id),
                        stripe_subscription_id = COALESCE(?, stripe_subscription_id)
                    WHERE id = ?
                    """,
                    (
                        plan, status,
                        customer if isinstance(customer, str) else None,
                        subscription_id if isinstance(subscription_id, str) else None,
                        normalized_id,
                    ),
                )

        elif event_type == "customer.subscription.deleted":
            customer = obj.get("customer")
            subscription_id = obj.get("id")
            row = None
            if isinstance(subscription_id, str):
                row = conn.execute(
                    "SELECT id FROM users WHERE stripe_subscription_id = ?", (subscription_id,)
                ).fetchone()
            if row is None and isinstance(customer, str):
                row = conn.execute("SELECT id FROM users WHERE stripe_customer_id = ?", (customer,)).fetchone()
            if row:
                conn.execute(
                    "UPDATE users SET plan = 'free', subscription_status = 'canceled', "
                    "stripe_subscription_id = NULL WHERE id = ?",
                    (int(row["id"]),),
                )

        elif event_type == "invoice.payment_failed":
            customer = obj.get("customer")
            if isinstance(customer, str):
                conn.execute(
                    "UPDATE users SET plan = 'free', subscription_status = 'past_due' "
                    "WHERE stripe_customer_id = ?",
                    (customer,),
                )

        conn.execute(
            "INSERT INTO stripe_webhook_events(event_id, processed_at) VALUES (?, ?)",
            (event_id, now),
        )
    return True
