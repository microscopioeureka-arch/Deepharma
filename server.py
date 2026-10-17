# -*- coding: utf-8 -*-
"""
server.py — servidor web (FastAPI). Versão "pasta única": este arquivo,
simulacao.py, index.html, style.css e app.js ficam todos na MESMA pasta.

Rode com:
    uvicorn server:app --reload --port 8000

Depois abra http://localhost:8000 no navegador.
"""

import os
import uuid
import threading
import secrets
import logging
from typing import Literal
from scientific_quality import build_quality_report, write_quality_report

from fastapi import FastAPI, HTTPException, Request, Form
from fastapi.staticfiles import StaticFiles
from fastapi.responses import HTMLResponse, RedirectResponse, JSONResponse, PlainTextResponse, FileResponse
from starlette.middleware.sessions import SessionMiddleware
from pydantic import BaseModel, Field, ConfigDict, model_validator

import simulacao  # mesmo diretório, import direto
import billing

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")   # criada em tempo de execução
WEB_DIR = os.path.join(BASE_DIR, "web")
os.makedirs(OUTPUT_DIR, exist_ok=True)

app = FastAPI(title="Simulador de Proteínas Híbridas")

# A senha e as chaves ficam nas variáveis secretas do servidor, nunca no JavaScript.
APP_PASSWORD = os.environ.get("APP_PASSWORD")
SESSION_SECRET = os.environ.get("SESSION_SECRET")
if not APP_PASSWORD or not SESSION_SECRET:
    raise RuntimeError(
        "Defina APP_PASSWORD e SESSION_SECRET antes de iniciar o servidor."
    )

billing.init_db()

LOGIN_HTML = """
<!doctype html><html lang="pt-BR"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Entrar — Deepharma</title><style>
*{box-sizing:border-box}body{margin:0;min-height:100vh;display:grid;place-items:center;padding:24px;background:#101714;color:#eef5f0;font-family:system-ui,sans-serif}
.auth{width:min(100%,430px);padding:28px;border:1px solid #285046;border-radius:18px;background:#12211c}
h1{margin-top:0}p{color:#b2c5bb;line-height:1.5}label{display:block;margin:14px 0 6px}
input{width:100%;padding:12px;border:1px solid #466b5b;border-radius:8px;background:#0d1814;color:white;font-size:1rem}
button{width:100%;margin-top:16px;padding:13px;border:0;border-radius:8px;background:#40d7a0;color:#062117;font-weight:800;cursor:pointer}
.erro{color:#ffaaaa;min-height:1.2em}a{color:#73e6ba}.hint{font-size:.85rem;color:#91a99b}
</style></head><body><form class="auth" method="post" action="/login">
<h1>Deepharma</h1><p>Entre na sua conta para acessar o simulador.</p>
<input type="hidden" name="csrf_token" value="__CSRF__">
<label for="email">E-mail da conta (deixe em branco para o acesso legado)</label>
<input id="email" name="email" type="email" autocomplete="username">
<label for="senha">Senha</label><input id="senha" name="senha" type="password" autocomplete="current-password" required>
<button type="submit">Entrar</button><div class="erro">__ERRO__</div>
<p>Primeira vez aqui? <a href="/cadastro">Criar conta</a></p>
<p class="hint">O acesso legado por senha compartilhada não pode contratar planos pagos.</p>
</form></body></html>
"""

CADASTRO_HTML = """
<!doctype html><html lang="pt-BR"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Criar conta — Deepharma</title><style>
*{box-sizing:border-box}body{margin:0;min-height:100vh;display:grid;place-items:center;padding:24px;background:#101714;color:#eef5f0;font-family:system-ui,sans-serif}
.auth{width:min(100%,430px);padding:28px;border:1px solid #285046;border-radius:18px;background:#12211c}
h1{margin-top:0}p{color:#b2c5bb;line-height:1.5}label{display:block;margin:14px 0 6px}
input{width:100%;padding:12px;border:1px solid #466b5b;border-radius:8px;background:#0d1814;color:white;font-size:1rem}
button{width:100%;margin-top:16px;padding:13px;border:0;border-radius:8px;background:#40d7a0;color:#062117;font-weight:800;cursor:pointer}
.erro{color:#ffaaaa;min-height:1.2em}a{color:#73e6ba}.hint{font-size:.85rem;color:#91a99b}
</style></head><body><form class="auth" method="post" action="/cadastro">
<h1>Criar conta</h1><p>Conta individual para guardar seu acesso e limite de simulações.</p>
<input type="hidden" name="csrf_token" value="__CSRF__">
<label for="email">E-mail</label><input id="email" name="email" type="email" autocomplete="email" required maxlength="254">
<label for="senha">Senha (mínimo 10 caracteres)</label><input id="senha" name="senha" type="password" autocomplete="new-password" required minlength="10" maxlength="128">
<label for="confirmar_senha">Confirmar senha</label><input id="confirmar_senha" name="confirmar_senha" type="password" autocomplete="new-password" required minlength="10" maxlength="128">
<button type="submit">Criar conta gratuita</button><div class="erro">__ERRO__</div>
<p>Já tem conta? <a href="/login">Entrar</a></p><p class="hint">Antes da abertura pública, faltam verificação de e-mail e banco persistente.</p>
</form></body></html>
"""


def pagina_login(erro="", csrf_token=""):
    return LOGIN_HTML.replace("__ERRO__", erro).replace("__CSRF__", csrf_token)


def pagina_cadastro(erro="", csrf_token=""):
    return CADASTRO_HTML.replace("__ERRO__", erro).replace("__CSRF__", csrf_token)


def _csrf_token(request: Request) -> str:
    token = request.session.get("csrf_token")
    if not token:
        token = secrets.token_urlsafe(32)
        request.session["csrf_token"] = token
    return token


def _csrf_valido(request: Request, token: str) -> bool:
    saved = request.session.get("csrf_token", "")
    return bool(saved and token and secrets.compare_digest(str(saved), str(token)))


def _iniciar_sessao_usuario(request: Request, user: dict) -> None:
    request.session.clear()
    request.session["logged_in"] = True
    request.session["user_id"] = int(user["id"])
    request.session["email"] = str(user["email"])
    request.session["csrf_token"] = secrets.token_urlsafe(32)


PUBLIC_PATHS = {
    "/", "/landing.css", "/landing.js", "/login", "/cadastro", "/register",
    "/healthz", "/favicon.ico", "/api/csrf", "/api/billing/status",
    "/api/billing/webhook", "/billing/success", "/billing/cancel",
}


@app.middleware("http")
async def exigir_login(request: Request, call_next):
    caminho = request.url.path
    if caminho not in PUBLIC_PATHS and not request.session.get("logged_in"):
        if caminho.startswith("/api/"):
            return JSONResponse({"detail": "Faça login para usar o Deepharma."}, status_code=401)
        return RedirectResponse("/login", status_code=303)
    return await call_next(request)


COOKIE_SECURE = os.environ.get("COOKIE_SECURE", "0").lower() in {"1", "true", "yes"}
app.add_middleware(SessionMiddleware, secret_key=SESSION_SECRET, same_site="lax",
                   https_only=COOKIE_SECURE, max_age=60 * 60 * 24 * 30)


@app.get("/login", response_class=HTMLResponse)
async def tela_login(request: Request, checkout: str = ""):
    if checkout in {"deeplatinum", "deepharmax"}:
        request.session["checkout_after_login"] = checkout
    if request.session.get("logged_in"):
        pending = request.session.get("checkout_after_login")
        if pending in {"deeplatinum", "deepharmax"}:
            return RedirectResponse("/?checkout=" + pending, status_code=303)
        return RedirectResponse("/simulador", status_code=303)
    return HTMLResponse(pagina_login("", _csrf_token(request)))


@app.post("/login", response_class=HTMLResponse)
async def fazer_login(request: Request, senha: str = Form(...), email: str = Form(""),
                      csrf_token: str = Form("")):
    if not _csrf_valido(request, csrf_token):
        return HTMLResponse(pagina_login("Sessão expirada. Recarregue a página.", _csrf_token(request)), status_code=403)
    if email.strip():
        user = billing.verify_user(email, senha)
        if not user:
            return HTMLResponse(pagina_login("E-mail ou senha incorretos.", _csrf_token(request)), status_code=401)
    elif secrets.compare_digest(senha.encode("utf-8"), APP_PASSWORD.encode("utf-8")):
        user = billing.get_legacy_user()
    else:
        return HTMLResponse(pagina_login("Credenciais inválidas.", _csrf_token(request)), status_code=401)
    pending = request.session.get("checkout_after_login")
    _iniciar_sessao_usuario(request, user)
    if pending in {"deeplatinum", "deepharmax"}:
        return RedirectResponse("/?checkout=" + pending, status_code=303)
    return RedirectResponse("/simulador", status_code=303)


@app.get("/cadastro", response_class=HTMLResponse)
async def tela_cadastro(request: Request, checkout: str = ""):
    if checkout in {"deeplatinum", "deepharmax"}:
        request.session["checkout_after_login"] = checkout
    if request.session.get("logged_in"):
        pending = request.session.get("checkout_after_login")
        if pending in {"deeplatinum", "deepharmax"}:
            return RedirectResponse("/?checkout=" + pending, status_code=303)
        return RedirectResponse("/simulador", status_code=303)
    return HTMLResponse(pagina_cadastro("", _csrf_token(request)))


@app.get("/register")
async def cadastro_compatibilidade():
    return RedirectResponse("/cadastro", status_code=303)


@app.post("/cadastro", response_class=HTMLResponse)
async def criar_conta(request: Request, email: str = Form(...), senha: str = Form(...),
                      confirmar_senha: str = Form(...), csrf_token: str = Form("")):
    if not _csrf_valido(request, csrf_token):
        return HTMLResponse(pagina_cadastro("Sessão expirada. Recarregue a página.", _csrf_token(request)), status_code=403)
    if senha != confirmar_senha:
        return HTMLResponse(pagina_cadastro("As senhas não coincidem.", _csrf_token(request)), status_code=400)
    try:
        user = billing.create_user(email, senha)
    except billing.EmailAlreadyExists as exc:
        return HTMLResponse(pagina_cadastro(str(exc), _csrf_token(request)), status_code=409)
    except ValueError as exc:
        return HTMLResponse(pagina_cadastro(str(exc), _csrf_token(request)), status_code=400)
    pending = request.session.get("checkout_after_login")
    _iniciar_sessao_usuario(request, user)
    if pending in {"deeplatinum", "deepharmax"}:
        return RedirectResponse("/?checkout=" + pending, status_code=303)
    return RedirectResponse("/simulador", status_code=303)


@app.get("/logout")
async def fazer_logout(request: Request):
    request.session.clear()
    return RedirectResponse("/login", status_code=303)


@app.get("/healthz")
async def healthz():
    return PlainTextResponse("ok")


@app.get("/api/csrf")
async def csrf_token(request: Request):
    return {"csrf_token": _csrf_token(request)}


@app.get("/api/billing/status")
def billing_status():
    return billing.public_billing_status()


@app.get("/api/usage")
def uso_mensal(request: Request):
    user_id = request.session.get("user_id")
    if not user_id:
        raise HTTPException(status_code=401, detail="Entre em uma conta individual.")
    try:
        return billing.usage_summary(user_id)
    except ValueError as exc:
        raise HTTPException(status_code=401, detail=str(exc))


class CheckoutRequest(BaseModel):
    plan: Literal["deeplatinum", "deepharmax"]


@app.post("/api/billing/checkout")
def iniciar_checkout(payload: CheckoutRequest, request: Request):
    user_id = request.session.get("user_id")
    if not user_id:
        raise HTTPException(status_code=401, detail="Entre em uma conta individual para assinar.")
    if not _csrf_valido(request, request.headers.get("x-csrf-token", "")):
        raise HTTPException(status_code=403, detail="Token de segurança expirado. Atualize a página.")
    user = billing.get_user(user_id)
    if not user:
        raise HTTPException(status_code=401, detail="Conta não encontrada.")
    try:
        url = billing.create_checkout_session(user, payload.plan)
    except billing.BillingConfigurationError as exc:
        status_code = 503 if not billing.billing_is_configured() else 400
        raise HTTPException(status_code=status_code, detail=str(exc))
    return {"url": url}


@app.post("/api/billing/portal")
def portal_assinatura(request: Request):
    user_id = request.session.get("user_id")
    if not user_id:
        raise HTTPException(status_code=401, detail="Entre em uma conta individual.")
    if not _csrf_valido(request, request.headers.get("x-csrf-token", "")):
        raise HTTPException(status_code=403, detail="Token de segurança expirado. Atualize a página.")
    user = billing.get_user(user_id)
    if not user:
        raise HTTPException(status_code=401, detail="Conta não encontrada.")
    try:
        return {"url": billing.create_customer_portal_session(user)}
    except billing.BillingConfigurationError as exc:
        status_code = 503 if not billing.billing_is_configured() else 400
        raise HTTPException(status_code=status_code, detail=str(exc))


@app.post("/api/billing/webhook")
async def stripe_webhook(request: Request):
    if not os.environ.get("STRIPE_WEBHOOK_SECRET", "").strip():
        raise HTTPException(status_code=503, detail="Webhook de pagamentos não configurado.")
    payload = await request.body()
    try:
        processed = billing.process_stripe_webhook(payload, request.headers.get("stripe-signature", ""))
    except billing.WebhookSignatureError as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    except billing.BillingConfigurationError as exc:
        raise HTTPException(status_code=503, detail=str(exc))
    return {"received": True, "duplicate": not processed}


@app.get("/billing/success", response_class=HTMLResponse)
def billing_success():
    return HTMLResponse("<!doctype html><html lang='pt-BR'><meta charset='utf-8'><title>Assinatura em verificação</title>"
        "<body style='font-family:system-ui;background:#101714;color:#eaf5ee;padding:3rem'><h1>Checkout concluído</h1>"
        "<p>A assinatura só será ativada após confirmação segura do provedor.</p><a href='/simulador'>Voltar ao simulador</a></body></html>")


@app.get("/billing/cancel", response_class=HTMLResponse)
def billing_cancel():
    return HTMLResponse("<!doctype html><html lang='pt-BR'><meta charset='utf-8'><title>Checkout cancelado</title>"
        "<body style='font-family:system-ui;background:#101714;color:#eaf5ee;padding:3rem'><h1>Checkout cancelado</h1>"
        "<p>Nenhum plano foi ativado por esta visita.</p><a href='/'>Voltar ao Deepharma</a></body></html>")


# Output files are served through an account-ownership-checked route below.

JOBS: dict[str, dict] = {}
_lock = threading.Lock()


@app.get("/outputs/{job_id}/{file_path:path}")
def arquivo_resultado(job_id: str, file_path: str, request: Request):
    with _lock:
        job = JOBS.get(job_id)
    if job is None or job.get("user_id") != request.session.get("user_id"):
        raise HTTPException(status_code=404, detail="arquivo não encontrado")

    job_dir = os.path.realpath(os.path.join(OUTPUT_DIR, job_id))
    caminho = os.path.realpath(os.path.join(job_dir, file_path))
    if os.path.commonpath([job_dir, caminho]) != job_dir or not os.path.isfile(caminho):
        raise HTTPException(status_code=404, detail="arquivo não encontrado")
    return FileResponse(caminho)


class VirusCustomizado(BaseModel):
    model_config = ConfigDict(allow_inf_nan=False)
    nome: str = Field(default="", max_length=200)
    preset: str = ""                  # "covid19" | "hiv" | "influenza" | "" (vazio = customizado)
    tipo_acido: Literal["RNA", "DNA"] = "RNA"           # "RNA" ou "DNA" (só usado se preset vazio)
    sequencia_acido: str = Field(default="", max_length=10000)          # texto cru digitado pelo usuário (só se preset vazio)
    raio_envelope: float = Field(default=6.0, ge=1, le=20)
    n_espiculas: int = Field(default=30, ge=1, le=100)
    usar_imagem_referencia: bool = False   # ajusta cor/alongamento com foto real (Wikimedia)
    nome_organismo_imagem: str = ""        # termo de busca da imagem (ex.: organismo do resultado escolhido)


class BacteriaGenomaCustomizado(BaseModel):
    tipo_acido: str = "DNA"          # "DNA" ou "RNA"
    sequencia_acido: str = ""         # texto cru (digitado ou vindo da busca no NCBI)
    usar_imagem_referencia: bool = False   # ajusta cor/proporção da 1ª bactéria com foto real
    nome_organismo_imagem: str = ""        # termo de busca da imagem (ex.: organismo do resultado escolhido)


class SimularRequest(BaseModel):
    model_config = ConfigDict(allow_inf_nan=False)
    # ---- proteína nova por combinação de aminoácidos ------------------
    comprimento_proteina: int = Field(default=30, ge=5, le=80)       # nº total de resíduos da cadeia
    aminoacidos_pool: list[str] = Field(default_factory=list)     # AAs permitidos ([] = todos os 20)
    composicao_fixa: dict[str, int] = Field(default_factory=dict) # {"ALA":4,"LEU":6} — ignora comprimento/pool
    n_proteinas: int = Field(default=1, ge=1, le=5)                 # quantas proteínas gerar (1-5)

    # ---- ambiente e velocidade ----------------------------------------
    ambiente: str = "citoplasma"
    modo_rapido: bool = True
    limite_energetico: float | None = None  # energia máxima aceita na formação; menor é melhor
    tipo_resposta: str = "nenhuma"          # estratégia didática contra vírus/bactérias

    # ---- compostos reais como ligantes na triagem --------------------
    elementos: list[str] = Field(default_factory=lambda: ["Paracetamol", "Cafeina", "Cloroquina"])
    incluir_compostos_reais: bool = False

    # ---- alvos -------------------------------------------------------
    gerar_virus: bool = True
    gerar_bacterias: bool = False
    n_bacterias: int = 2                 # 1-6 — cada bactéria gerada é diferente das outras
    virus_customizado: VirusCustomizado | None = None
    bacteria_genoma_customizado: BacteriaGenomaCustomizado | None = None

    # ---- previsão de mutações nos vírus -------------------------------
    prever_mutacao: bool = False
    n_mutacoes: int = Field(default=3, ge=1, le=10)


    @model_validator(mode="after")
    def validar_configuracao(self):
        if self.ambiente not in simulacao.ENVIRONMENTS:
            raise ValueError("Ambiente desconhecido.")
        if any(aa not in simulacao.AMINOACIDOS for aa in self.aminoacidos_pool):
            raise ValueError("Aminoácido desconhecido.")
        if self.composicao_fixa:
            if any(aa not in simulacao.AMINOACIDOS or q < 0
                   for aa, q in self.composicao_fixa.items()):
                raise ValueError("Composição inválida: use aminoácidos conhecidos e contagens não negativas.")
            if not 5 <= sum(self.composicao_fixa.values()) <= 80:
                raise ValueError("A composição deve totalizar de 5 a 80 resíduos.")
        if len(self.elementos) > len(simulacao.COMPOSTOS_REAIS) or any(
                e not in simulacao.COMPOSTOS_REAIS for e in self.elementos):
            raise ValueError("Composto desconhecido ou lista excessiva.")
        return self


@app.get("/api/opcoes")
def opcoes():
    return {
        "compostos": list(simulacao.COMPOSTOS_REAIS.keys()),
        "ambientes": [
            {"chave": k, "label": v["label"]} for k, v in simulacao.ENVIRONMENTS.items()
        ],
        "aminoacidos": list(simulacao.AMINOACIDOS.keys()),
        "virus_presets": [
            {"chave": k, "nome": v["nome"], "fonte": v["fonte"]}
            for k, v in simulacao.VIRUS_REAIS_PRESETS.items()
        ],
        "tipos_resposta": [
            {"chave": k, "nome": v["label"], "descricao": v["descricao"]}
            for k, v in simulacao.TIPOS_RESPOSTA.items()
        ],
    }


def _rodar_job(job_id: str, config: dict):
    job_dir = os.path.join(OUTPUT_DIR, job_id)
    os.makedirs(job_dir, exist_ok=True)

    def progresso(msg):
        with _lock:
            JOBS[job_id]["log"].append(msg)

    try:
        resultado = simulacao.executar_simulacao(config, job_dir, progresso=progresso)
        quality = build_quality_report(config, resultado, BASE_DIR)
        write_quality_report(quality, os.path.join(job_dir, "qualidade.json"))
        resultado["qualidade_cientifica"] = quality
        resultado.setdefault("arquivos", {})["qualidade"] = "qualidade.json"
        if quality["integridade_numerica"] != "aprovada":
            raise ValueError("Resultado contém valores não finitos; consulte o registro do servidor.")
        try:
            billing.finish_simulation(job_id, "succeeded")
        except Exception:
            logging.exception("Não foi possível finalizar registro de uso do job %s", job_id)
        with _lock:
            JOBS[job_id]["estado"] = "concluido"
            JOBS[job_id]["resultado"] = resultado
    except Exception as e:
        try:
            billing.finish_simulation(job_id, "failed")
        except Exception:
            logging.exception("Não foi possível registrar falha de uso do job %s", job_id)
        with _lock:
            JOBS[job_id]["estado"] = "erro"
            JOBS[job_id]["erro"] = "A execução falhou. Informe o identificador do job ao administrador."
        logging.exception("Falha no job %s", job_id)


@app.get("/api/buscar/acido_nucleico")
def buscar_acido_nucleico(termo: str, max_bases: int = 90, tipo: str = "RNA"):
    """Busca uma sequência de DNA/RNA REAL no NCBI (qualquer organismo —
    vírus, bactéria, gene específico, ou accession) para usar como
    genoma customizado de um vírus ou nucleoide/plasmídeo de bactéria.
    `tipo=RNA` (padrão, usado para vírus) converte T->U; `tipo=DNA`
    (usado para bactérias) mantém as bases como estão depositadas.
    Atalho de compatibilidade: baixa direto o 1º resultado. Prefira o
    par /lista + /detalhe abaixo, que deixa o usuário escolher."""
    try:
        return simulacao.buscar_acido_nucleico_ncbi(
            termo, max_bases=max_bases, preferir_rna=(tipo.upper() != "DNA"))
    except RuntimeError as e:
        raise HTTPException(status_code=404, detail=str(e))


@app.get("/api/buscar/acido_nucleico/lista")
def buscar_acido_nucleico_lista(termo: str, max_resultados: int = 8):
    """Lista candidatos do NCBI Nucleotide (título, organismo, tamanho)
    SEM baixar a sequência inteira ainda — para o usuário escolher qual
    accession usar antes de baixar (ver /detalhe)."""
    try:
        return {"resultados": simulacao.buscar_acido_nucleico_ncbi_lista(termo, max_resultados=max_resultados)}
    except RuntimeError as e:
        raise HTTPException(status_code=404, detail=str(e))


@app.get("/api/buscar/acido_nucleico/detalhe")
def buscar_acido_nucleico_detalhe(accession: str, max_bases: int = 90, tipo: str = "RNA"):
    """Baixa a sequência completa de UM accession específico, já
    escolhido pelo usuário na lista de /lista acima."""
    try:
        return simulacao.obter_acido_nucleico_por_accession(
            accession, max_bases=max_bases, preferir_rna=(tipo.upper() != "DNA"))
    except RuntimeError as e:
        raise HTTPException(status_code=404, detail=str(e))


@app.get("/api/buscar/proteina")
def buscar_proteina(termo: str, max_residuos: int = 40):
    """Busca uma proteína REAL no UniProt (gene, proteína ou organismo)
    para usar como composição exata na geração da proteína.
    Atalho de compatibilidade: baixa direto o 1º resultado. Prefira o
    par /lista + /detalhe abaixo, que deixa o usuário escolher."""
    try:
        return simulacao.buscar_proteina_uniprot(termo, max_residuos=max_residuos)
    except RuntimeError as e:
        raise HTTPException(status_code=404, detail=str(e))


@app.get("/api/buscar/proteina/lista")
def buscar_proteina_lista(termo: str, max_resultados: int = 8):
    """Lista candidatos do UniProt (nome, organismo, tamanho) SEM baixar
    a sequência completa ainda — para o usuário escolher qual accession
    usar antes de baixar (ver /detalhe)."""
    try:
        return {"resultados": simulacao.buscar_proteina_uniprot_lista(termo, max_resultados=max_resultados)}
    except RuntimeError as e:
        raise HTTPException(status_code=404, detail=str(e))


@app.get("/api/buscar/proteina/detalhe")
def buscar_proteina_detalhe(accession: str, max_residuos: int = 40):
    """Baixa a sequência completa de UMA proteína específica, já
    escolhida pelo usuário na lista de /lista acima."""
    try:
        return simulacao.obter_proteina_por_accession(accession, max_residuos=max_residuos)
    except RuntimeError as e:
        raise HTTPException(status_code=404, detail=str(e))


@app.post("/api/simular")
def simular(req: SimularRequest, request: Request):
    user_id = request.session.get("user_id")
    if not user_id:
        raise HTTPException(status_code=401, detail="Entre em uma conta para usar o simulador.")

    job_id = uuid.uuid4().hex
    config = req.model_dump() if hasattr(req, "model_dump") else req.dict()
    # Mantém a configuração antiga do frontend compatível com o backend.
    config["gerar_bacterias"] = False
    config["n_bacterias"] = 0
    config["bacteria_genoma_customizado"] = None

    with _lock:
        if any(job["estado"] == "executando" for job in JOBS.values()):
            raise HTTPException(status_code=429, detail="Há uma simulação em execução. Aguarde sua conclusão.")
        try:
            quota = billing.reserve_simulation(user_id, job_id)
        except billing.SimulationLimitReached as exc:
            raise HTTPException(
                status_code=429,
                detail={"message": str(exc), "usage": exc.summary},
            )
        except ValueError as exc:
            raise HTTPException(status_code=401, detail=str(exc))

        while len(JOBS) >= 100:
            del JOBS[next(iter(JOBS))]
        JOBS[job_id] = {
            "estado": "executando", "log": [], "resultado": None, "erro": None,
            "user_id": int(user_id), "quota": quota,
        }

    t = threading.Thread(target=_rodar_job, args=(job_id, config), daemon=True)
    try:
        t.start()
    except Exception:
        try:
            billing.finish_simulation(job_id, "failed")
        except Exception:
            logging.exception("Falha ao finalizar reserva de cota do job %s", job_id)
        with _lock:
            JOBS.pop(job_id, None)
        logging.exception("Falha ao iniciar job %s", job_id)
        raise HTTPException(status_code=503, detail="Não foi possível iniciar a execução.")
    return {"job_id": job_id}

@app.get("/api/status/{job_id}")
def status(job_id: str, request: Request):
    with _lock:
        job = JOBS.get(job_id)
    if job is None or job.get("user_id") != request.session.get("user_id"):
        raise HTTPException(status_code=404, detail="job não encontrado")
    resultado = job.get("resultado")
    resultado_web = _prefixar_urls(resultado, job_id) if resultado else None
    return {
        "estado": job["estado"], "log": job["log"],
        "resultado": resultado_web, "erro": job.get("erro"),
    }


def _prefixar_urls(resultado: dict, job_id: str) -> dict:
    def full(rel):
        return f"/outputs/{job_id}/{rel}" if rel else None

    out = dict(resultado)
    out["arquivos"] = {k: full(v) for k, v in resultado.get("arquivos", {}).items()}
    out["proteinas"] = [{**p, "glb": full(p.get("glb"))} for p in resultado.get("proteinas", [])]
    out["virus"] = [{**v, "glb": full(v.get("glb"))} for v in resultado.get("virus", [])]
    out["bacterias"] = []  # removido da interface e da execução pública
    out["docking"] = [{**d, "glb": full(d.get("glb"))} for d in resultado.get("docking", [])]
    out["mutacoes"] = [{**m, "glb": full(m.get("glb"))} for m in resultado.get("mutacoes", [])]
    mt = resultado.get("mutacao_top")
    out["mutacao_top"] = {**mt, "glb": full(mt.get("glb"))} if mt else None
    return out





# Serve o frontend (index.html, style.css, app.js) na raiz "/" —
# como tudo está na mesma pasta, servimos a própria BASE_DIR.
# O StaticFiles precisa ser montado por ÚLTIMO, senão ele "engole" as
# rotas /api/... definidas acima.
# Homepage pública + simulador autenticado.
@app.get("/", response_class=HTMLResponse)
def homepage():
    return FileResponse(os.path.join(WEB_DIR, "landing.html"), media_type="text/html")


@app.get("/simulador", response_class=HTMLResponse)
def pagina_simulador():
    return FileResponse(os.path.join(WEB_DIR, "index.html"), media_type="text/html")


@app.get("/style.css")
def frontend_css():
    return FileResponse(os.path.join(WEB_DIR, "style.css"), media_type="text/css")


@app.get("/design.css")
def frontend_design_css():
    return FileResponse(os.path.join(WEB_DIR, "design.css"), media_type="text/css")


@app.get("/app.js")
def frontend_js():
    return FileResponse(os.path.join(WEB_DIR, "app.js"), media_type="application/javascript")


@app.get("/landing.css")
def homepage_css():
    return FileResponse(os.path.join(WEB_DIR, "landing.css"), media_type="text/css")


@app.get("/landing.js")
def homepage_js():
    return FileResponse(os.path.join(WEB_DIR, "landing.js"), media_type="application/javascript")
