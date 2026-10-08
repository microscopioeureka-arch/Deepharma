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
import traceback
import secrets

from fastapi import FastAPI, HTTPException, Request, Form
from fastapi.staticfiles import StaticFiles
from fastapi.responses import HTMLResponse, RedirectResponse, JSONResponse, PlainTextResponse
from starlette.middleware.sessions import SessionMiddleware
from pydantic import BaseModel, Field

import simulacao  # mesmo diretório, import direto

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

LOGIN_HTML = """
<!doctype html>
<html lang="pt-BR">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width,initial-scale=1">
  <title>Acesso restrito — Simulador Molecular</title>
  <style>
    * { box-sizing: border-box; }
    body { margin: 0; min-height: 100vh; display: grid; place-items: center;
           background: #0d1117; color: #eef2f7; font-family: Arial, sans-serif; }
    .login { width: min(92vw, 390px); padding: 30px; border-radius: 18px;
             background: #1b2330; box-shadow: 0 18px 60px #0009; }
    h1 { margin: 0 0 10px; font-size: 1.45rem; }
    p { color: #aeb9c8; line-height: 1.5; }
    label { display: block; margin: 18px 0 7px; color: #c8d0db; }
    input { width: 100%; padding: 13px; border: 1px solid #46566d;
            border-radius: 9px; background: #101721; color: white; font-size: 1rem; }
    button { width: 100%; margin-top: 18px; padding: 13px; border: 0;
             border-radius: 9px; background: #50c878; color: #07130b;
             font-weight: bold; cursor: pointer; font-size: 1rem; }
    .erro { min-height: 1.3em; margin-top: 12px; color: #ff9292; }
  </style>
</head>
<body>
  <form class="login" method="post" action="/login">
    <h1>🧬 Simulador Molecular</h1>
    <p>Este site é privado. Digite a senha para continuar.</p>
    <label for="senha">Senha</label>
    <input id="senha" name="senha" type="password" required autofocus>
    <div class="erro">__ERRO__</div>
    <button type="submit">Entrar</button>
  </form>
</body>
</html>
"""


def pagina_login(erro=""):
    return LOGIN_HTML.replace("__ERRO__", erro)


PUBLIC_PATHS = {"/login", "/healthz", "/favicon.ico"}


@app.middleware("http")
async def exigir_login(request: Request, call_next):
    caminho = request.url.path
    if caminho not in PUBLIC_PATHS and not request.session.get("logged_in"):
        if caminho.startswith("/api/"):
            return JSONResponse(
                {"detail": "Faça login para usar o simulador."}, status_code=401
            )
        return RedirectResponse("/login", status_code=303)
    return await call_next(request)


# Deve ser adicionado DEPOIS do middleware de autenticação para que a sessão
# seja o middleware externo e esteja disponível em request.session.
COOKIE_SECURE = os.environ.get("COOKIE_SECURE", "0").lower() in {"1", "true", "yes"}
app.add_middleware(
    SessionMiddleware,
    secret_key=SESSION_SECRET,
    same_site="lax",
    https_only=COOKIE_SECURE,
    max_age=60 * 60 * 24 * 30,
)


@app.get("/login", response_class=HTMLResponse)
async def tela_login(request: Request):
    if request.session.get("logged_in"):
        return RedirectResponse("/", status_code=303)
    return HTMLResponse(pagina_login())


@app.post("/login", response_class=HTMLResponse)
async def fazer_login(request: Request, senha: str = Form(...)):
    if secrets.compare_digest(senha, APP_PASSWORD):
        request.session["logged_in"] = True
        return RedirectResponse("/", status_code=303)
    return HTMLResponse(pagina_login("Senha incorreta."), status_code=401)


@app.get("/logout")
async def fazer_logout(request: Request):
    request.session.clear()
    return RedirectResponse("/login", status_code=303)


@app.get("/healthz")
async def healthz():
    return PlainTextResponse("ok")

app.mount("/outputs", StaticFiles(directory=OUTPUT_DIR), name="outputs")

JOBS: dict[str, dict] = {}
_lock = threading.Lock()


class VirusCustomizado(BaseModel):
    nome: str = ""
    preset: str = ""                  # "covid19" | "hiv" | "influenza" | "" (vazio = customizado)
    tipo_acido: str = "RNA"           # "RNA" ou "DNA" (só usado se preset vazio)
    sequencia_acido: str = ""          # texto cru digitado pelo usuário (só se preset vazio)
    raio_envelope: float = 6.0
    n_espiculas: int = 30
    usar_imagem_referencia: bool = False   # ajusta cor/alongamento com foto real (Wikimedia)
    nome_organismo_imagem: str = ""        # termo de busca da imagem (ex.: organismo do resultado escolhido)


class BacteriaGenomaCustomizado(BaseModel):
    tipo_acido: str = "DNA"          # "DNA" ou "RNA"
    sequencia_acido: str = ""         # texto cru (digitado ou vindo da busca no NCBI)
    usar_imagem_referencia: bool = False   # ajusta cor/proporção da 1ª bactéria com foto real
    nome_organismo_imagem: str = ""        # termo de busca da imagem (ex.: organismo do resultado escolhido)


class SimularRequest(BaseModel):
    # ---- proteína nova por combinação de aminoácidos ------------------
    comprimento_proteina: int = 30       # nº total de resíduos da cadeia
    aminoacidos_pool: list[str] = Field(default_factory=list)     # AAs permitidos ([] = todos os 20)
    composicao_fixa: dict[str, int] = Field(default_factory=dict) # {"ALA":4,"LEU":6} — ignora comprimento/pool
    n_proteinas: int = 1                 # quantas proteínas gerar (1-5)

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
    n_mutacoes: int = 3


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
        with _lock:
            JOBS[job_id]["estado"] = "concluido"
            JOBS[job_id]["resultado"] = resultado
    except Exception as e:
        with _lock:
            JOBS[job_id]["estado"] = "erro"
            JOBS[job_id]["erro"] = f"{e}\n{traceback.format_exc()}"


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
def simular(req: SimularRequest):
    job_id = uuid.uuid4().hex[:12]
    with _lock:
        JOBS[job_id] = {"estado": "executando", "log": [], "resultado": None, "erro": None}
    config = req.model_dump() if hasattr(req, "model_dump") else req.dict()
    t = threading.Thread(target=_rodar_job, args=(job_id, config), daemon=True)
    t.start()
    return {"job_id": job_id}


@app.get("/api/status/{job_id}")
def status(job_id: str):
    with _lock:
        job = JOBS.get(job_id)
    if job is None:
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
    out["bacterias"] = [{**b, "glb": full(b.get("glb"))} for b in resultado.get("bacterias", [])]
    out["docking"] = [{**d, "glb": full(d.get("glb"))} for d in resultado.get("docking", [])]
    out["mutacoes"] = [{**m, "glb": full(m.get("glb"))} for m in resultado.get("mutacoes", [])]
    mt = resultado.get("mutacao_top")
    out["mutacao_top"] = {**mt, "glb": full(mt.get("glb"))} if mt else None
    return out





# Serve o frontend (index.html, style.css, app.js) na raiz "/" —
# como tudo está na mesma pasta, servimos a própria BASE_DIR.
# O StaticFiles precisa ser montado por ÚLTIMO, senão ele "engole" as
# rotas /api/... definidas acima.
app.mount("/", StaticFiles(directory=WEB_DIR, html=True), name="frontend")