# Simulador Molecular protegido por senha

## Rodar localmente

Use Python 3.11 ou mais recente. No terminal, dentro desta pasta:

### Windows PowerShell

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
$env:APP_PASSWORD = "uma-senha-local"
$env:SESSION_SECRET = "gere-uma-chave-longa-e-aleatoria"
$env:COOKIE_SECURE = "0"
uvicorn server:app --reload --port 8000
```

### Linux/macOS

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
export APP_PASSWORD='uma-senha-local'
export SESSION_SECRET="$(python -c 'import secrets; print(secrets.token_urlsafe(32))')"
export COOKIE_SECURE=0
uvicorn server:app --reload --port 8000
```

Abra http://127.0.0.1:8000. A tela de senha deve aparecer antes da interface.

## Publicar no Render

Envie esta pasta para um repositório privado no GitHub. No Render, crie um **Web Service** e use:

```text
Build command: pip install -r requirements.txt
Start command: uvicorn server:app --host 0.0.0.0 --port $PORT
Health check path: /healthz
```

Crie no Render as variáveis secretas:

```text
APP_PASSWORD = a senha escolhida por você
SESSION_SECRET = uma chave aleatória longa
COOKIE_SECURE = 1
```

Não salve essas variáveis no GitHub. O arquivo `server.py` protege a tela, as rotas `/api/...` e os arquivos gerados em `/outputs`.

## Novos controles do simulador

- **Limite energético:** na interface, ative o corte e informe a energia máxima. O modelo usa unidades arbitrárias; energia menor ou igual ao corte é considerada dentro do limite. São tentados alguns candidatos e o de menor energia é escolhido.
- **Tipo de resposta:** a interface oferece análogo de nucleotídeo/ácido nucleico, inibidor de capsômero/proteína de superfície e coquetéis virais, bacterianos ou amplos. A seleção organiza a interpretação do docking para vírus e bactérias; não é um preditor clínico.
- Para que moléculas reais também participem do ranking, ative **Incluir compostos reais**. Caso contrário, o candidato principal é a proteína gerada pelo próprio simulador.

## Observações

- Use apenas um worker/processo, pois os jobs são armazenados em memória (`JOBS`).
- `scipy` é necessário pelo `trimesh` para corrigir as normais das malhas 3D e já está incluído no `requirements.txt`.
- O armazenamento local do Render pode ser apagado quando o serviço reiniciar. Para guardar resultados permanentemente, será necessário um disco persistente ou armazenamento externo.
- A simulação depende de `torch`, `trimesh` e, para fármacos, `rdkit`; a instalação pode ser pesada. O Dockerfile incluído pode ser usado em plataformas que aceitem Docker.
- O site e a API estão no mesmo serviço; por isso `const API = ""` permanece no `app.js`.
