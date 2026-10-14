"""Audit of numerical output and provenance; never a biological certificate."""
import hashlib
import importlib.metadata
import json
import math
import platform
from datetime import datetime, timezone
from pathlib import Path


def nonfinite_paths(value, path="resultado"):
    if isinstance(value, float) and not math.isfinite(value):
        return [path]
    if isinstance(value, dict):
        return [p for k, v in value.items() for p in nonfinite_paths(v, f"{path}.{k}")]
    if isinstance(value, (list, tuple)):
        return [p for i, v in enumerate(value) for p in nonfinite_paths(v, f"{path}[{i}]")]
    return []


def build_quality_report(config, result, source_dir):
    invalid = nonfinite_paths(result)
    versions = {}
    for package in ("numpy", "torch", "trimesh", "rdkit", "fastapi", "pydantic"):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = None
    hashes = {}
    for name in ("server.py", "simulacao.py", "scientific_quality.py"):
        path = Path(source_dir) / name
        if path.is_file():
            hashes[name] = hashlib.sha256(path.read_bytes()).hexdigest()
    return {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "integridade_numerica": "reprovada" if invalid else "aprovada",
        "campos_nao_finitos": invalid,
        "validade_biologica": "nao_demonstrada",
        "metodo": "modelo geometrico e escore heuristico",
        "unidade_escore": "arbitraria; nao kcal/mol nem probabilidade",
        "configuracao": config,
        "python": platform.python_version(),
        "plataforma": platform.platform(),
        "dependencias": versions,
        "sha256_codigo": hashes,
        "limitacoes": [
            "Integridade numérica avalia apenas valores finitos, não correção física.",
            "A checagem final não examina todas as coordenadas intermediárias ou a qualidade das malhas.",
            "Não há validação experimental, campo de força atomístico validado, solvente explícito ou análise de convergência.",
            "Não há garantia de reprodução bit a bit: sementes e fontes externas não são integralmente controladas.",
            "Scores relativos não medem eficácia, afinidade experimental ou efeito de variantes.",
            "Modelos visuais e cores de micrografias não validam estrutura molecular.",
            "Composição de aminoácidos não preserva sequência; fragmentos importados não equivalem ao genoma completo.",
        ],
        "referencias": [
            "https://www.wwpdb.org/validation/validation-reports",
            "https://manual.gromacs.org/current/user-guide/",
        ],
    }


def write_quality_report(report, path):
    # Config validation rejects NaN; nonfinite result values are recorded as paths only.
    Path(path).write_text(json.dumps(report, ensure_ascii=False, indent=2,
                                     allow_nan=False), encoding="utf-8")
