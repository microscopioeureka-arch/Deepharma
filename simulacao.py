# -*- coding: utf-8 -*-
"""
simulacao.py
==============================================================================
Núcleo do SIMULADOR DE DOCKING PROTEICO + MODELOS VIRAIS FICTÍCIOS +
TRIAGEM DE FÁRMACO, refatorado a partir do script original em forma de
FUNÇÃO reutilizável (executar_simulacao), para poder ser chamado tanto:

  - localmente no VS Code (veja rodar_local.py), ou
  - a partir de um servidor web (veja server.py / app FastAPI),

que expõe ao usuário três controles principais para a PROTEÍNA HÍBRIDA:

  1. DOMÍNIOS   -> quantos "domínios" de aminoácidos (sequências geradas
                   aleatoriamente a partir da tabela AMINOACIDOS) compõem
                   a proteína híbrida, e o tamanho de cada domínio.
                   Isso é ANÁLOGO a uma "fusion protein" real: dois ou
                   mais domínios de proteínas diferentes unidos por um
                   peptídeo "linker" (aqui usamos o motivo clássico
                   glicina-serina, "GGGGS", muito usado em fusion
                   proteins de verdade por ser curto e flexível).
  2. QUANTIDADE -> quantas proteínas híbridas gerar (n_hibridos).
  3. AMBIENTE   -> em qual "ambiente" a proteína híbrida está sendo
                   dobrada/produzida (ver ENVIRONMENTS abaixo) — muda os
                   pesos físicos (eletrostático x hidrofóbico x contato
                   com água) usados na física de dobramento, simulando,
                   por exemplo, uma membrana lipídica (mais apolar) vs.
                   citoplasma aquoso vs. meio ácido/alcalino/salino.

Os COMPOSTOS REAIS (via RDKit) continuam disponíveis na simulação, mas
agora entram apenas como LIGANTES/FÁRMACOS na etapa de triagem/docking
contra os vírus fictícios (`incluir_compostos_reais=True`) — não são
mais usados como "matéria-prima" da proteína híbrida. Isso é mais fiel
à biologia real: uma fusion protein é feita de domínios de aminoácidos,
enquanto fármacos pequenos (Paracetamol, Cafeína etc.) são ligantes que
se encaixam num bolso de ligação, não blocos de construção covalentes
do esqueleto da proteína.

Tudo o mais (as 20 tabelas de aminoácidos, os 5 modelos de vírus
fictícios, os 8 compostos reais via RDKit, a física de dobramento e o
docking de corpo rígido) é o MESMO do script original — só foi
organizado em funções e ganhou parâmetros em vez de rodar direto no
nível do módulo.

DIFERENÇAS DELIBERADAS EM RELAÇÃO AO SCRIPT ORIGINAL (para uso web):
  - `modo_rapido=True` (padrão no servidor web) reduz nº de resíduos,
    iterações de dobramento e nº de starts do docking, e pula o vídeo
    de fusão (frame a frame) e o ESMFold (chamada de API externa lenta)
    por padrão — para responder em segundos/poucos minutos em vez de
    dezenas de minutos. Tudo isso é configurável.
  - A triagem (seção 14 original) agora testa, por padrão, só a(s)
    proteína(s) híbrida(s) geradas pelo usuário contra os 5 vírus
    fictícios (em vez de também rodar os 8 compostos reais "puros" —
    isso pode ser reativado com `incluir_compostos_reais=True`).
"""

from __future__ import annotations

import os
import re
import json
import time
import shutil
import datetime
import unicodedata
from urllib.parse import quote
import numpy as np
import torch
import trimesh

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    MATPLOTLIB_DISPONIVEL = True
except ImportError:
    MATPLOTLIB_DISPONIVEL = False

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, Crippen, Lipinski, rdMolDescriptors
    RDKIT_DISPONIVEL = True
except ImportError:
    RDKIT_DISPONIVEL = False

try:
    import requests
except ImportError:
    requests = None

try:
    import io
    from PIL import Image
    PIL_DISPONIVEL = True
except ImportError:
    PIL_DISPONIVEL = False


# =====================================================================
# 1. TABELAS FÍSICO-QUÍMICAS (idênticas ao script original)
# =====================================================================
AMINOACIDOS = {
    "ALA": {"hidrofobicidade": 1.8,  "carga": 0.0,  "raio_vdw": 2.30, "massa": 71.08},
    "ARG": {"hidrofobicidade": -4.5, "carga": 1.0,  "raio_vdw": 3.35, "massa": 156.19},
    "ASN": {"hidrofobicidade": -3.5, "carga": 0.0,  "raio_vdw": 2.72, "massa": 114.10},
    "ASP": {"hidrofobicidade": -3.5, "carga": -1.0, "raio_vdw": 2.66, "massa": 115.09},
    "CYS": {"hidrofobicidade": 2.5,  "carga": 0.0,  "raio_vdw": 2.57, "massa": 103.14},
    "GLN": {"hidrofobicidade": -3.5, "carga": 0.0,  "raio_vdw": 2.94, "massa": 128.13},
    "GLU": {"hidrofobicidade": -3.5, "carga": -1.0, "raio_vdw": 2.87, "massa": 129.12},
    "GLY": {"hidrofobicidade": -0.4, "carga": 0.0,  "raio_vdw": 2.10, "massa": 57.05},
    "HIS": {"hidrofobicidade": -3.2, "carga": 0.1,  "raio_vdw": 3.03, "massa": 137.14},
    "ILE": {"hidrofobicidade": 4.5,  "carga": 0.0,  "raio_vdw": 3.09, "massa": 113.16},
    "LEU": {"hidrofobicidade": 3.8,  "carga": 0.0,  "raio_vdw": 3.09, "massa": 113.16},
    "LYS": {"hidrofobicidade": -3.9, "carga": 1.0,  "raio_vdw": 3.18, "massa": 128.17},
    "MET": {"hidrofobicidade": 1.9,  "carga": 0.0,  "raio_vdw": 3.09, "massa": 131.19},
    "PHE": {"hidrofobicidade": 2.8,  "carga": 0.0,  "raio_vdw": 3.28, "massa": 147.18},
    "PRO": {"hidrofobicidade": -1.6, "carga": 0.0,  "raio_vdw": 2.78, "massa": 97.12},
    "SER": {"hidrofobicidade": -0.8, "carga": 0.0,  "raio_vdw": 2.46, "massa": 87.08},
    "THR": {"hidrofobicidade": -0.7, "carga": 0.0,  "raio_vdw": 2.69, "massa": 101.10},
    "TRP": {"hidrofobicidade": -0.9, "carga": 0.0,  "raio_vdw": 3.54, "massa": 186.21},
    "TYR": {"hidrofobicidade": -1.3, "carga": 0.0,  "raio_vdw": 3.34, "massa": 163.18},
    "VAL": {"hidrofobicidade": 4.2,  "carga": 0.0,  "raio_vdw": 2.93, "massa": 99.13},
}
LISTA_AA = list(AMINOACIDOS.keys())
TRES_PARA_UMA = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
DIST_CA_CA = 3.8

COR1_HIDROFOBICO = np.array([240, 80, 80, 255], dtype=float)
COR1_HIDROFILICO = np.array([80, 240, 150, 255], dtype=float)
COR_LIGACAO = np.array([210, 210, 210, 255], dtype=np.uint8)

# Constantes usadas pelas funções antigas de proteína híbrida. Elas não
# eram necessárias no caminho novo da simulação web, mas ficam definidas
# para que essas funções também possam ser chamadas sem NameError.
LINKER_GS = ["GLY", "GLY", "GLY", "GLY", "SER"]
COR_LINKER = np.array([170, 175, 185, 255], dtype=np.uint8)
PALETA_DOMINIOS = [
    np.array([100, 180, 255, 255], dtype=np.uint8),
    np.array([255, 160, 90, 255], dtype=np.uint8),
    np.array([170, 110, 240, 255], dtype=np.uint8),
    np.array([100, 220, 150, 255], dtype=np.uint8),
]

PALETA_ENVELOPES = [
    np.array([80, 200, 160, 50], dtype=np.uint8),
    np.array([90, 140, 220, 50], dtype=np.uint8),
    np.array([220, 140, 80, 50], dtype=np.uint8),
    np.array([180, 90, 200, 50], dtype=np.uint8),
    np.array([220, 200, 65, 50], dtype=np.uint8),
]
PALETA_INTERNA = [
    {"spike": np.array([255, 100, 120, 255], dtype=np.uint8),
     "cabeca": np.array([255, 170, 185, 255], dtype=np.uint8),
     "glicano": np.array([215, 255, 210, 235], dtype=np.uint8),
     "nucleo": np.array([120, 40, 90, 255], dtype=np.uint8)},
    {"spike": np.array([120, 220, 255, 255], dtype=np.uint8),
     "cabeca": np.array([170, 235, 255, 255], dtype=np.uint8),
     "glicano": np.array([255, 240, 180, 235], dtype=np.uint8),
     "nucleo": np.array([20, 70, 110, 255], dtype=np.uint8)},
    {"spike": np.array([255, 210, 90, 255], dtype=np.uint8),
     "cabeca": np.array([255, 230, 150, 255], dtype=np.uint8),
     "glicano": np.array([200, 255, 230, 235], dtype=np.uint8),
     "nucleo": np.array([120, 70, 10, 255], dtype=np.uint8)},
    {"spike": np.array([190, 120, 255, 255], dtype=np.uint8),
     "cabeca": np.array([220, 175, 255, 255], dtype=np.uint8),
     "glicano": np.array([255, 220, 245, 255], dtype=np.uint8),
     "nucleo": np.array([55, 15, 90, 255], dtype=np.uint8)},
    {"spike": np.array([120, 255, 170, 255], dtype=np.uint8),
     "cabeca": np.array([180, 255, 210, 255], dtype=np.uint8),
     "glicano": np.array([255, 245, 200, 255], dtype=np.uint8),
     "nucleo": np.array([10, 90, 55, 255], dtype=np.uint8)},
]

CORES_CPK = {
    "C": np.array([120, 120, 120, 255], dtype=np.uint8),
    "N": np.array([50, 90, 240, 255], dtype=np.uint8),
    "O": np.array([230, 50, 50, 255], dtype=np.uint8),
    "H": np.array([235, 235, 235, 255], dtype=np.uint8),
    "S": np.array([230, 220, 50, 255], dtype=np.uint8),
    "F": np.array([120, 220, 120, 255], dtype=np.uint8),
    "Cl": np.array([60, 200, 60, 255], dtype=np.uint8),
    "P": np.array([255, 140, 0, 255], dtype=np.uint8),
}
COR_CPK_PADRAO = np.array([200, 100, 220, 255], dtype=np.uint8)

COMPOSTOS_REAIS = {
    "Paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "Ibuprofeno": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "Cafeina": "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
    "Aspirina": "CC(=O)Oc1ccccc1C(=O)O",
    "Favipiravir": "NC(=O)c1nc(F)cnc1O",
    "Ribavirina": "C1=NC(=NN1C2C(C(C(O2)CO)O)O)C(=O)N",
    "Oseltamivir": "CCC(CC)OC1C=C(CC(C1NC(=O)C)N)C(=O)OCC",
    "Cloroquina": "CCN(CC)CCCC(C)Nc1ccnc2cc(Cl)ccc12",
}
RAIO_COVALENTE = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05,
                   "F": 0.57, "Cl": 1.02, "P": 1.07}
MASSA_ATOMICA = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "S": 32.06,
                  "F": 18.998, "Cl": 35.45, "P": 30.974}
HIDRO_PROXY_ELEMENTO = {"C": 1.0, "S": 0.8, "H": 0.2, "N": -1.0, "O": -1.2,
                         "F": -0.6, "Cl": -0.4, "P": -0.3}

# Cores para visualização de resíduos: hidrofóbico (quente) / hidrofílico (frio)
COR_HIDROFOBO   = np.array([240, 80,  80,  255], dtype=float)   # vermelho
COR_HIDROFILICO = np.array([80,  200, 140, 255], dtype=float)   # verde

# ---------------------------------------------------------------------
# AMBIENTES disponíveis para a produção/dobramento da proteína híbrida.
# Cada ambiente escala os pesos da física de dobramento (mola + VdW +
# Coulomb + colapso hidrofóbico + afinidade por água), simulando
# condições fisico-químicas diferentes.
# ---------------------------------------------------------------------
ENVIRONMENTS = {
    "citoplasma": {
        "label": "Citoplasma (aquoso, pH neutro)",
        "k_mola": 1.0, "k_vdw": 1.0, "k_eletrostatico": 3.0,
        "k_hidrofobico": 1.5, "k_hidrofilico_agua": 0.4,
    },
    "membrana": {
        "label": "Membrana lipídica (ambiente apolar)",
        "k_mola": 1.0, "k_vdw": 1.2, "k_eletrostatico": 1.5,
        "k_hidrofobico": 3.0, "k_hidrofilico_agua": 0.15,
    },
    "acido": {
        "label": "Meio ácido (ex. endossomo, pH ~5)",
        "k_mola": 1.0, "k_vdw": 1.0, "k_eletrostatico": 4.5,
        "k_hidrofobico": 1.2, "k_hidrofilico_agua": 0.5,
    },
    "alcalino": {
        "label": "Meio alcalino (pH ~9)",
        "k_mola": 1.0, "k_vdw": 1.0, "k_eletrostatico": 2.0,
        "k_hidrofobico": 1.7, "k_hidrofilico_agua": 0.35,
    },
    "salino": {
        "label": "Solução salina concentrada (blindagem eletrostática)",
        "k_mola": 1.0, "k_vdw": 1.0, "k_eletrostatico": 0.8,
        "k_hidrofobico": 1.8, "k_hidrofilico_agua": 0.3,
    },
}

# Estratégias de resposta são categorias DIDÁTICAS para interpretar o
# docking. Elas não transformam o simulador em um preditor clínico nem
# afirmam que um composto selecionado é realmente eficaz contra um
# patógeno. A mesma seleção pode ser aplicada a vírus e bactérias:
# para vírus, "capsômero" descreve a estrutura viral; para bactérias,
# a categoria é interpretada como proteína/envelope de superfície.
TIPOS_RESPOSTA = {
    "nenhuma": {
        "label": "Somente gerar e ranquear a proteína",
        "descricao": "Não aplicar uma estratégia de resposta; mostra apenas o docking ilustrativo.",
        "familia": "nenhuma",
    },
    "analogo_nucleotideo": {
        "label": "Análogo de nucleotídeo / bloqueio de ácido nucleico",
        "descricao": "Interpreta a proteína como candidata a interferir na replicação de RNA/DNA.",
        "familia": "acido_nucleico",
    },
    "inibidor_capsomero_superficie": {
        "label": "Inibidor de capsômero / proteína de superfície",
        "descricao": "Para vírus, representa um alvo estrutural do capsômero; para bactérias, uma proteína/envelope de superfície.",
        "familia": "estrutura_superficie",
    },
    "coquetel_antiviral": {
        "label": "Coquetel antiviral (ácido nucleico + estrutura)",
        "descricao": "Compara uma resposta combinada contra alvos virais de ácido nucleico e superfície.",
        "familia": "coquetel_viral",
    },
    "coquetel_antibacteriano": {
        "label": "Coquetel antibacteriano (ácido nucleico + superfície)",
        "descricao": "Compara uma resposta combinada contra alvos estruturais e de ácido nucleico bacterianos.",
        "familia": "coquetel_bacteriano",
    },
    "coquetel_amplo": {
        "label": "Coquetel amplo (vírus e bactérias)",
        "descricao": "Aplica a mesma leitura combinada aos alvos virais e bacterianos presentes na simulação.",
        "familia": "coquetel_amplo",
    },
}

VIRUS_MODELOS_BASE = [
    {"nome": "Virus-Ficticio-A", "seed": 101, "cor_env_idx": 0,
     "raio_envelope": 6.0, "n_res_spike": 15, "n_bases_rna": 22,
     "n_espiculas": 55, "raio_espicula": 0.30, "altura_espicula": 0.9,
     "raio_cabeca_espicula": 0.34, "n_glicanos_por_espicula": 3,
     "amplitude_organica": 0.06, "forma_alongada": False, "fator_alongamento": 1.0},
    {"nome": "Virus-Ficticio-B", "seed": 202, "cor_env_idx": 1,
     "raio_envelope": 4.8, "n_res_spike": 20, "n_bases_rna": 30,
     "n_espiculas": 15, "raio_espicula": 0.45, "altura_espicula": 1.6,
     "raio_cabeca_espicula": 0.55, "n_glicanos_por_espicula": 4,
     "amplitude_organica": 0.08, "forma_alongada": False, "fator_alongamento": 1.0},
    {"nome": "Virus-Ficticio-C", "seed": 303, "cor_env_idx": 2,
     "raio_envelope": 3.2, "n_res_spike": 10, "n_bases_rna": 16,
     "n_espiculas": 0, "raio_espicula": 0.0, "altura_espicula": 0.0,
     "raio_cabeca_espicula": 0.0, "n_glicanos_por_espicula": 0,
     "amplitude_organica": 0.10, "forma_alongada": True, "fator_alongamento": 2.6},
    {"nome": "Virus-Ficticio-D", "seed": 404, "cor_env_idx": 3,
     "raio_envelope": 7.5, "n_res_spike": 24, "n_bases_rna": 34,
     "n_espiculas": 90, "raio_espicula": 0.22, "altura_espicula": 0.6,
     "raio_cabeca_espicula": 0.24, "n_glicanos_por_espicula": 2,
     "amplitude_organica": 0.05, "forma_alongada": False, "fator_alongamento": 1.0},
    {"nome": "Virus-Ficticio-E", "seed": 505, "cor_env_idx": 4,
     "raio_envelope": 5.5, "n_res_spike": 15, "n_bases_rna": 24,
     "n_espiculas": 8, "raio_espicula": 0.55, "altura_espicula": 2.0,
     "raio_cabeca_espicula": 0.62, "n_glicanos_por_espicula": 5,
     "amplitude_organica": 0.09, "forma_alongada": True, "fator_alongamento": 1.4},
]


# =====================================================================
# 2. HELPERS
# =====================================================================
def slug(texto: str) -> str:
    texto = texto.split(" (")[0]
    for orig, rep in [("á", "a"), ("â", "a"), ("ã", "a"), ("é", "e"), ("ê", "e"),
                       ("í", "i"), ("ó", "o"), ("ô", "o"), ("õ", "o"), ("ú", "u"), ("ç", "c")]:
        texto = texto.replace(orig, rep)
    return re.sub(r"[^a-zA-Z0-9]+", "_", texto).strip("_").lower()


def gerar_sequencia(n: int, seed=None):
    gen = torch.Generator().manual_seed(seed) if seed is not None else None
    idx = torch.randint(0, len(LISTA_AA), (n,), generator=gen)
    return [LISTA_AA[i] for i in idx]


def sequencia_para_fasta(seq):
    return "".join(TRES_PARA_UMA[a] for a in seq)


def _matriz_rotacao(v_origem, v_destino):
    a = v_origem / np.linalg.norm(v_origem)
    b = v_destino / np.linalg.norm(v_destino)
    cos_ang = float(np.clip(np.dot(a, b), -1.0, 1.0))
    if cos_ang > 1.0 - 1e-8:
        return np.eye(4)
    if cos_ang < -1.0 + 1e-8:
        eixo = np.cross(a, [1.0, 0.0, 0.0])
        if np.linalg.norm(eixo) < 1e-6:
            eixo = np.cross(a, [0.0, 1.0, 0.0])
        eixo = eixo / np.linalg.norm(eixo)
        return trimesh.transformations.rotation_matrix(np.pi, eixo)
    eixo = np.cross(a, b) / np.linalg.norm(np.cross(a, b))
    return trimesh.transformations.rotation_matrix(np.arccos(cos_ang), eixo)


def criar_cilindro_ligacao(p1, p2, cor=COR_LIGACAO, raio=0.18):
    vetor = p2 - p1
    h = float(np.linalg.norm(vetor))
    if h < 1e-6:
        return None
    cil = trimesh.creation.cylinder(radius=raio, height=h, sections=12)
    cil.apply_transform(_matriz_rotacao(np.array([0.0, 0.0, 1.0]), vetor))
    cil.apply_translation((p1 + p2) / 2)
    cil.visual.vertex_colors = cor
    return cil


def raio_por_aa(nome_aa):
    return 0.55 + (AMINOACIDOS[nome_aa]["raio_vdw"] / 3.54) * 0.5


def aa_e_hidrofobico(nome_aa):
    return AMINOACIDOS[nome_aa]["hidrofobicidade"] > 0.0


def criar_geometria_aminoacido(p, nome_aa, cor_hb, cor_hf, destaque=False, cor_destaque=None):
    # destaque=True amplia o resíduo e usa uma cor de aviso — usado para
    # marcar visualmente posições MUTADAS num vírus (ver prever_mutacao_virus).
    raio = raio_por_aa(nome_aa) * (1.45 if destaque else 1.0)
    esfera = trimesh.creation.icosphere(radius=raio, subdivisions=3)
    if destaque:
        cor = cor_destaque if cor_destaque is not None else np.array([255, 225, 0, 255], dtype=np.uint8)
    else:
        cor = cor_hb if aa_e_hidrofobico(nome_aa) else cor_hf
    esfera.visual.vertex_colors = cor.astype(np.uint8)
    esfera.apply_translation(p)
    return esfera


def _catmull_rom(pontos, pontos_por_segmento=6):
    """Interpola uma curva suave (spline Catmull-Rom) passando por
    TODOS os pontos de entrada, na ordem. Usada para desenhar o
    esqueleto (backbone) de uma cadeia de aminoácidos como um tubo
    curvo e contínuo, em vez de segmentos retos entre cada resíduo —
    visualmente muito mais parecido com o traçado de fita/backbone
    usado em visualizadores moleculares reais (PyMOL, ChimeraX) do
    que uma sequência de bolinhas ligadas por palitos retos."""
    pts = np.asarray(pontos, dtype=float)
    n = len(pts)
    if n < 3:
        return pts
    pts_ext = np.vstack([pts[0] + (pts[0] - pts[1]), pts, pts[-1] + (pts[-1] - pts[-2])])
    curva = []
    for i in range(1, n):
        p0, p1, p2, p3 = pts_ext[i - 1], pts_ext[i], pts_ext[i + 1], pts_ext[i + 2]
        ts = np.linspace(0, 1, pontos_por_segmento, endpoint=(i == n - 1))
        for t in ts:
            t2, t3 = t * t, t * t * t
            ponto = 0.5 * ((2 * p1) + (-p0 + p2) * t
                            + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t2
                            + (-p0 + 3 * p1 - 3 * p2 + p3) * t3)
            curva.append(ponto)
    return np.array(curva)


COR_BACKBONE = np.array([196, 200, 210, 255], dtype=np.uint8)


def criar_backbone_suave(pontos, cor=None, raio=0.13):
    """Desenha o backbone (esqueleto carbono-alfa) da cadeia como um
    tubo curvo contínuo — spline Catmull-Rom + cilindros/esferas de
    junção — em vez de segmentos retos ligando cada resíduo. Retorna
    uma lista de geometrias trimesh a serem somadas à cena."""
    if len(pontos) < 2:
        return []
    if cor is None:
        cor = COR_BACKBONE
    curva = _catmull_rom(pontos, pontos_por_segmento=6) if len(pontos) >= 3 else np.asarray(pontos, dtype=float)
    geos = []
    for i in range(len(curva) - 1):
        cil = criar_cilindro_ligacao(curva[i], curva[i + 1], cor=cor, raio=raio)
        if cil is not None:
            geos.append(cil)
        if 0 < i < len(curva) - 1:
            junta = trimesh.creation.icosphere(radius=raio, subdivisions=1)
            junta.visual.vertex_colors = cor
            junta.apply_translation(curva[i])
            geos.append(junta)
    return geos


def montar_cadeia(pontos, sequencia, cor_hb, cor_hf, indices_destaque=None, cor_destaque=None):
    indices_destaque = set(indices_destaque or [])
    geos = criar_backbone_suave(pontos, raio=0.12)
    geos.extend(criar_geometria_aminoacido(p, sequencia[i], cor_hb, cor_hf,
                                            destaque=(i in indices_destaque), cor_destaque=cor_destaque)
                for i, p in enumerate(pontos))
    return geos


def gerar_pontos_enrolados(num_res, raio_interno, voltas=3.5, fase=0.0, offset=np.zeros(3)):
    u = np.linspace(0, 1, num_res)
    theta = fase + voltas * 2 * np.pi * u
    z = raio_interno * (1 - 2 * u)
    r_xy = np.sqrt(np.maximum(raio_interno ** 2 - z ** 2, 0.0))
    return np.stack([r_xy * np.cos(theta), r_xy * np.sin(theta), z], axis=1) + offset


def exportar_glb(scene, caminho):
    try:
        import pygltflib  # noqa: F401
    except ImportError:
        return None
    try:
        scene.export(caminho)
        if os.path.exists(caminho) and os.path.getsize(caminho) > 0:
            return caminho
        return None
    except Exception as e:
        print(f"[ERRO] Falha ao exportar .glb: {e}")
        return None


# =====================================================================
# 3. FÍSICA DE DOBRAMENTO
# =====================================================================
def _loss_backbone(p, k_mola, dist_alvo=DIST_CA_CA):
    d = torch.norm(p[1:] - p[:-1], dim=1)
    return k_mola * torch.mean((d - dist_alvo) ** 2)


def _loss_vdw(p, raios, k_vdw):
    d = torch.cdist(p, p)
    sr = raios.unsqueeze(0) + raios.unsqueeze(1)
    mask = ~torch.eye(p.shape[0], dtype=torch.bool)
    return k_vdw * torch.sum((torch.relu(sr - d) * mask) ** 2)


def _loss_coulomb(p, cargas, k_ele):
    d = torch.cdist(p, p) + 1e-2
    q = cargas.unsqueeze(0) * cargas.unsqueeze(1)
    mask = ~torch.eye(p.shape[0], dtype=torch.bool)
    return k_ele * torch.sum((q / d) * mask)


def _loss_hidrofobico(p, hidro, k_hidro):
    centro = p.mean(dim=0)
    dist_c = torch.norm(p - centro, dim=1)
    return k_hidro * torch.mean(torch.relu(hidro) * dist_c)


def _loss_hidrofilico_agua(p, hidro, k_agua):
    centro = p.mean(dim=0)
    dist_c = torch.norm(p - centro, dim=1)
    return k_agua * torch.mean(torch.relu(-hidro) * (1.0 / (dist_c + 1e-3)))


def loss_cadeia(p, raios, cargas, hidro, k):
    """k: dict com k_mola, k_vdw, k_eletrostatico, k_hidrofobico, k_hidrofilico_agua

    Usa DIST_CA_CA (3.8 Å) como distância-alvo do "backbone" — a mesma
    distância real entre carbonos-alfa consecutivos numa cadeia
    polipeptídica. É a física usada tanto para as espículas dos vírus
    fictícios quanto (agora) para a proteína híbrida baseada em
    domínios de aminoácidos.
    """
    return (
        _loss_backbone(p, k["k_mola"])
        + _loss_vdw(p, raios, k["k_vdw"])
        + _loss_coulomb(p, cargas, k["k_eletrostatico"])
        + _loss_hidrofobico(p, hidro, k["k_hidrofobico"])
        + _loss_hidrofilico_agua(p, hidro, k["k_hidrofilico_agua"])
    )


def _loss_repulsao_inter(p1, p2, r1, r2):
    d = torch.cdist(p1, p2)
    sr = r1.unsqueeze(1) + r2.unsqueeze(0)
    return torch.sum(torch.relu(sr - d) ** 2)


def _loss_atracao_inter(p1, p2, c1, c2, h1, h2):
    d = torch.cdist(p1, p2) + 1e-2
    eletro = torch.sum((c1.unsqueeze(1) * c2.unsqueeze(0)) / d)
    hidro = -0.1 * torch.sum((torch.relu(h1).unsqueeze(1) * torch.relu(h2).unsqueeze(0)) / d)
    return eletro + hidro


def _matriz_rotacao_torch(eixo_angulo):
    theta = torch.norm(eixo_angulo) + 1e-8
    k = eixo_angulo / theta
    zero = torch.zeros((), dtype=eixo_angulo.dtype)
    K = torch.stack([
        torch.stack([zero, -k[2], k[1]]),
        torch.stack([k[2], zero, -k[0]]),
        torch.stack([-k[1], k[0], zero]),
    ])
    I = torch.eye(3, dtype=eixo_angulo.dtype)
    return I + torch.sin(theta) * K + (1 - torch.cos(theta)) * (K @ K)


def dockar_molecula_corpo_rigido(pontos_mol, cargas_mol, raios_mol, hidro_mol,
                                  p_virus, carga_virus, raio_virus, hidro_virus,
                                  n_iter=200, posicao_inicial_x=9.0, seed=0, n_starts=3):
    gerador = torch.Generator().manual_seed(seed)
    pontos_local = torch.tensor(pontos_mol, dtype=torch.float32)
    pontos_local = pontos_local - pontos_local.mean(dim=0)
    cargas_t = torch.tensor(cargas_mol, dtype=torch.float32)
    raios_t = torch.tensor(raios_mol, dtype=torch.float32)
    hidro_t = torch.tensor(hidro_mol, dtype=torch.float32)

    melhor_loss_global = float("inf")
    melhor_pontos_global = None

    for i_start in range(max(1, n_starts)):
        if i_start == 0:
            translacao = torch.tensor([posicao_inicial_x, 0.0, 0.0], requires_grad=True)
            rotacao = torch.zeros(3, requires_grad=True)
        else:
            direcao = torch.randn(3, generator=gerador)
            direcao = direcao / (torch.norm(direcao) + 1e-8)
            translacao = (direcao * posicao_inicial_x).clone().detach().requires_grad_(True)
            rotacao = ((torch.rand(3, generator=gerador) - 0.5) * 2 * np.pi).clone().detach().requires_grad_(True)

        opt = torch.optim.Adam([translacao, rotacao], lr=0.05)
        melhor_loss = float("inf")
        melhor_pontos = None
        sem_melhora = 0
        for _ in range(n_iter):
            opt.zero_grad()
            R = _matriz_rotacao_torch(rotacao)
            p_mundo = pontos_local @ R.T + translacao
            loss = (_loss_repulsao_inter(p_mundo, p_virus, raios_t, raio_virus)
                    + _loss_atracao_inter(p_mundo, p_virus, cargas_t, carga_virus, hidro_t, hidro_virus))
            loss.backward()
            opt.step()
            loss_atual = float(loss.item())
            if loss_atual < melhor_loss - 1e-4:
                melhor_loss = loss_atual
                melhor_pontos = p_mundo.detach().numpy().copy()
                sem_melhora = 0
            else:
                sem_melhora += 1
                if sem_melhora >= 20:
                    break
        if melhor_pontos is None:
            with torch.no_grad():
                R = _matriz_rotacao_torch(rotacao)
                melhor_pontos = (pontos_local @ R.T + translacao).numpy()
                melhor_loss = 0.0
        if melhor_loss < melhor_loss_global:
            melhor_loss_global = melhor_loss
            melhor_pontos_global = melhor_pontos

    return melhor_pontos_global, melhor_loss_global


# =====================================================================
# 4. VÍRUS FICTÍCIOS (geometria)
# =====================================================================
def _perturbar_membrana_organica(mesh, seed, amplitude=0.07, n_ondas=4):
    rng = np.random.default_rng(seed)
    verts = mesh.vertices.copy()
    raio_base = np.linalg.norm(verts, axis=1)
    normed = verts / raio_base[:, None]
    theta = np.arccos(np.clip(normed[:, 2], -1.0, 1.0))
    phi = np.arctan2(normed[:, 1], normed[:, 0])
    fator = np.zeros(len(verts))
    for _ in range(n_ondas):
        l = rng.integers(2, 5)
        m = rng.integers(1, 4)
        fase = rng.uniform(0, 2 * np.pi)
        amp_i = amplitude * rng.uniform(0.35, 1.0) / n_ondas
        fator += amp_i * np.cos(l * theta + fase) * np.cos(m * phi + fase)
    novo_raio = raio_base * (1.0 + fator)
    mesh.vertices = normed * novo_raio[:, None]
    mesh.fix_normals()
    return mesh


def criar_envelope_viral_organico(raio, cor, seed, amplitude=0.07):
    envelope = trimesh.creation.icosphere(radius=raio, subdivisions=3)
    envelope = _perturbar_membrana_organica(envelope, seed, amplitude=amplitude)
    rgba_float = cor.astype(float) / 255.0
    material = trimesh.visual.material.PBRMaterial(
        baseColorFactor=tuple(rgba_float), alphaMode="BLEND",
        doubleSided=True, metallicFactor=0.0, roughnessFactor=0.7,
    )
    envelope.visual = trimesh.visual.TextureVisuals(material=material)
    return envelope


def criar_bicamada_lipidica(raio, cor, seed, amplitude=0.07, espessura_relativa=0.055):
    externa = criar_envelope_viral_organico(raio, cor, seed, amplitude=amplitude)
    raio_interno = raio * (1.0 - espessura_relativa)
    cor_interna = cor.copy().astype(float)
    cor_interna[:3] *= 0.6
    cor_interna[3] = min(255, cor_interna[3] * 1.6 + 15)
    cor_interna = cor_interna.astype(np.uint8)
    interna = criar_envelope_viral_organico(raio_interno, cor_interna, seed + 1, amplitude=amplitude * 0.9)
    return [externa, interna]


BASES_RNA = ["A", "U", "G", "C"]
CORES_RNA = {
    "A": np.array([255, 99, 71, 255], dtype=np.uint8),
    "U": np.array([65, 105, 225, 255], dtype=np.uint8),
    "G": np.array([60, 179, 113, 255], dtype=np.uint8),
    "C": np.array([238, 210, 2, 255], dtype=np.uint8),
}
COR_BACKBONE_RNA = np.array([200, 200, 200, 255], dtype=np.uint8)


def gerar_rna_virus(seed, raio_nucleo, num_bases=24):
    rng = np.random.default_rng(seed + 9999)
    bases = list(rng.choice(BASES_RNA, size=num_bases))
    raio_rna = raio_nucleo * rng.uniform(0.35, 0.55)
    voltas = rng.uniform(2.0, 4.5)
    fase = rng.uniform(0, 2 * np.pi)
    offset = rng.normal(scale=raio_nucleo * 0.15, size=3)
    pontos = gerar_pontos_enrolados(num_bases, raio_rna, voltas=voltas, fase=fase, offset=offset)
    return bases, pontos


def criar_geometria_rna(bases, pontos):
    geos = []
    for p, base in zip(pontos, bases):
        esfera = trimesh.creation.icosphere(radius=0.32, subdivisions=1)
        esfera.visual.vertex_colors = CORES_RNA[base]
        esfera.apply_translation(p)
        geos.append(esfera)
    for i in range(len(pontos) - 1):
        cil = criar_cilindro_ligacao(pontos[i], pontos[i + 1], cor=COR_BACKBONE_RNA, raio=0.06)
        if cil is not None:
            geos.append(cil)
    return geos


# ---------------------------------------------------------------------
# ÁCIDO NUCLEICO CUSTOMIZADO — o usuário digita a própria sequência de
# RNA (bases A/U/G/C) ou DNA (bases A/T/G/C) em vez de usar uma sorteada
# aleatoriamente. Usado tanto para "vírus customizados" quanto para o
# nucleoide das bactérias.
# ---------------------------------------------------------------------
CORES_DNA = {
    "A": np.array([255, 99, 71, 255], dtype=np.uint8),
    "T": np.array([65, 105, 225, 255], dtype=np.uint8),
    "G": np.array([60, 179, 113, 255], dtype=np.uint8),
    "C": np.array([238, 210, 2, 255], dtype=np.uint8),
}


def parse_sequencia_acido_nucleico(texto: str, tipo: str = "RNA") -> list:
    """Converte o texto digitado pelo usuário (ex.: 'AUGCAUGC...' para
    RNA, 'ATGCATGC...' para DNA) numa lista de bases válidas, ignorando
    espaços, quebras de linha e caracteres inválidos."""
    texto = (texto or "").upper().strip()
    bases_validas = {"A", "U", "G", "C"} if tipo == "RNA" else {"A", "T", "G", "C"}
    bases = [b for b in texto if b in bases_validas]
    if len(bases) < 3:
        raise ValueError(f"Sequência de {tipo} muito curta ou inválida (mínimo 3 bases).")
    return bases


def criar_geometria_acido_nucleico(bases, pontos, tipo="RNA"):
    paleta = CORES_RNA if tipo == "RNA" else CORES_DNA
    geos = []
    for p, base in zip(pontos, bases):
        esfera = trimesh.creation.icosphere(radius=0.32, subdivisions=1)
        esfera.visual.vertex_colors = paleta.get(base, COR_BACKBONE_RNA)
        esfera.apply_translation(p)
        geos.append(esfera)
    for i in range(len(pontos) - 1):
        cil = criar_cilindro_ligacao(pontos[i], pontos[i + 1], cor=COR_BACKBONE_RNA, raio=0.06)
        if cil is not None:
            geos.append(cil)
    return geos


def pontos_fibonacci_esfera(n, seed=0, jitter=0.03):
    if n <= 0:
        return np.zeros((0, 3))
    rng = np.random.default_rng(seed)
    indices = np.arange(n) + 0.5
    phi = np.arccos(np.clip(1 - 2 * indices / n, -1.0, 1.0))
    golden_angle = np.pi * (3 - np.sqrt(5))
    theta = golden_angle * indices
    pts = np.stack([np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi)], axis=1)
    pts = pts + rng.normal(scale=jitter, size=pts.shape)
    pts = pts / np.linalg.norm(pts, axis=1, keepdims=True)
    return pts


def relaxar_pontos_esfera(pontos_unitarios, n_iter=35, taxa=0.05):
    n = len(pontos_unitarios)
    if n <= 1:
        return pontos_unitarios
    pts = pontos_unitarios.copy()
    for _ in range(n_iter):
        diffs = pts[:, None, :] - pts[None, :, :]
        d2 = np.sum(diffs ** 2, axis=-1)
        np.fill_diagonal(d2, np.inf)
        d2 = np.maximum(d2, 1e-4)
        forca = np.sum(diffs / (d2[..., None] ** 1.5), axis=1)
        pts = pts + taxa * forca
        pts = pts / np.linalg.norm(pts, axis=1, keepdims=True)
    return pts


def criar_espicula_realista(direcao_unitaria, raio_esfera, altura_total, raio_haste,
                             raio_cabeca, cor_haste, cor_cabeca, escala_posicao=None):
    if escala_posicao is None:
        escala_posicao = np.array([1.0, 1.0, 1.0])
    inicio = direcao_unitaria * raio_esfera * escala_posicao
    fim_haste = direcao_unitaria * (raio_esfera + altura_total * 0.62) * escala_posicao
    centro_cabeca = direcao_unitaria * (raio_esfera + altura_total * 0.62 + raio_cabeca * 0.55) * escala_posicao

    partes = []
    haste = criar_cilindro_ligacao(inicio, fim_haste, cor=cor_haste, raio=raio_haste)
    if haste is not None:
        partes.append(haste)

    cabeca = trimesh.creation.icosphere(radius=raio_cabeca, subdivisions=1)
    S = np.eye(4)
    S[0, 0] = S[1, 1] = 1.15
    S[2, 2] = 0.85
    cabeca.apply_transform(S)
    direcao_real = centro_cabeca - fim_haste
    norma = np.linalg.norm(direcao_real)
    direcao_real = direcao_unitaria if norma < 1e-6 else direcao_real / norma
    cabeca.apply_transform(_matriz_rotacao(np.array([0.0, 0.0, 1.0]), direcao_real))
    cabeca.apply_translation(centro_cabeca)
    cabeca.visual.vertex_colors = cor_cabeca.astype(np.uint8)
    partes.append(cabeca)
    return trimesh.util.concatenate(partes)


def criar_glicanos(centro_cabeca, raio_cabeca, n_glicanos, seed, cor):
    if n_glicanos <= 0:
        return []
    rng = np.random.default_rng(seed)
    geos = []
    for _ in range(n_glicanos):
        offset = rng.normal(scale=raio_cabeca * 0.85, size=3)
        p = centro_cabeca + offset
        s = trimesh.creation.icosphere(radius=raio_cabeca * 0.24, subdivisions=1)
        s.visual.vertex_colors = cor.astype(np.uint8)
        s.apply_translation(p)
        geos.append(s)
    return geos


def criar_nucleo_interno(raio, cor):
    nucleo = trimesh.creation.icosphere(radius=raio, subdivisions=1)
    nucleo.visual.vertex_colors = cor.astype(np.uint8)
    return nucleo


def montar_virus_3d(sequencia, seed, estilo, centro=np.zeros(3), genoma_customizado=None,
                     indices_mutados=None):
    """
    genoma_customizado: opcional, {"tipo": "RNA"|"DNA", "bases": [...]}.
    Quando fornecido, usa a sequência DIGITADA PELO USUÁRIO no lugar do
    RNA sorteado aleatoriamente (gerar_rna_virus) — assim dá pra montar
    um vírus fictício com genoma "não pré-definido".

    indices_mutados: opcional, lista de índices (0-based) da sequência da
    espícula (spike) a destacar visualmente em amarelo/laranja — usado
    para mostrar ONDE uma mutação prevista (ver prever_mutacao_virus)
    caiu na proteína de superfície.
    """
    e = estilo
    raio_envelope = e["raio_envelope"]
    raio_nucleo = e["raio_nucleo"]
    amplitude_organica = e.get("amplitude_organica", 0.07)
    geos = []
    geos.extend(criar_bicamada_lipidica(raio_envelope, e["cor_envelope"], seed, amplitude=amplitude_organica))

    cor_nucleo_interno = e.get("cor_nucleo_interno", np.array([90, 90, 90, 255], dtype=np.uint8))
    geos.append(criar_nucleo_interno(raio_nucleo * 0.5, cor_nucleo_interno))

    cor_spike = e.get("cor_spike", np.array([255, 100, 120, 255], dtype=np.uint8))
    cor_cabeca_espicula = e.get("cor_cabeca_espicula_visual", np.array([255, 170, 185, 255], dtype=np.uint8))
    cor_glicano = e.get("cor_glicano", np.array([215, 255, 210, 235], dtype=np.uint8))

    pontos_spike = gerar_pontos_enrolados(len(sequencia), raio_nucleo)
    cor_mutacao = np.array([255, 210, 0, 255], dtype=np.uint8)
    geos.extend(montar_cadeia(pontos_spike, sequencia, cor_spike, cor_spike,
                               indices_destaque=indices_mutados, cor_destaque=cor_mutacao))

    if genoma_customizado:
        bases_ac = genoma_customizado["bases"]
        tipo_ac = genoma_customizado.get("tipo", "RNA")
        raio_ac = raio_nucleo * 0.45
        voltas_ac = max(2.0, len(bases_ac) / 12)
        pontos_ac = gerar_pontos_enrolados(len(bases_ac), raio_ac, voltas=voltas_ac)
        geos.extend(criar_geometria_acido_nucleico(bases_ac, pontos_ac, tipo=tipo_ac))
    else:
        bases_rna, pontos_rna = gerar_rna_virus(seed, raio_nucleo, e["num_bases_rna"])
        geos.extend(criar_geometria_rna(bases_rna, pontos_rna))

    n_esp = e.get("n_espiculas", 0)
    if n_esp > 0:
        direcoes = relaxar_pontos_esfera(pontos_fibonacci_esfera(n_esp, seed=seed))
        raio_haste = e.get("raio_haste_espicula", e["raio_espicula"] * 0.35)
        raio_cabeca = e.get("raio_cabeca_espicula", e["raio_espicula"])
        n_glicanos = e.get("n_glicanos_por_espicula", 3)
        for i, d in enumerate(direcoes):
            geos.append(criar_espicula_realista(d, raio_envelope, e["altura_espicula"],
                                                 raio_haste, raio_cabeca, cor_spike, cor_cabeca_espicula))
            centro_cabeca = d * (raio_envelope + e["altura_espicula"] * 0.62 + raio_cabeca * 0.55)
            geos.extend(criar_glicanos(centro_cabeca, raio_cabeca, n_glicanos, seed=seed * 1000 + i, cor=cor_glicano))

    if e.get("forma_alongada", False):
        fator = e.get("fator_alongamento", 1.6)
        S = np.eye(4)
        S[2, 2] = fator
        for g in geos:
            g.apply_transform(S)
    for g in geos:
        g.apply_translation(centro)
    return geos, pontos_spike


# =====================================================================
# 4B. BACTÉRIAS FICTÍCIAS — corpo (bacilo/coco) + parede + nucleoide +
#     flagelo + uma proteína de superfície (adesina/pilina) usada como
#     alvo de docking, igual ao spike dos vírus.
# =====================================================================
def mutar_sequencia(sequencia, n_mutacoes, seed=None):
    """
    Sorteia `n_mutacoes` posições (sem repetição) na sequência de
    aminoácidos e troca cada uma por um resíduo diferente, simulando
    mutações pontuais (substituições, análogas a SNPs não-sinônimos
    num gene viral real). Retorna (nova_sequencia, lista_de_mutacoes),
    onde cada mutação é {"posicao","de","para","de_1letra","para_1letra"}.
    """
    rng = np.random.default_rng(seed)
    seq = list(sequencia)
    n = len(seq)
    n_mut = max(1, min(int(n_mutacoes), n))
    posicoes = rng.choice(n, size=n_mut, replace=False)
    mutacoes = []
    for idx in sorted(int(i) for i in posicoes):
        aa_orig = seq[idx]
        candidatos = [aa for aa in LISTA_AA if aa != aa_orig]
        aa_novo = str(rng.choice(candidatos))
        seq[idx] = aa_novo
        mutacoes.append({
            "posicao": idx + 1, "de": aa_orig, "para": aa_novo,
            "de_1letra": TRES_PARA_UMA[aa_orig], "para_1letra": TRES_PARA_UMA[aa_novo],
        })
    return seq, mutacoes


def classificar_mutacao_ponto(aa_orig, aa_novo):
    """
    Índice heurístico (0-100) de "impacto estrutural estimado" de uma
    mutação pontual, combinando a variação de hidrofobicidade, carga e
    raio de van der Waals entre o resíduo original e o novo — um proxy
    simples (não um preditor biológico validado) para o quanto a
    mutação pode alterar a superfície de ligação da proteína.
    """
    d_hidro = abs(AMINOACIDOS[aa_novo]["hidrofobicidade"] - AMINOACIDOS[aa_orig]["hidrofobicidade"])
    d_carga = abs(AMINOACIDOS[aa_novo]["carga"] - AMINOACIDOS[aa_orig]["carga"])
    d_raio = abs(AMINOACIDOS[aa_novo]["raio_vdw"] - AMINOACIDOS[aa_orig]["raio_vdw"])
    indice = min(100.0, (d_hidro / 8.7) * 45.0 + d_carga * 35.0 + (d_raio / 1.5) * 20.0)
    if indice < 15:
        classe = "conservativa"
    elif indice < 40:
        classe = "moderada"
    else:
        classe = "radical"
    return round(float(indice), 1), classe


def avaliar_mutacoes(mutacoes):
    """Preenche cada mutação com seu índice/classe e retorna o resumo
    geral (índice médio + classificação) do conjunto de mutações."""
    if not mutacoes:
        return 0.0, "sem mutações"
    for m in mutacoes:
        indice, classe = classificar_mutacao_ponto(m["de"], m["para"])
        m["indice_impacto"] = indice
        m["classe"] = classe
    media = round(sum(m["indice_impacto"] for m in mutacoes) / len(mutacoes), 1)
    if media < 15:
        geral = "conservativa — baixo impacto estimado na ligação"
    elif media < 40:
        geral = "moderada — pode alterar a afinidade de ligação"
    else:
        geral = "radical — alto risco estimado de escape/resistência"
    return media, geral


def prever_mutacao_virus(v, n_mutacoes, seed_mutacao, dir_saida, run_id, modo_rapido=True):
    """
    Gera uma VARIANTE MUTANTE do vírus `v` (dict de runtime, já com
    'sequencia', 'seed' e 'estilo' preenchidos por executar_simulacao):
    sorteia mutações pontuais na proteína de espícula, remonta o
    modelo 3D destacando os resíduos mutados (amarelo) e calcula um
    índice heurístico de impacto. NÃO refaz o docking (isso é feito à
    parte, só para o par vencedor, em executar_simulacao) — este índice
    é um proxy físico-químico rápido, ilustrativo, sem validade
    biológica real.
    """
    seq_mutada, mutacoes = mutar_sequencia(v["sequencia"], n_mutacoes, seed=seed_mutacao)
    indice_medio, classificacao = avaliar_mutacoes(mutacoes)
    indices_pos = [m["posicao"] - 1 for m in mutacoes]

    genoma_customizado = v.get("_genoma_customizado")
    geos, pontos_spike_mut = montar_virus_3d(
        seq_mutada, v["seed"], v["estilo"], genoma_customizado=genoma_customizado,
        indices_mutados=indices_pos,
    )
    scene = trimesh.Scene()
    for g in geos:
        scene.add_geometry(g)
    scene.apply_translation(-scene.centroid)
    glb = exportar_glb(scene, os.path.join(dir_saida, f"mutante_{slug(v['nome'])}_{run_id}.glb"))

    return {
        "virus": v["nome"],
        "mutacoes": mutacoes,
        "indice_impacto_medio": indice_medio,
        "classificacao": classificacao,
        "sequencia_original_fasta": sequencia_para_fasta(v["sequencia"]),
        "sequencia_mutada_fasta": sequencia_para_fasta(seq_mutada),
        "glb": glb,
        "_seq_mutada": seq_mutada,
        "_pontos_spike_mutados": pontos_spike_mut,
    }


def gerar_pontos_helice(n, raio, passo, comprimento, offset=np.zeros(3)):
    """Pontos de uma hélice ao longo do eixo Z local (usada no flagelo)."""
    voltas = comprimento / max(passo, 1e-3)
    t = np.linspace(0, voltas * 2 * np.pi, n)
    x = raio * np.cos(t)
    y = raio * np.sin(t)
    z = np.linspace(0, comprimento, n)
    return np.stack([x, y, z], axis=1) + offset


def criar_corpo_basal(origem, direcao_unitaria, cor):
    """
    Corpo basal simplificado do flagelo bacteriano: anéis empilhados
    ancorados na parede celular (representação estilizada dos anéis
    L/P (membrana externa/peptideoglicano) e MS/C (membrana interna) do
    motor flagelar real), que giram para produzir o movimento do
    filamento. Retorna (geometrias, ponto_de_saida) — o segundo valor é
    onde o gancho/filamento deve começar.
    """
    geos = []
    raios_aneis = [0.24, 0.19, 0.14]
    alturas = [0.05, 0.09, 0.07]
    z_acum = 0.0
    for raio_anel, altura in zip(raios_aneis, alturas):
        centro = origem + direcao_unitaria * (z_acum + altura / 2)
        disco = trimesh.creation.cylinder(radius=raio_anel, height=altura, sections=16)
        disco.apply_transform(_matriz_rotacao(np.array([0.0, 0.0, 1.0]), direcao_unitaria))
        disco.apply_translation(centro)
        disco.visual.vertex_colors = cor
        geos.append(disco)
        z_acum += altura
    return geos, origem + direcao_unitaria * z_acum


def criar_flagelo(origem, direcao_unitaria, comprimento, raio_helice, raio_tubo, cor,
                   n_pontos=48, passo=1.3):
    """
    Flagelo bacteriano com anatomia simplificada mas mais fiel à real:
    corpo basal (anéis ancorados na parede) -> gancho curto (curvatura
    mais fechada, faz a transição de eixo) -> filamento helicoidal
    afunilado (mais fino perto da ponta), com juntas arredondadas entre
    segmentos em vez de cilindros "facetados" soltos.
    """
    cor_basal = np.array([175, 175, 190, 255], dtype=np.uint8)
    geos, origem_gancho = criar_corpo_basal(origem, direcao_unitaria, cor_basal)

    # gancho: trecho curto com curvatura mais acentuada (raio maior,
    # passo mais apertado) que conecta o motor basal ao filamento
    R = _matriz_rotacao(np.array([0.0, 0.0, 1.0]), direcao_unitaria)[:3, :3]
    n_gancho = 8
    comp_gancho = max(comprimento * 0.10, 0.3)
    pts_gancho = gerar_pontos_helice(n_gancho, raio_helice * 1.35, passo * 0.5, comp_gancho) @ R.T + origem_gancho
    for i in range(len(pts_gancho) - 1):
        cil = criar_cilindro_ligacao(pts_gancho[i], pts_gancho[i + 1], cor=cor_basal, raio=raio_tubo * 1.3)
        if cil is not None:
            geos.append(cil)

    # filamento: hélice longa e afunilada (mais fina perto da ponta),
    # com uma pequena esfera em cada junta pra suavizar o traçado.
    comp_filamento = comprimento * 0.86
    origem_filamento = pts_gancho[-1] if len(pts_gancho) else origem_gancho
    pontos_mundo = gerar_pontos_helice(n_pontos, raio_helice, passo, comp_filamento) @ R.T + origem_filamento
    for i in range(len(pontos_mundo) - 1):
        frac = i / max(1, len(pontos_mundo) - 2)
        raio_local = raio_tubo * (1.0 - 0.35 * frac)
        cil = criar_cilindro_ligacao(pontos_mundo[i], pontos_mundo[i + 1], cor=cor, raio=raio_local)
        if cil is not None:
            geos.append(cil)
        if i > 0:
            junta = trimesh.creation.icosphere(radius=raio_local, subdivisions=1)
            junta.visual.vertex_colors = cor
            junta.apply_translation(pontos_mundo[i])
            geos.append(junta)
    return geos


def gerar_pontos_esparsos_no_corpo(n, raio_max_xy, meia_altura, forma, seed,
                                    excluir_centro_raio=0.0):
    """
    Amostra `n` pontos aleatórios DENTRO do volume do corpo bacteriano
    (esfera para 'coco', cilindro com raio constante para 'bacilo'),
    opcionalmente evitando uma zona central cilíndrica (usada para não
    sobrepor o nucleoide) — usado para espalhar ribossomos e plasmídeos
    pelo citoplasma de forma orgânica, não numa grade.
    """
    rng = np.random.default_rng(seed)
    pontos = []
    tentativas = 0
    while len(pontos) < n and tentativas < n * 40:
        tentativas += 1
        if forma == "coco":
            v = rng.normal(size=3)
            v = v / (np.linalg.norm(v) + 1e-8)
            r = raio_max_xy * (rng.uniform(0.12, 0.92) ** (1 / 3))
            p = v * r
        else:
            z = rng.uniform(-meia_altura, meia_altura)
            ang = rng.uniform(0, 2 * np.pi)
            r = raio_max_xy * np.sqrt(rng.uniform(0.05, 0.90))
            p = np.array([r * np.cos(ang), r * np.sin(ang), z])
        if excluir_centro_raio > 0 and np.linalg.norm(p[:2]) < excluir_centro_raio and abs(p[2]) < meia_altura * 0.6:
            continue
        pontos.append(p)
    return np.array(pontos) if pontos else np.zeros((0, 3))


def criar_ribossomos(centros, raio, cor):
    """Ribossomos livres no citoplasma — pequenas partículas escuras.
    Bactérias (procariontes) não têm retículo endoplasmático: os
    ribossomos ficam soltos no citoplasma, por isso são só pontos
    dispersos, sem organela associada."""
    geos = []
    for c in centros:
        esf = trimesh.creation.icosphere(radius=raio, subdivisions=0)
        esf.visual.vertex_colors = cor
        esf.apply_translation(c)
        geos.append(esf)
    return geos


def gerar_pontos_anel(n, raio, seed=0, offset=np.zeros(3), deform=0.10):
    """Pontos aproximadamente circulares (com leve ondulação orgânica)
    no plano local — usado para o DNA circular dos plasmídeos, que na
    biologia real são mesmo anéis fechados de DNA."""
    rng = np.random.default_rng(seed)
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    x = raio * np.cos(t)
    y = raio * np.sin(t)
    z = rng.normal(scale=raio * deform, size=n)
    return np.stack([x, y, z], axis=1) + offset


def criar_plasmideo(centro, raio, n_bases, seed, cor):
    """Plasmídeo: pequeno anel EXTRA de DNA circular, independente do
    nucleoide principal — em bactérias reais costuma carregar genes
    acessórios (ex.: resistência a antibióticos)."""
    pontos = gerar_pontos_anel(n_bases, raio, seed=seed, offset=centro)
    geos = []
    for p in pontos:
        esf = trimesh.creation.icosphere(radius=raio * 0.11, subdivisions=0)
        esf.visual.vertex_colors = cor
        esf.apply_translation(p)
        geos.append(esf)
    n = len(pontos)
    for i in range(n):
        cil = criar_cilindro_ligacao(pontos[i], pontos[(i + 1) % n], cor=cor, raio=raio * 0.045)
        if cil is not None:
            geos.append(cil)
    return geos


PALETA_BACTERIAS = [
    {"parede": np.array([120, 200, 140, 60], dtype=np.uint8),
     "citoplasma": np.array([70, 120, 90, 220], dtype=np.uint8),
     "proteina_superficie": np.array([255, 150, 90, 255], dtype=np.uint8)},
    {"parede": np.array([200, 150, 90, 60], dtype=np.uint8),
     "citoplasma": np.array([140, 90, 50, 220], dtype=np.uint8),
     "proteina_superficie": np.array([120, 190, 255, 255], dtype=np.uint8)},
    {"parede": np.array([170, 120, 210, 60], dtype=np.uint8),
     "citoplasma": np.array([100, 60, 130, 220], dtype=np.uint8),
     "proteina_superficie": np.array([210, 255, 130, 255], dtype=np.uint8)},
    {"parede": np.array([90, 170, 200, 60], dtype=np.uint8),
     "citoplasma": np.array([50, 100, 130, 220], dtype=np.uint8),
     "proteina_superficie": np.array([255, 210, 90, 255], dtype=np.uint8)},
    {"parede": np.array([210, 110, 130, 60], dtype=np.uint8),
     "citoplasma": np.array([140, 60, 75, 220], dtype=np.uint8),
     "proteina_superficie": np.array([150, 255, 200, 255], dtype=np.uint8)},
    {"parede": np.array([190, 190, 110, 60], dtype=np.uint8),
     "citoplasma": np.array([120, 120, 60, 220], dtype=np.uint8),
     "proteina_superficie": np.array([180, 150, 255, 255], dtype=np.uint8)},
]

NOMES_BACTERIAS = ["Alfa", "Beta", "Gama", "Delta", "Epsilon", "Zeta", "Eta", "Theta"]


def gerar_estilo_bacteria(indice: int, seed_base: int = 600) -> dict:
    """
    Gera os parâmetros de UMA bactéria fictícia de forma procedural e
    determinística (mesmo índice + mesma seed_base => sempre a mesma
    bactéria), para que as N bactérias pedidas pelo usuário sejam
    REALMENTE N modelos diferentes entre si — forma (bacilo/coco
    alternados), tamanho do corpo, nº de flagelos, nº de plasmídeos,
    quantidade de ribossomos e paleta de cores — em vez de sempre
    reciclar um número fixo de modelos pré-definidos.

    Antes desta função, `BACTERIA_MODELOS_BASE` só tinha 2 entradas
    fixas, então pedir 2+ bactérias sempre devolvia os MESMOS 2
    modelos; agora dá pra gerar quantas forem pedidas, cada uma com
    sua própria "personalidade" morfológica.
    """
    seed = seed_base + indice * 101
    rng = np.random.default_rng(seed)
    forma = "bacilo" if indice % 2 == 0 else "coco"
    if forma == "bacilo":
        raio_corpo = round(float(rng.uniform(1.0, 1.75)), 2)
        comprimento = round(float(rng.uniform(4.0, 8.5)), 2)
    else:
        raio_corpo = round(float(rng.uniform(1.6, 2.7)), 2)
        comprimento = 0.0
    nome_sufixo = NOMES_BACTERIAS[indice % len(NOMES_BACTERIAS)]
    return {
        "nome": f"Bacteria-Ficticia-{nome_sufixo}",
        "seed": seed,
        "forma": forma,
        "comprimento": comprimento,
        "raio_corpo": raio_corpo,
        "n_flagelos": int(rng.integers(1, 5)),
        "comprimento_flagelo": round(float(rng.uniform(3.0, 7.0)), 2),
        "raio_helice_flagelo": round(float(rng.uniform(0.30, 0.50)), 2),
        "n_bases_dna": int(rng.integers(24, 48)),
        "n_res_proteina_superficie": int(rng.integers(8, 17)),
        "paleta_idx": indice % len(PALETA_BACTERIAS),
        "n_ribossomos": int(rng.integers(50, 115)),
        "n_plasmideos": int(rng.integers(0, 4)),
    }


def calcular_estatisticas_bacteria(b: dict) -> dict:
    """Resumo numérico de UMA bactéria já gerada — usado pelo frontend
    para mostrar os detalhes visuais/estruturais (forma, tamanho,
    flagelos, plasmídeos, ribossomos etc.) igual já é feito para as
    proteínas."""
    return {
        "forma": b["forma"],
        "raio_corpo": b["raio_corpo"],
        "comprimento": b["comprimento"],
        "n_flagelos": b["n_flagelos"],
        "n_plasmideos": b["n_plasmideos"],
        "n_ribossomos": b["n_ribossomos"],
        "n_bases_dna": b.get("n_bases_dna"),
        "n_residuos_proteina_superficie": len(b.get("sequencia_proteina_superficie", [])),
        "genoma_customizado": bool(b.get("_genoma_customizado")),
    }

# =====================================================================
# PRESETS DE VÍRUS REAIS — trechos curtos e REAIS (não o genoma inteiro,
# que tem dezenas de milhares de bases e inviabilizaria a visualização)
# das primeiras ~90 bases das sequências de referência públicas do
# NCBI/GenBank, convertidas de DNA (como estão depositadas) para RNA
# (T -> U), já que os 3 são vírus de RNA. Fonte / accession de cada uma:
#   - COVID-19 (SARS-CoV-2): NCBI RefSeq NC_045512.2 (Wuhan-Hu-1)
#   - HIV-1: NCBI RefSeq NC_001802.1 (isolado HXB2)
#   - Influenza A H1N1: GenBank CY083910.1 (gene da hemaglutinina, HA)
# Uso educacional/ilustrativo — não é o genoma completo do vírus.
# =====================================================================
VIRUS_REAIS_PRESETS = {
    "covid19": {
        "nome": "SARS-CoV-2 (COVID-19)",
        "tipo_acido": "RNA",
        "sequencia_acido": (
            "AUUAAAGGUUUAUACCUUCCCAGGUAACAAACCAACCAACUUUCGAUCUCUUGUAGAUCU"
            "GUUCUCUAAACGAACUUUAAAAUCUGUGUG"
        ),
        "fonte": "NCBI RefSeq NC_045512.2 (primeiras ~90 bases)",
    },
    "hiv": {
        "nome": "HIV-1",
        "tipo_acido": "RNA",
        "sequencia_acido": (
            "GGUCUCUCUGGUUAGACCAGAUCUGAGCCUGGGAGCUCUCUGGCUAACUAGGGAACCCA"
            "CUGCUUAAGCCUCAAUAAAGCUUGCCUUGAG"
        ),
        "fonte": "NCBI RefSeq NC_001802.1, isolado HXB2 (primeiras ~90 bases)",
    },
    "influenza": {
        "nome": "Influenza A (H1N1)",
        "tipo_acido": "RNA",
        "sequencia_acido": (
            "GGAAAACAAAAGCAACAAAAAUGAAGGCAAUACUAGUAGUUCUGCUAUAUACAUUUGUA"
            "ACCGCAAAUGCAGACACAUUAUGUAUAGGUU"
        ),
        "fonte": "GenBank CY083910.1, gene da hemaglutinina (primeiras ~90 bases)",
    },
}

BACTERIA_MODELOS_BASE = [
    {"nome": "Bacteria-Ficticia-A", "seed": 601, "forma": "bacilo",
     "comprimento": 6.0, "raio_corpo": 1.4, "n_flagelos": 1,
     "comprimento_flagelo": 6.0, "raio_helice_flagelo": 0.4,
     "n_bases_dna": 40, "n_res_proteina_superficie": 12, "paleta_idx": 0,
     "n_ribossomos": 70, "n_plasmideos": 2},
    {"nome": "Bacteria-Ficticia-B", "seed": 602, "forma": "coco",
     "comprimento": 0.0, "raio_corpo": 2.2, "n_flagelos": 4,
     "comprimento_flagelo": 4.0, "raio_helice_flagelo": 0.35,
     "n_bases_dna": 30, "n_res_proteina_superficie": 10, "paleta_idx": 1,
     "n_ribossomos": 90, "n_plasmideos": 3},
]


def montar_bacteria_3d(sequencia_proteina_superficie, seed, estilo, centro=np.zeros(3),
                        genoma_customizado=None):
    """
    Monta uma bactéria fictícia: corpo em cápsula (bacilo) ou esfera
    (coco), parede celular semitransparente, citoplasma interno,
    nucleoide (DNA circular — pode ser customizado via
    genoma_customizado, {"tipo": "DNA"|"RNA", "bases": [...]}), flagelo(s)
    helicoidais e uma proteína de superfície (equivalente ao spike do
    vírus) usada como alvo no docking/triagem.

    Retorna (geos, pontos_proteina_superficie) — o segundo valor é
    usado do mesmo jeito que pontos_spike dos vírus.
    """
    e = estilo
    forma = e.get("forma", "bacilo")
    raio_corpo = e["raio_corpo"]
    comprimento = e["comprimento"]
    geos = []

    if forma == "coco":
        corpo = trimesh.creation.icosphere(radius=raio_corpo, subdivisions=3)
    else:
        corpo = trimesh.creation.capsule(radius=raio_corpo, height=comprimento, count=[16, 16])
    corpo = _perturbar_membrana_organica(corpo, seed, amplitude=0.025, n_ondas=3)
    rgba = e["cor_parede"].astype(float) / 255.0
    material = trimesh.visual.material.PBRMaterial(
        baseColorFactor=tuple(rgba), alphaMode="BLEND", doubleSided=True, roughnessFactor=0.85)
    corpo.visual = trimesh.visual.TextureVisuals(material=material)
    geos.append(corpo)

    if forma == "coco":
        interno = trimesh.creation.icosphere(radius=raio_corpo * 0.85, subdivisions=2)
    else:
        interno = trimesh.creation.capsule(radius=raio_corpo * 0.85, height=comprimento * 0.9, count=[12, 12])
    interno.visual.vertex_colors = e["cor_citoplasma"]
    geos.append(interno)

    # nucleoide: DNA (ou RNA) circular — customizado se fornecido, senão
    # reaproveita a mesma geometria de fita usada no RNA viral
    if genoma_customizado:
        bases_ac = genoma_customizado["bases"]
        tipo_ac = genoma_customizado.get("tipo", "DNA")
        raio_ac = raio_corpo * 0.45
        voltas_ac = max(2.0, len(bases_ac) / 12)
        pontos_ac = gerar_pontos_enrolados(len(bases_ac), raio_ac, voltas=voltas_ac)
        geos.extend(criar_geometria_acido_nucleico(bases_ac, pontos_ac, tipo=tipo_ac))
    else:
        bases_dna, pontos_dna = gerar_rna_virus(seed, raio_corpo * 0.5, e.get("n_bases_dna", 40))
        geos.extend(criar_geometria_rna(bases_dna, pontos_dna))

    # ribossomos livres no citoplasma. Bactérias são PROCARIONTES: têm
    # ribossomos (menores que os eucarióticos, 70S) soltos no citoplasma,
    # mas NÃO têm mitocôndria, complexo de Golgi, núcleo ou retículo
    # endoplasmático — essas são organelas membranosas exclusivas de
    # células eucarióticas, por isso não entram no modelo.
    n_ribo = int(e.get("n_ribossomos", 60))
    meia_altura_corpo = comprimento / 2 if forma != "coco" else 0.0
    if n_ribo > 0:
        centros_ribo = gerar_pontos_esparsos_no_corpo(
            n_ribo, raio_corpo * 0.80, meia_altura_corpo, forma, seed + 4242,
            excluir_centro_raio=raio_corpo * 0.16)
        cor_ribo = e.get("cor_ribossomo", np.array([55, 40, 72, 255], dtype=np.uint8))
        geos.extend(criar_ribossomos(centros_ribo, raio_corpo * 0.045, cor_ribo))

    # plasmídeos: pequenos anéis extras de DNA circular, independentes
    # do nucleoide principal (estrutura real de muitas bactérias).
    n_plasm = int(e.get("n_plasmideos", 2))
    if n_plasm > 0:
        rng_plasm = np.random.default_rng(seed + 8181)
        cor_plasmideo = e.get("cor_plasmideo", np.array([255, 105, 180, 255], dtype=np.uint8))
        for i in range(n_plasm):
            raio_pl = raio_corpo * rng_plasm.uniform(0.14, 0.22)
            ang = rng_plasm.uniform(0, 2 * np.pi)
            dist = raio_corpo * rng_plasm.uniform(0.30, 0.58)
            z_pl = rng_plasm.uniform(-1, 1) * (meia_altura_corpo * 0.5 if forma != "coco" else raio_corpo * 0.25)
            centro_pl = np.array([dist * np.cos(ang), dist * np.sin(ang), z_pl])
            geos.extend(criar_plasmideo(centro_pl, raio_pl, 10, seed + 8181 + i, cor_plasmideo))

    # proteína de superfície (adesina/pilina) — alvo de docking.
    # FIX: antes essa fita era colocada numa esfera de raio quase igual
    # ao raio do corpo e ainda deslocada no eixo Z (offset) — perto dos
    # polos da cápsula (onde ela é arredondada, mais estreita) a fita
    # acabava ficando PARCIALMENTE FORA da parede celular ("pedaço
    # flutuando"). Agora ela fica centrada na origem (meio do corpo,
    # onde a cápsula é cilíndrica e a "sobra" de raio é constante) e com
    # um raio nitidamente menor que o raio do corpo, com folga de sobra.
    raio_ancoragem = raio_corpo * 0.55
    pontos_proteina = gerar_pontos_enrolados(len(sequencia_proteina_superficie), raio_ancoragem, voltas=2.0)
    cor_prot = e.get("cor_proteina_superficie", np.array([255, 150, 90, 255], dtype=np.uint8))
    geos.extend(montar_cadeia(pontos_proteina, sequencia_proteina_superficie, cor_prot, cor_prot))

    # flagelo(s).
    # FIX: antes a origem do flagelo era "direcao_aleatoria_3D * (comprimento/2)"
    # — isso só dá um ponto na SUPERFÍCIE do corpo quando a forma é uma
    # esfera perfeita (coco). Para o bacilo (cápsula alongada no eixo Z),
    # uma direção aleatória qualquer multiplicada por comprimento/2 cai
    # LONGE da superfície real na maioria dos ângulos (por isso os
    # "pedaços voando" longe do corpo). Agora: no coco, mantém direção
    # aleatória (a esfera é simétrica em qualquer direção, então
    # direcao*raio_corpo cai exatamente na casca). No bacilo, o(s)
    # flagelo(s) nascem perto dos POLOS de verdade da cápsula (topo/base,
    # ao longo do eixo Z), com uma leve inclinação aleatória só pra não
    # ficarem 100% retos.
    n_flag = e.get("n_flagelos", 1)
    rng = np.random.default_rng(seed + 77)
    cor_flagelo = np.array([210, 210, 225, 255], dtype=np.uint8)
    if forma == "coco":
        for _ in range(n_flag):
            direcao = rng.normal(size=3)
            direcao = direcao / (np.linalg.norm(direcao) + 1e-8)
            origem = direcao * raio_corpo  # ponto exato na superfície da esfera
            geos.extend(criar_flagelo(
                origem, direcao, e.get("comprimento_flagelo", 5.0),
                e.get("raio_helice_flagelo", 0.4), 0.08, cor_flagelo))
    else:
        polos = [1.0, -1.0] if n_flag > 1 else [1.0]
        for i in range(n_flag):
            polo = polos[i % len(polos)]
            tilt = rng.normal(scale=0.10, size=3)
            tilt[2] = 0.0
            direcao = np.array([0.0, 0.0, polo]) + tilt
            direcao = direcao / (np.linalg.norm(direcao) + 1e-8)
            # ponto exato no polo da cápsula: metade do comprimento
            # (parte cilíndrica) + raio (calota hemisférica da ponta)
            origem = np.array([0.0, 0.0, polo * (comprimento / 2 + raio_corpo)])
            geos.extend(criar_flagelo(
                origem, direcao, e.get("comprimento_flagelo", 5.0),
                e.get("raio_helice_flagelo", 0.4), 0.08, cor_flagelo))

    for g in geos:
        g.apply_translation(centro)
    return geos, pontos_proteina


# =====================================================================
# 5. COMPOSTOS REAIS (RDKit) — usados como LIGANTES na triagem/docking
# =====================================================================
def calcular_descritores_moleculares(mol):
    return {
        "peso_molecular_da": Descriptors.MolWt(mol),
        "logp": Crippen.MolLogP(mol),
        "tpsa_a2": rdMolDescriptors.CalcTPSA(mol),
        "doadores_h": Lipinski.NumHDonors(mol),
        "aceitadores_h": Lipinski.NumHAcceptors(mol),
        "ligacoes_rotacionaveis": Lipinski.NumRotatableBonds(mol),
        "n_aneis": rdMolDescriptors.CalcNumRings(mol),
        "n_aneis_aromaticos": rdMolDescriptors.CalcNumAromaticRings(mol),
        "carga_formal": Chem.GetFormalCharge(mol),
    }


def carregar_composto_real(nome):
    mol = Chem.MolFromSmiles(COMPOSTOS_REAIS[nome])
    if mol is None:
        raise ValueError(f"SMILES inválido para {nome}.")
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(mol, params) != 0:
        AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    conf = mol.GetConformer()
    pontos = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
    elementos = [atom.GetSymbol() for atom in mol.GetAtoms()]
    ligacoes = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]
    AllChem.ComputeGasteigerCharges(mol)
    cargas = np.array([float(a.GetProp("_GasteigerCharge")) for a in mol.GetAtoms()])
    cargas = np.nan_to_num(cargas, nan=0.0)
    raios = np.array([RAIO_COVALENTE.get(e, 0.75) for e in elementos])
    hidro = np.array([HIDRO_PROXY_ELEMENTO.get(e, 0.0) for e in elementos])
    descritores = calcular_descritores_moleculares(mol)
    return pontos, elementos, ligacoes, cargas, raios, hidro, descritores


def criar_geometria_composto_real(pontos, elementos, ligacoes, escala_raio=0.45):
    geos = []
    for p, el in zip(pontos, elementos):
        cor = CORES_CPK.get(el, COR_CPK_PADRAO)
        raio = RAIO_COVALENTE.get(el, 0.75) * escala_raio
        esfera = trimesh.creation.icosphere(radius=raio, subdivisions=2)
        esfera.visual.vertex_colors = cor.astype(np.uint8)
        esfera.apply_translation(p)
        geos.append(esfera)
    for (i, j) in ligacoes:
        cil = criar_cilindro_ligacao(pontos[i], pontos[j], cor=COR_LIGACAO, raio=0.11)
        if cil is not None:
            geos.append(cil)
    return geos


# =====================================================================
# 6. GERADOR DE PROTEÍNAS — combina aminoácidos e dobra a cadeia
# =====================================================================
def gerar_dominios_hibridos(n_dominios=2, tamanho_dominio=(8, 20), seed=None):  # mantido p/ compatibilidade interna
    """Gera N domínios, cada um uma sequência de aminoácidos aleatória
    (mesma lógica de gerar_sequencia usada nas espículas dos vírus),
    com tamanho sorteado dentro do intervalo tamanho_dominio."""
    rng = np.random.default_rng(seed)
    tam_min, tam_max = tamanho_dominio
    dominios = []
    for _ in range(n_dominios):
        tam = int(rng.integers(tam_min, tam_max + 1))
        seed_dom = int(rng.integers(0, 1_000_000))
        seq = gerar_sequencia(tam, seed=seed_dom)
        dominios.append(seq)
    return dominios


def montar_dominio_por_composicao(composicao: dict, seed=None) -> list:
    """
    composicao: {"ALA": 3, "LEU": 5, "LYS": 2, ...} — o usuário ESCOLHE
    quais aminoácidos entram e QUANTOS de cada um (equivalente a
    "elementos + quantidade", só que agora os elementos são os próprios
    aminoácidos, não compostos). A ORDEM final é embaralhada (a física
    de dobramento decide a forma 3D), mas a COMPOSIÇÃO exata é
    respeitada.
    """
    seq = []
    for aa, qtd in composicao.items():
        aa = aa.upper()
        if aa not in AMINOACIDOS:
            continue
        seq.extend([aa] * int(qtd))
    if len(seq) < 2:
        raise ValueError("A composição customizada precisa ter pelo menos 2 aminoácidos válidos.")
    rng = np.random.default_rng(seed)
    rng.shuffle(seq)
    return seq


def montar_sequencia_hibrida(dominios, repeticoes_linker=1):
    """Concatena os domínios intercalando o linker GS clássico entre
    eles (não é adicionado antes do primeiro nem depois do último
    domínio). Retorna a sequência completa e as "fronteiras" de cada
    domínio (índice inicial, final, índice do domínio) — usadas depois
    só para colorir a visualização 3D."""
    sequencia_final = []
    fronteiras = []
    cursor = 0
    for i, dominio in enumerate(dominios):
        sequencia_final.extend(dominio)
        fronteiras.append((cursor, cursor + len(dominio), i))
        cursor += len(dominio)
        if i < len(dominios) - 1 and repeticoes_linker > 0:
            linker = LINKER_GS * repeticoes_linker
            sequencia_final.extend(linker)
            cursor += len(linker)
    return sequencia_final, fronteiras


def gerar_proteina_hibrida(n_dominios, tamanho_dominio, k_ambiente,
                            repeticoes_linker=1, n_iter_dobra=250, seed=None,
                            dominio_customizado=None):
    """Gera e dobra uma proteína híbrida (fusion protein) composta por
    n_dominios domínios de aminoácidos unidos por linker GS, usando a
    MESMA física de dobramento (loss_cadeia, com DIST_CA_CA real) usada
    nas proteínas do resto da simulação — em vez de encadear átomos
    soltos de fármacos, que não é como fusion proteins reais existem.

    dominio_customizado: opcional, lista de códigos de 3 letras (ex.
    saída de montar_dominio_por_composicao) — se fornecido, vira o
    PRIMEIRO domínio da proteína (com composição exata escolhida pelo
    usuário); os domínios restantes (se n_dominios > 1) continuam
    sendo sorteados aleatoriamente como antes.
    """
    dominios = []
    if dominio_customizado:
        dominios.append(list(dominio_customizado))
    n_aleatorios = max(0, n_dominios - len(dominios))
    if n_aleatorios > 0:
        dominios.extend(gerar_dominios_hibridos(n_dominios=n_aleatorios, tamanho_dominio=tamanho_dominio, seed=seed))
    sequencia, fronteiras = montar_sequencia_hibrida(dominios, repeticoes_linker=repeticoes_linker)
    n = len(sequencia)

    raios = torch.tensor([raio_por_aa(a) for a in sequencia], dtype=torch.float32)
    cargas = torch.tensor([AMINOACIDOS[a]["carga"] for a in sequencia], dtype=torch.float32)
    hidro = torch.tensor([AMINOACIDOS[a]["hidrofobicidade"] for a in sequencia], dtype=torch.float32)

    gen = torch.Generator().manual_seed(seed) if seed is not None else None
    raio_inicial = max(4.0, n * 0.3)
    pontos_iniciais = gerar_pontos_enrolados(n, raio_inicial)
    ruido = (torch.rand(n, 3, generator=gen) - 0.5) * 0.5 if gen is not None else (torch.rand(n, 3) - 0.5) * 0.5
    p = (torch.tensor(pontos_iniciais, dtype=torch.float32) + ruido).clone().requires_grad_(True)

    opt = torch.optim.Adam([p], lr=0.05)
    for _ in range(n_iter_dobra):
        opt.zero_grad()
        loss_cadeia(p, raios, cargas, hidro, k_ambiente).backward()
        opt.step()

    pontos_final = p.detach().numpy()
    ligacoes = [(i, i + 1) for i in range(n - 1)]
    return {
        "pontos": pontos_final,
        "sequencia": sequencia,
        "ligacoes": ligacoes,
        "cargas": cargas.numpy(),
        "raios": raios.numpy(),
        "hidro": hidro.numpy(),
        "dominios": [f"Dominio{i+1}({len(d)}aa)" for i, d in enumerate(dominios)],
        "fronteiras": fronteiras,
    }


def criar_geometria_hibrida_visual(pontos, sequencia, fronteiras, escala_visual=0.72):
    """Mesmo estilo visual das outras cadeias (esfera por resíduo +
    cilindro de ligação), mas colorindo cada domínio com uma cor
    diferente (via PALETA_DOMINIOS) e o linker sempre em cinza —assim
    dá pra visualmente identificar onde termina um domínio e começa
    o próximo.

    escala_visual: fator só de RENDERIZAÇÃO (não afeta a física de
    dobramento, que já rodou antes desta função ser chamada) — reduz
    o tamanho das esferas/cilindros para a proteína híbrida ficar
    visualmente menor/mais compacta no visor 3D."""
    cores_residuo = [COR_LINKER] * len(sequencia)
    for (ini, fim, idx_dom) in fronteiras:
        cor = PALETA_DOMINIOS[idx_dom % len(PALETA_DOMINIOS)]
        for i in range(ini, fim):
            cores_residuo[i] = cor

    geos = criar_backbone_suave(pontos, raio=0.15 * escala_visual)
    for i, p in enumerate(pontos):
        raio = raio_por_aa(sequencia[i]) * escala_visual * 0.88
        esfera = trimesh.creation.icosphere(radius=raio, subdivisions=3)
        esfera.visual.vertex_colors = cores_residuo[i].astype(np.uint8)
        esfera.apply_translation(p)
        geos.append(esfera)
    return geos


def calcular_estatisticas_hibrida(sequencia):  # mantido p/ compatibilidade
    return calcular_estatisticas_proteina(sequencia)


# ------------------------------------------------------------------
# NOVAS FUNÇÕES — gerador simples de proteína por composição de AAs
# ------------------------------------------------------------------

def gerar_proteina_nova(
    comprimento: int = 30,
    aminoacidos_pool: list = None,
    composicao_fixa: dict = None,
    k_ambiente: dict = None,
    n_iter_dobra: int = 250,
    seed=None,
) -> dict:
    """
    Gera uma proteína nova combinando aminoácidos e dobra a cadeia com
    física de potencial (mola + VdW + Coulomb + hidrofóbico + água).

    Dois modos:
      - composicao_fixa vazio/None → sorteia `comprimento` resíduos
        aleatoriamente, usando só os AAs de `aminoacidos_pool`
        (padrão: todos os 20).
      - composicao_fixa={"ALA": 4, "LEU": 6, ...} → usa EXATAMENTE
        esses aminoácidos nessas quantidades e embaralha a ordem.
    """
    if k_ambiente is None:
        k_ambiente = ENVIRONMENTS["citoplasma"]

    rng = np.random.default_rng(seed)

    if composicao_fixa:
        seq = []
        for aa, qtd in composicao_fixa.items():
            aa = aa.upper()
            if aa in AMINOACIDOS:
                seq.extend([aa] * int(qtd))
        if len(seq) < 2:
            raise ValueError("Composição precisa ter pelo menos 2 aminoácidos válidos.")
        rng.shuffle(seq)
    else:
        pool = [aa for aa in (aminoacidos_pool or LISTA_AA) if aa in AMINOACIDOS]
        if not pool:
            pool = LISTA_AA
        comprimento = max(5, min(int(comprimento), 80))
        idx = rng.integers(0, len(pool), size=comprimento)
        seq = [pool[i] for i in idx]

    n = len(seq)
    raios  = torch.tensor([raio_por_aa(a)                        for a in seq], dtype=torch.float32)
    cargas = torch.tensor([AMINOACIDOS[a]["carga"]               for a in seq], dtype=torch.float32)
    hidro  = torch.tensor([AMINOACIDOS[a]["hidrofobicidade"]     for a in seq], dtype=torch.float32)

    raio_inicial = max(4.0, n * 0.3)
    pontos_ini = gerar_pontos_enrolados(n, raio_inicial)
    gen = torch.Generator().manual_seed(int(seed)) if seed is not None else None
    ruido = (torch.rand(n, 3, generator=gen) - 0.5) * 0.5 if gen is not None \
            else (torch.rand(n, 3) - 0.5) * 0.5
    p = (torch.tensor(pontos_ini, dtype=torch.float32) + ruido).clone().requires_grad_(True)

    opt = torch.optim.Adam([p], lr=0.05)
    for _ in range(n_iter_dobra):
        opt.zero_grad()
        loss_cadeia(p, raios, cargas, hidro, k_ambiente).backward()
        opt.step()

    with torch.no_grad():
        energia_final = float(loss_cadeia(p, raios, cargas, hidro, k_ambiente).item())
    pts = p.detach().numpy()
    return {
        "pontos":    pts,
        "sequencia": seq,
        "ligacoes":  [(i, i + 1) for i in range(n - 1)],
        "cargas":    cargas.numpy(),
        "raios":     raios.numpy(),
        "hidro":     hidro.numpy(),
        "energia_final": energia_final,
    }


def criar_geometria_proteina(pontos, sequencia, escala_visual=0.72):
    """
    Visualiza a proteína com um estilo mais próximo de um visualizador
    molecular real: um backbone (esqueleto) em tubo curvo e contínuo
    (spline suave, não segmentos retos) representando a cadeia
    principal, com as cadeias laterais dos resíduos como esferas
    coloridas por hidrofobicidade:
    - 🔴 Vermelho  → resíduo hidrofóbico (quer se esconder da água)
    - 🟢 Verde     → resíduo hidrofílico (fica na superfície, interage com água)
    """
    geos = criar_backbone_suave(pontos, raio=0.15 * escala_visual)
    for i, p in enumerate(pontos):
        raio = raio_por_aa(sequencia[i]) * escala_visual * 0.88
        esfera = trimesh.creation.icosphere(radius=raio, subdivisions=3)
        cor = COR_HIDROFOBO if aa_e_hidrofobico(sequencia[i]) else COR_HIDROFILICO
        esfera.visual.vertex_colors = cor.astype(np.uint8)
        esfera.apply_translation(p)
        geos.append(esfera)
    return geos


def calcular_estatisticas_proteina(sequencia):
    massa = sum(AMINOACIDOS[a]["massa"] for a in sequencia)
    carga = sum(AMINOACIDOS[a]["carga"] for a in sequencia)
    n_hb  = sum(1 for a in sequencia if AMINOACIDOS[a]["hidrofobicidade"] > 0)
    cont  = {}
    for a in sequencia:
        cont[a] = cont.get(a, 0) + 1
    return {
        "n_residuos":    len(sequencia),
        "massa_total_da": massa,
        "carga_liquida":  float(carga),
        "n_hidrofobo":    n_hb,
        "n_hidrofilico":  len(sequencia) - n_hb,
        "composicao":     cont,
    }


# =====================================================================
# 6. BUSCA DE DADOS REAIS NA INTERNET (NCBI / UniProt)
# =====================================================================
# Os 3 presets em VIRUS_REAIS_PRESETS acima são trechos FIXOS, escolhidos
# à mão. As duas funções abaixo permitem buscar QUALQUER outro vírus,
# bactéria (ou gene/proteína) disponível publicamente, em tempo real,
# usando as APIs REST públicas e sem chave de acesso do NCBI (ácidos
# nucleicos) e do UniProt (proteínas). Só funcionam se o servidor tiver
# acesso à internet e o pacote `requests` instalado.
_UM_PARA_TRES = {v: k for k, v in TRES_PARA_UMA.items()}


def buscar_acido_nucleico_ncbi(termo: str, max_bases: int = 90, preferir_rna: bool = True) -> dict:
    """
    Busca uma sequência de ácido nucleico REAL no NCBI Nucleotide (banco
    `nuccore`, via E-utilities, sem necessidade de chave de API) a partir
    de um termo livre — nome do organismo/vírus/bactéria, gene, ou
    accession (ex.: "SARS-CoV-2 spike", "Escherichia coli plasmid",
    "NC_045512.2") — e devolve um trecho inicial (as primeiras
    `max_bases` bases), pronto para usar como genoma customizado de um
    vírus OU nucleoide/plasmídeo de uma bactéria fictícia.

    Levanta RuntimeError com uma mensagem amigável se `requests` não
    estiver disponível, se a busca não encontrar nada, ou se a resposta
    do NCBI vier vazia/inesperada.
    """
    if requests is None:
        raise RuntimeError("O pacote 'requests' não está instalado no servidor — busca online indisponível.")
    termo = (termo or "").strip()
    if not termo:
        raise RuntimeError("Digite um termo de busca (nome do organismo, gene ou accession).")

    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    try:
        r = requests.get(f"{base}/esearch.fcgi", params={
            "db": "nuccore", "term": termo, "retmode": "json", "retmax": 1, "sort": "relevance",
        }, timeout=15)
        r.raise_for_status()
        ids = r.json().get("esearchresult", {}).get("idlist", [])
    except Exception as e:
        raise RuntimeError(f"Falha ao consultar o NCBI: {e}")
    if not ids:
        raise RuntimeError(f"Nenhum resultado encontrado no NCBI para '{termo}'.")

    seq_id = ids[0]
    try:
        r2 = requests.get(f"{base}/efetch.fcgi", params={
            "db": "nuccore", "id": seq_id, "rettype": "fasta", "retmode": "text",
        }, timeout=15)
        r2.raise_for_status()
    except Exception as e:
        raise RuntimeError(f"Falha ao baixar a sequência do NCBI: {e}")

    linhas = r2.text.strip().splitlines()
    if not linhas or not linhas[0].startswith(">"):
        raise RuntimeError("Resposta inesperada do NCBI (sem cabeçalho FASTA).")
    cabecalho = linhas[0][1:].strip()
    bases_brutas = "".join(linhas[1:]).upper()
    bases_dna = re.sub(r"[^ACGT]", "", bases_brutas)
    if not bases_dna:
        raise RuntimeError("A sequência retornada pelo NCBI não contém bases DNA/RNA reconhecíveis.")

    trecho_dna = bases_dna[:max(6, int(max_bases))]
    if preferir_rna:
        trecho_final, tipo_final = trecho_dna.replace("T", "U"), "RNA"
    else:
        trecho_final, tipo_final = trecho_dna, "DNA"
    return {
        "nome": cabecalho[:90],
        "tipo_acido": tipo_final,
        "sequencia_acido": trecho_final,
        "fonte": f"NCBI Nucleotide, accession {seq_id} — primeiras {len(trecho_final)} de {len(bases_dna)} bases",
        "accession": seq_id,
        "total_bases_disponiveis": len(bases_dna),
    }


def buscar_proteina_uniprot(termo: str, max_residuos: int = 40) -> dict:
    """
    Busca uma proteína REAL no UniProt (API REST pública, sem chave) a
    partir de um termo livre — nome do gene/proteína/organismo (ex.:
    "human insulin", "spike protein SARS-CoV-2", "lysozyme") — e devolve
    um trecho inicial da sequência de aminoácidos já convertido para os
    códigos de 3 letras usados internamente por este simulador, mais uma
    composição pronta (contagem por aminoácido) para preencher o modo
    "composição exata" da geração de proteína.
    """
    if requests is None:
        raise RuntimeError("O pacote 'requests' não está instalado no servidor — busca online indisponível.")
    termo = (termo or "").strip()
    if not termo:
        raise RuntimeError("Digite um termo de busca (nome do gene, proteína ou organismo).")

    try:
        r = requests.get("https://rest.uniprot.org/uniprotkb/search", params={
            "query": termo, "fields": "accession,id,protein_name,organism_name,sequence",
            "format": "json", "size": 1,
        }, timeout=15)
        r.raise_for_status()
        resultados = r.json().get("results", [])
    except Exception as e:
        raise RuntimeError(f"Falha ao consultar o UniProt: {e}")
    if not resultados:
        raise RuntimeError(f"Nenhum resultado encontrado no UniProt para '{termo}'.")

    entrada = resultados[0]
    seq_completa = (entrada.get("sequence") or {}).get("value", "")
    if not seq_completa:
        raise RuntimeError("O UniProt não retornou uma sequência para esse resultado.")

    trecho = seq_completa[:max(6, int(max_residuos))]
    seq_3letras = [_UM_PARA_TRES[c] for c in trecho if c in _UM_PARA_TRES]
    if len(seq_3letras) < 2:
        raise RuntimeError("A sequência encontrada não pôde ser convertida (só aminoácidos não-padrão).")

    composicao = {}
    for aa in seq_3letras:
        composicao[aa] = composicao.get(aa, 0) + 1

    nome_proteina = (((entrada.get("proteinDescription") or {}).get("recommendedName") or {})
                      .get("fullName") or {}).get("value")
    nome_organismo = (entrada.get("organism") or {}).get("scientificName")
    accession = entrada.get("primaryAccession", "")
    nome_final = nome_proteina or entrada.get("uniProtkbId") or accession or termo

    return {
        "nome": nome_final,
        "organismo": nome_organismo,
        "accession": accession,
        "sequencia_fasta": "".join(TRES_PARA_UMA[a] for a in seq_3letras),
        "composicao_fixa": composicao,
        "fonte": f"UniProt {accession} — primeiros {len(seq_3letras)} de {len(seq_completa)} resíduos",
        "total_residuos_disponiveis": len(seq_completa),
    }


def buscar_acido_nucleico_ncbi_lista(termo: str, max_resultados: int = 8) -> list:
    """
    Busca no NCBI Nucleotide e devolve uma LISTA de candidatos (título,
    organismo, tamanho em bases, accession) para o usuário ESCOLHER qual
    usar — em vez de baixar direto o primeiro resultado como
    `buscar_acido_nucleico_ncbi` faz. Depois que o usuário escolher um
    item da lista, chame `obter_acido_nucleico_por_accession` com o
    accession escolhido para baixar a sequência completa.

    Usa dois passos de E-utilities do NCBI: esearch (pega os IDs que
    batem com o termo) + esummary (pega título/organismo/tamanho de
    cada ID em uma única chamada, SEM baixar a sequência inteira de
    cada um — isso mantém a busca rápida mesmo com vários candidatos).
    """
    if requests is None:
        raise RuntimeError("O pacote 'requests' não está instalado no servidor — busca online indisponível.")
    termo = (termo or "").strip()
    if not termo:
        raise RuntimeError("Digite um termo de busca (nome do organismo, gene ou accession).")

    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    try:
        r = requests.get(f"{base}/esearch.fcgi", params={
            "db": "nuccore", "term": termo, "retmode": "json",
            "retmax": max(1, min(20, int(max_resultados))), "sort": "relevance",
        }, timeout=15)
        r.raise_for_status()
        ids = r.json().get("esearchresult", {}).get("idlist", [])
    except Exception as e:
        raise RuntimeError(f"Falha ao consultar o NCBI: {e}")
    if not ids:
        raise RuntimeError(f"Nenhum resultado encontrado no NCBI para '{termo}'.")

    try:
        r2 = requests.get(f"{base}/esummary.fcgi", params={
            "db": "nuccore", "id": ",".join(ids), "retmode": "json",
        }, timeout=15)
        r2.raise_for_status()
        resumo = r2.json().get("result", {})
    except Exception as e:
        raise RuntimeError(f"Falha ao consultar resumos no NCBI: {e}")

    candidatos = []
    for seq_id in ids:
        info = resumo.get(seq_id) or {}
        if not info:
            continue
        tamanho = info.get("slen")
        candidatos.append({
            "accession": info.get("accessionversion") or seq_id,
            "titulo": (info.get("title") or "(sem título)")[:140],
            "organismo": info.get("organism", ""),
            "tamanho_bases": int(tamanho) if tamanho else None,
        })
    if not candidatos:
        raise RuntimeError(f"Nenhum resultado utilizável encontrado no NCBI para '{termo}'.")
    return candidatos


def obter_acido_nucleico_por_accession(accession: str, max_bases: int = 90, preferir_rna: bool = True) -> dict:
    """
    Baixa a sequência (FASTA) de UM accession específico já escolhido
    pelo usuário na lista de `buscar_acido_nucleico_ncbi_lista`, e
    devolve no mesmo formato que `buscar_acido_nucleico_ncbi` (pronta
    para usar como genoma customizado de um vírus/bactéria).
    """
    if requests is None:
        raise RuntimeError("O pacote 'requests' não está instalado no servidor — busca online indisponível.")
    accession = (accession or "").strip()
    if not accession:
        raise RuntimeError("Accession não informado.")

    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    try:
        r = requests.get(f"{base}/efetch.fcgi", params={
            "db": "nuccore", "id": accession, "rettype": "fasta", "retmode": "text",
        }, timeout=15)
        r.raise_for_status()
    except Exception as e:
        raise RuntimeError(f"Falha ao baixar a sequência do NCBI: {e}")

    linhas = r.text.strip().splitlines()
    if not linhas or not linhas[0].startswith(">"):
        raise RuntimeError("Resposta inesperada do NCBI (sem cabeçalho FASTA).")
    cabecalho = linhas[0][1:].strip()
    bases_brutas = "".join(linhas[1:]).upper()
    bases_dna = re.sub(r"[^ACGT]", "", bases_brutas)
    if not bases_dna:
        raise RuntimeError("A sequência retornada pelo NCBI não contém bases DNA/RNA reconhecíveis.")

    trecho_dna = bases_dna[:max(6, int(max_bases))]
    if preferir_rna:
        trecho_final, tipo_final = trecho_dna.replace("T", "U"), "RNA"
    else:
        trecho_final, tipo_final = trecho_dna, "DNA"
    return {
        "nome": cabecalho[:90],
        "tipo_acido": tipo_final,
        "sequencia_acido": trecho_final,
        "fonte": f"NCBI Nucleotide, accession {accession} — primeiras {len(trecho_final)} de {len(bases_dna)} bases",
        "accession": accession,
        "total_bases_disponiveis": len(bases_dna),
    }


def buscar_proteina_uniprot_lista(termo: str, max_resultados: int = 8) -> list:
    """
    Busca no UniProt e devolve uma LISTA de candidatos (nome, organismo,
    accession, tamanho em resíduos) para o usuário ESCOLHER qual usar —
    em vez de baixar direto o primeiro resultado como
    `buscar_proteina_uniprot` faz. Depois que o usuário escolher, chame
    `obter_proteina_por_accession` com o accession escolhido.
    """
    if requests is None:
        raise RuntimeError("O pacote 'requests' não está instalado no servidor — busca online indisponível.")
    termo = (termo or "").strip()
    if not termo:
        raise RuntimeError("Digite um termo de busca (nome do gene, proteína ou organismo).")

    try:
        r = requests.get("https://rest.uniprot.org/uniprotkb/search", params={
            "query": termo, "fields": "accession,id,protein_name,organism_name,length",
            "format": "json", "size": max(1, min(20, int(max_resultados))),
        }, timeout=15)
        r.raise_for_status()
        resultados = r.json().get("results", [])
    except Exception as e:
        raise RuntimeError(f"Falha ao consultar o UniProt: {e}")
    if not resultados:
        raise RuntimeError(f"Nenhum resultado encontrado no UniProt para '{termo}'.")

    candidatos = []
    for entrada in resultados:
        nome_proteina = (((entrada.get("proteinDescription") or {}).get("recommendedName") or {})
                          .get("fullName") or {}).get("value")
        nome_organismo = (entrada.get("organism") or {}).get("scientificName")
        accession = entrada.get("primaryAccession", "")
        seq_info = entrada.get("sequence") or {}
        candidatos.append({
            "accession": accession,
            "nome": nome_proteina or entrada.get("uniProtkbId") or accession,
            "organismo": nome_organismo,
            "tamanho_residuos": seq_info.get("length") if isinstance(seq_info, dict) else None,
        })
    return candidatos


def obter_proteina_por_accession(accession: str, max_residuos: int = 40) -> dict:
    """
    Baixa a sequência completa de UMA proteína (por accession do
    UniProt) já escolhida pelo usuário na lista de
    `buscar_proteina_uniprot_lista`, e devolve no mesmo formato que
    `buscar_proteina_uniprot`.
    """
    if requests is None:
        raise RuntimeError("O pacote 'requests' não está instalado no servidor — busca online indisponível.")
    accession = (accession or "").strip()
    if not accession:
        raise RuntimeError("Accession não informado.")

    try:
        r = requests.get(f"https://rest.uniprot.org/uniprotkb/{accession}.json", timeout=15)
        r.raise_for_status()
        entrada = r.json()
    except Exception as e:
        raise RuntimeError(f"Falha ao consultar o UniProt: {e}")

    seq_completa = (entrada.get("sequence") or {}).get("value", "")
    if not seq_completa:
        raise RuntimeError("O UniProt não retornou uma sequência para esse resultado.")

    trecho = seq_completa[:max(6, int(max_residuos))]
    seq_3letras = [_UM_PARA_TRES[c] for c in trecho if c in _UM_PARA_TRES]
    if len(seq_3letras) < 2:
        raise RuntimeError("A sequência encontrada não pôde ser convertida (só aminoácidos não-padrão).")

    composicao = {}
    for aa in seq_3letras:
        composicao[aa] = composicao.get(aa, 0) + 1

    nome_proteina = (((entrada.get("proteinDescription") or {}).get("recommendedName") or {})
                      .get("fullName") or {}).get("value")
    nome_organismo = (entrada.get("organism") or {}).get("scientificName")
    nome_final = nome_proteina or entrada.get("uniProtkbId") or accession

    return {
        "nome": nome_final,
        "organismo": nome_organismo,
        "accession": accession,
        "sequencia_fasta": "".join(TRES_PARA_UMA[a] for a in seq_3letras),
        "composicao_fixa": composicao,
        "fonte": f"UniProt {accession} — primeiros {len(seq_3letras)} de {len(seq_completa)} resíduos",
        "total_residuos_disponiveis": len(seq_completa),
    }


# =====================================================================
# 6B. IMAGEM REAL DE REFERÊNCIA (Wikimedia Commons) — NÃO reconstrói
# geometria 3D a partir da imagem (isso exigiria um pipeline de visão
# computacional bem mais pesado, fora do escopo daqui). Em vez disso,
# usa uma foto/micrografia REAL do organismo só para ler dois números
# simples — a cor média e a proporção largura/altura — e usa esses
# números para "puxar" a paleta de cor e o alongamento do corpo do
# modelo procedural na direção do organismo real. Sempre best-effort:
# se não achar imagem, ou faltar `requests`/Pillow, ou o Commons estiver
# fora do ar, o modelo simplesmente segue 100% fictício (sem travar).
# =====================================================================
# ── Banco local de imagens de referência (cache em disco) ────────────
# Guarda, por organismo já consultado, a cor média, a proporção
# largura/altura e a URL pública da foto — assim, se o mesmo vírus ou
# bactéria real for usado de novo (em outra simulação, ou pelo mesmo
# usuário testando parâmetros diferentes), o Wikimedia Commons não
# precisa ser consultado de novo. É um cache best-effort: se não der
# pra ler/gravar em disco por qualquer motivo, a simulação simplesmente
# segue sem cache (busca direto na internet toda vez).
CACHE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "dados_cache")
CACHE_IMAGENS_PATH = os.path.join(CACHE_DIR, "banco_imagens_referencia.json")


def _carregar_banco_imagens() -> dict:
    try:
        with open(CACHE_IMAGENS_PATH, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {}


def _salvar_banco_imagens(banco: dict) -> None:
    try:
        os.makedirs(CACHE_DIR, exist_ok=True)
        tmp_path = CACHE_IMAGENS_PATH + ".tmp"
        with open(tmp_path, "w", encoding="utf-8") as f:
            json.dump(banco, f, ensure_ascii=False, indent=2)
        os.replace(tmp_path, CACHE_IMAGENS_PATH)
    except Exception:
        pass


def _texto_sem_html(valor) -> str:
    """Converte metadados do Wikimedia em texto simples para filtragem."""
    if isinstance(valor, dict):
        valor = valor.get("value", "")
    valor = str(valor or "")
    valor = re.sub(r"<[^>]+>", " ", valor)
    return re.sub(r"\s+", " ", valor).strip()


def _normalizar_busca(valor: str) -> str:
    valor = _texto_sem_html(valor).lower()
    valor = unicodedata.normalize("NFKD", valor).encode("ascii", "ignore").decode("ascii")
    return re.sub(r"\s+", " ", valor).strip()


def _meta_imagem(extmetadata: dict, chave: str) -> str:
    return _texto_sem_html((extmetadata or {}).get(chave, ""))


def _inferir_forma_bacteria(texto: str, proporcao_lxa: float | None = None) -> str | None:
    """Infere apenas uma forma procedural simples a partir do texto científico.

    Não tenta afirmar que uma fotografia 2D fornece uma reconstrução 3D.
    """
    t = _normalizar_busca(texto)
    if re.search(r"\b(coccus|cocci|coco|cocos|spherical|sphere|esferic)", t):
        return "coco"
    if re.search(r"\b(bacillus|bacilli|bacilo|bacilos|rod[- ]shaped|rod bacteria|bastonete)", t):
        return "bacilo"
    if re.search(r"\b(spirillum|spirochete|spiral|helical|espiral)", t):
        # O renderer atual não possui uma bactéria helicoidal dedicada;
        # usamos o corpo alongado como aproximação segura.
        return "bacilo"
    return None


def _pontuar_candidato_imagem(titulo: str, texto: str, termo: str, tipo: str) -> tuple[int, str]:
    """Pontua uma imagem científica e rejeita arte/ilustrações."""
    combinado = _normalizar_busca(f"{titulo} {texto}")
    organismo = _normalizar_busca(termo)
    negativos = [
        "painting", "paintings", "oil on canvas", "watercolor", "watercolour",
        "illustration", "drawing", "sketch", "artwork", "digital art", "poster",
        "logo", "cartoon", "comic", "sculpture", "statue", "museum", "toy",
        "figurine", "diagram", "schematic", "map", "flag", "album cover",
        "book cover", "stamp", "3d render", "rendering", "computer generated",
    ]
    if any(palavra in combinado for palavra in negativos):
        return -1000, "arte/ilustração rejeitada"

    score = 0
    palavras_organismo = [p for p in re.findall(r"[a-z0-9-]{4,}", organismo)]
    score += min(30, 8 * sum(1 for p in palavras_organismo if p in combinado))

    comuns = [
        "micrograph", "micrograf", "microscopy", "microscop", "electron",
        "sem", "tem", "fluorescen", "phase contrast", "gram stain",
        "stained", "culture", "colony", "cell", "cells", "pathogen",
        "laboratory", "laboratory", "scientific",
    ]
    score += 5 * sum(1 for palavra in comuns if palavra in combinado)

    if tipo == "bacteria":
        score += 10 * sum(1 for palavra in ("bacter", "bacill", "coccus", "cocci", "gram") if palavra in combinado)
    elif tipo == "virus":
        score += 10 * sum(1 for palavra in ("virus", "virion", "viral", "phage", "electron micrograph") if palavra in combinado)

    # Sem nenhum indício de organismo ou microscopia, não aceitamos uma
    # imagem genérica que possa ser uma pintura ou uma ilustração.
    if score < 8:
        return -100, "sem indícios científicos suficientes"
    return score, "candidato científico"


def obter_referencia_visual(termo: str, tipo: str = "organismo", log=None) -> dict | None:
    """Obtém uma imagem científica filtrada de vírus ou bactéria.

    A chave do cache inclui o tipo e uma versão implícita do novo filtro;
    assim, imagens antigas escolhidas pela busca permissiva não são
    reutilizadas.
    """
    def _log(msg):
        if log:
            log(msg)

    tipo = "bacteria" if str(tipo).lower().startswith("bact") else "virus" if str(tipo).lower().startswith("vir") else "organismo"
    termo_limpo = re.sub(r"\s+", " ", (termo or "").strip()).strip('"')
    if not termo_limpo:
        return None

    chave = f"v2:{tipo}:{termo_limpo.lower()}"
    banco = _carregar_banco_imagens()
    if chave in banco:
        _log(f"  🗄️ Referência científica de '{termo_limpo}' já está no banco local.")
        return banco[chave]

    dados_img = buscar_imagem_referencia_wikimedia(termo_limpo, tipo=tipo, log=log)
    if not dados_img:
        return None
    cor_media = analisar_imagem_referencia(dados_img["thumb_bytes"])
    if cor_media is None:
        _log("  ⚠ Não foi possível analisar a imagem científica; seguindo sem ajuste visual.")
        return None

    proporcao = dados_img["largura_original"] / max(dados_img["altura_original"], 1)
    texto_classificacao = f"{dados_img.get('titulo', '')} {dados_img.get('descricao', '')} {dados_img.get('categorias', '')}"
    forma_sugerida = _inferir_forma_bacteria(texto_classificacao, proporcao) if tipo == "bacteria" else None
    entrada = {
        "tipo": tipo,
        "termo_consultado": termo_limpo,
        "cor_media": [round(float(c), 1) for c in cor_media],
        "proporcao": round(float(proporcao), 3),
        "forma_sugerida": forma_sugerida,
        "thumb_url": dados_img.get("thumb_url_publica") or dados_img["fonte_url"],
        "fonte_url": dados_img["fonte_url"],
        "titulo": dados_img["titulo"],
        "autor": dados_img.get("autor", ""),
        "licenca": dados_img.get("licenca", ""),
        "consulta": dados_img.get("consulta", ""),
    }
    banco[chave] = entrada
    _salvar_banco_imagens(banco)
    _log(
        f"  🖼️ Referência científica selecionada para '{termo_limpo}': "
        f"{dados_img['titulo']} ({dados_img['fonte_url']})."
    )
    return entrada


def buscar_imagem_referencia_wikimedia(termo: str, tipo: str = "organismo", log=None) -> dict | None:
    """Busca uma micrografia científica no Wikimedia Commons.

    A versão anterior pegava o primeiro resultado textual e podia escolher
    pinturas. Agora coleta vários candidatos, lê título/categorias/descrição,
    rejeita arte e só baixa uma imagem com sinais de microscopia ou de
    material biológico. Se não encontrar uma candidata confiável, retorna
    None em vez de escolher uma imagem aleatória.
    """
    def _log(msg):
        if log:
            log(msg)

    if requests is None:
        _log("  ⚠ 'requests' não instalado — pulando referência científica.")
        return None
    termo = re.sub(r"\s+", " ", (termo or "").strip()).strip('"')
    if not termo:
        return None

    tipo = "bacteria" if str(tipo).lower().startswith("bact") else "virus" if str(tipo).lower().startswith("vir") else "organismo"
    tipo_label = "bacterial" if tipo == "bacteria" else "viral" if tipo == "virus" else "biological"
    termo_busca = termo.replace('"', "")[:120]
    if tipo == "bacteria":
        consultas = [
            f'"{termo_busca}" bacteria micrograph',
            f'"{termo_busca}" bacterial microscopy',
            f'"{termo_busca}" electron micrograph',
        ]
    elif tipo == "virus":
        consultas = [
            f'"{termo_busca}" virus micrograph',
            f'"{termo_busca}" virion electron microscopy',
            f'"{termo_busca}" viral microscopy',
        ]
    else:
        consultas = [f'"{termo_busca}" {tipo_label} micrograph']

    base = "https://commons.wikimedia.org/w/api.php"
    headers = {"User-Agent": "SimuladorMolecular/1.1 (scientific-reference-filter)"}
    titulos = []
    consulta_por_titulo = {}
    try:
        for consulta in consultas:
            r = requests.get(base, params={
                "action": "query", "list": "search", "srsearch": consulta,
                "srnamespace": 6, "srwhat": "text", "format": "json",
                "srlimit": 12,
            }, timeout=15, headers=headers)
            r.raise_for_status()
            hits = r.json().get("query", {}).get("search", [])
            for hit in hits:
                titulo = hit.get("title", "")
                if titulo and titulo not in consulta_por_titulo:
                    titulos.append(titulo)
                    consulta_por_titulo[titulo] = consulta
    except Exception as e:
        _log(f"  ⚠ Falha ao consultar o Wikimedia Commons: {e}")
        return None

    if not titulos:
        _log(f"  (nenhuma imagem científica encontrada para '{termo}')")
        return None

    candidatos = []
    try:
        for inicio in range(0, len(titulos), 40):
            lote = titulos[inicio:inicio + 40]
            r2 = requests.get(base, params={
                "action": "query", "prop": "imageinfo|categories",
                "titles": "|".join(lote),
                "iiprop": "url|size|mime|extmetadata", "iiurlwidth": 480,
                "cllimit": "max", "format": "json", "formatversion": 2,
            }, timeout=15, headers=headers)
            r2.raise_for_status()
            paginas = r2.json().get("query", {}).get("pages", [])
            if isinstance(paginas, dict):
                paginas = list(paginas.values())
            for pagina in paginas:
                infos = pagina.get("imageinfo") or []
                if not infos:
                    continue
                info = infos[0]
                mime = str(info.get("mime", "")).lower()
                if not mime.startswith("image/") or mime in {"image/svg+xml", "image/gif"}:
                    continue
                largura = int(info.get("width") or 0)
                altura = int(info.get("height") or 0)
                if largura < 100 or altura < 100:
                    continue
                ext = info.get("extmetadata") or {}
                descricao = " ".join(_meta_imagem(ext, chave) for chave in (
                    "ImageDescription", "ObjectName", "Categories", "Credit", "Artist"
                ))
                categorias = " ".join(
                    _texto_sem_html(c.get("title", "")) for c in (pagina.get("categories") or [])
                )
                titulo = pagina.get("title", "").replace("File:", "", 1).strip()
                texto = f"{descricao} {categorias}"
                score, motivo = _pontuar_candidato_imagem(titulo, texto, termo, tipo)
                if score < 0:
                    continue
                candidatos.append({
                    "score": score,
                    "titulo": titulo,
                    "info": info,
                    "descricao": descricao,
                    "categorias": categorias,
                    "consulta": consulta_por_titulo.get(pagina.get("title", ""), consultas[0]),
                    "autor": _meta_imagem(ext, "Artist") or _meta_imagem(ext, "Credit"),
                    "licenca": _meta_imagem(ext, "LicenseShortName") or _meta_imagem(ext, "UsageTerms"),
                })
    except Exception as e:
        _log(f"  ⚠ Falha ao ler metadados do Wikimedia Commons: {e}")
        return None

    if not candidatos:
        _log(
            f"  ⚠ Nenhuma micrografia confiável foi encontrada para '{termo}'. "
            "Nenhuma pintura ou ilustração será usada."
        )
        return None

    candidatos.sort(key=lambda item: item["score"], reverse=True)
    for candidato in candidatos[:8]:
        info = candidato["info"]
        thumb_url = info.get("thumburl") or info.get("url")
        if not thumb_url:
            continue
        try:
            r3 = requests.get(thumb_url, timeout=15, headers=headers)
            r3.raise_for_status()
            # Valida que o conteúdo é uma imagem que Pillow consegue abrir;
            # isso evita guardar respostas HTML ou arquivos corrompidos.
            if PIL_DISPONIVEL:
                with Image.open(io.BytesIO(r3.content)) as imagem:
                    imagem.verify()
            return {
                "thumb_bytes": r3.content,
                "thumb_url_publica": thumb_url,
                "largura_original": int(info.get("width") or 1),
                "altura_original": int(info.get("height") or 1),
                "fonte_url": "https://commons.wikimedia.org/wiki/" + quote(
                    candidato["titulo"].replace(" ", "_"), safe="/_:(),.-"
                ),
                "titulo": candidato["titulo"],
                "descricao": candidato["descricao"],
                "categorias": candidato["categorias"],
                "autor": candidato["autor"],
                "licenca": candidato["licenca"],
                "consulta": candidato["consulta"],
            }
        except Exception:
            continue

    _log(f"  ⚠ As imagens encontradas para '{termo}' não puderam ser validadas.")
    return None


def analisar_imagem_referencia(thumb_bytes: bytes) -> np.ndarray | None:
    """Abre a imagem baixada e devolve sua cor MÉDIA como [R,G,B]
    (float, 0-255) — uma leitura simples e barata, o suficiente pra
    'puxar' a paleta do modelo procedural na direção certa (ex.: uma
    bactéria com foto meio azulada puxa a parede pro azul), sem
    pretender ser uma análise de imagem sofisticada."""
    if not PIL_DISPONIVEL:
        return None
    try:
        img = Image.open(io.BytesIO(thumb_bytes)).convert("RGB")
        img = img.resize((48, 48))
        arr = np.asarray(img, dtype=float).reshape(-1, 3)
        return arr.mean(axis=0)
    except Exception:
        return None


def _blend_cor(cor_original: np.ndarray, cor_imagem: np.ndarray, peso_imagem: float = 0.5) -> np.ndarray:
    """Mistura a cor RGB original (uint8, com canal alfa opcional) com a
    cor extraída da imagem real, preservando o canal alfa (opacidade)
    do original — só a cor muda, a transparência da parede não."""
    cor_original = np.array(cor_original, dtype=float)
    nova = cor_original.copy()
    nova[:3] = cor_original[:3] * (1 - peso_imagem) + np.clip(cor_imagem, 0, 255) * peso_imagem
    return np.clip(nova, 0, 255).astype(np.uint8)


def aplicar_referencia_no_virus(estilo: dict, cor_imagem: np.ndarray, proporcao_lxa: float) -> dict:
    """Ajusta a paleta e o alongamento do estilo de UM vírus com base
    numa imagem real: a cor do envelope é puxada em direção à cor média
    da foto; a proporção largura/altura da foto (um proxy simples e
    imperfeito da forma real — não uma medição precisa) empurra o
    vírus para uma forma mais alongada ou mais esférica dentro de
    limites razoáveis, pra nunca deformar demais o modelo."""
    novo = dict(estilo)
    novo["cor_envelope"] = _blend_cor(estilo["cor_envelope"], cor_imagem, peso_imagem=0.45)
    if proporcao_lxa > 1.2:
        fator = float(np.clip(1.0 + (proporcao_lxa - 1.0) * 0.6, 1.0, 1.9))
        novo["forma_alongada"] = True
        novo["fator_alongamento"] = fator
    elif proporcao_lxa < 0.85:
        novo["forma_alongada"] = False
        novo["fator_alongamento"] = 1.0
    return novo


def aplicar_referencia_na_bacteria(
    estilo: dict,
    cor_imagem: np.ndarray,
    proporcao_lxa: float,
    forma_sugerida: str | None = None,
) -> dict:
    """Igual a `aplicar_referencia_no_virus`, mas para uma bactéria: a
    cor da parede E do citoplasma são puxadas em direção à cor média da
    foto real, e — só quando a forma é 'bacilo' (tem um eixo alongado
    de verdade) — a proporção largura/altura da foto ajusta a razão
    comprimento/raio dentro de limites que preservam uma forma de
    bastonete plausível."""
    novo = dict(estilo)
    novo["cor_parede"] = _blend_cor(estilo["cor_parede"], cor_imagem, peso_imagem=0.45)
    novo["cor_citoplasma"] = _blend_cor(estilo["cor_citoplasma"], cor_imagem, peso_imagem=0.35)
    if forma_sugerida in {"bacilo", "coco"}:
        novo["forma"] = forma_sugerida
    if novo.get("forma") == "bacilo" and novo.get("raio_corpo"):
        # A proporção da foto é apenas um ajuste de referência; não é uma
        # reconstrução 3D. Mantemos limites biologicamente plausíveis.
        razao_alvo = float(np.clip(proporcao_lxa, 1.8, 6.0))
        novo["comprimento"] = round(max(novo["raio_corpo"] * 2.2,
                                         novo["raio_corpo"] * razao_alvo), 2)
    elif novo.get("forma") == "coco":
        novo["comprimento"] = 0.0
    return novo


def construir_resposta_simulada(tipo_resposta, alvos_runtime, resultados_triagem,
                                limite_energetico=None):
    """Converte o ranking de docking em uma resposta legível para a UI.

    A resposta é uma classificação computacional/educacional. Ela não
    estima dose, segurança, resistência ou eficácia clínica. Quando o
    usuário informa um limite, considera-se que uma energia menor ou
    igual ao limite atende ao critério escolhido.
    """
    info = TIPOS_RESPOSTA.get(tipo_resposta, TIPOS_RESPOSTA["nenhuma"])
    candidatos_por_alvo = {}
    for item in resultados_triagem or []:
        candidatos_por_alvo.setdefault(item["virus"], []).append(item)

    alvos = []
    n_dentro = 0
    for alvo in alvos_runtime or []:
        nome = alvo["nome"]
        candidatos = candidatos_por_alvo.get(nome, [])
        melhor = min(candidatos, key=lambda x: float(x["energia"])) if candidatos else None
        if melhor is None:
            alvos.append({
                "nome": nome,
                "tipo": alvo.get("tipo_alvo", "alvo"),
                "energia_melhor": None,
                "score_local": None,
                "dentro_limite": False if limite_energetico is not None else None,
                "status": "sem resultado de docking",
                "ligante_melhor": None,
            })
            continue

        energias = [float(x["energia"]) for x in candidatos]
        energia = float(melhor["energia"])
        minimo, maximo = min(energias), max(energias)
        score_local = 50.0 if abs(maximo - minimo) < 1e-9 else 100.0 * (maximo - energia) / (maximo - minimo)
        if limite_energetico is None:
            dentro = None
            status = "melhor candidato ranqueado"
        else:
            dentro = energia <= limite_energetico
            if dentro:
                n_dentro += 1
                status = "dentro do limite energético"
            else:
                status = "acima do limite energético"
        alvos.append({
            "nome": nome,
            "tipo": alvo.get("tipo_alvo", "alvo"),
            "energia_melhor": round(energia, 4),
            "score_local": round(float(score_local), 1),
            "dentro_limite": dentro,
            "status": status,
            "ligante_melhor": melhor.get("composto"),
        })

    total = len(alvos)
    cobertura = (n_dentro / total * 100.0) if limite_energetico is not None and total else None
    if limite_energetico is None:
        resumo = "Ranking ilustrativo concluído; nenhum corte energético foi aplicado."
    elif total:
        resumo = f"{n_dentro} de {total} alvo(s) dentro do limite (cobertura ilustrativa de {cobertura:.1f}%)."
    else:
        resumo = "Nenhum alvo foi selecionado para avaliar a resposta."

    return {
        "tipo": tipo_resposta if tipo_resposta in TIPOS_RESPOSTA else "nenhuma",
        "label": info["label"],
        "descricao": info["descricao"],
        "familia": info["familia"],
        "limite_energetico": round(float(limite_energetico), 4) if limite_energetico is not None else None,
        "regra": "energia menor ou igual ao limite" if limite_energetico is not None else "sem corte energético",
        "resumo": resumo,
        "cobertura_percentual": round(cobertura, 1) if cobertura is not None else None,
        "alvos": alvos,
        "observacao": (
            "Modelo físico/procedural ilustrativo: a categoria escolhida organiza a leitura "
            "do docking e vale para vírus e bactérias, mas não representa uma recomendação "
            "terapêutica nem comprova atividade antiviral ou antibacteriana."
        ),
    }


# =====================================================================
# 7. FUNÇÃO PRINCIPAL — chamada pelo servidor web ou por rodar_local.py
# =====================================================================
def executar_simulacao(config: dict, output_dir: str, progresso=None) -> dict:
    """
    config esperado (todos os campos são opcionais, com defaults):
      {
        "n_dominios": 2,             # nº de domínios de aminoácidos na proteína híbrida
        "tamanho_dominio_min": 8,    # tamanho mínimo (em resíduos) de cada domínio
        "tamanho_dominio_max": 20,   # tamanho máximo (em resíduos) de cada domínio
        "repeticoes_linker": 1,      # quantas vezes repetir o motivo GGGGS entre domínios
        "n_hibridos": 1,             # quantas proteínas híbridas gerar
        "ambiente": "citoplasma",    # chave de ENVIRONMENTS
        "modo_rapido": True,         # reduz iterações/resíduos p/ resposta web rápida
        "limite_energetico": None,   # energia máxima aceita para a formação (unidade do modelo)
        "tipo_resposta": "nenhuma",  # estratégia didática contra vírus/bactérias
        "elementos": ["Paracetamol", "Cafeina", "Cloroquina"],  # ligantes p/ triagem
        "incluir_compostos_reais": False,  # roda triagem também com os compostos reais como ligantes
        "gerar_virus": True,         # gera os 5 modelos virais fictícios
        "seed": None,

        # ---- NOVO: aminoácidos escolhidos pelo usuário (em vez de sorteados) ----
        "aminoacidos_customizados": {},   # ex.: {"ALA": 4, "LEU": 6, "LYS": 2}
                                            # vira o 1º domínio da proteína híbrida,
                                            # com essa composição EXATA (ordem embaralhada)

        # ---- NOVO: bactérias fictícias ----
        "gerar_bacterias": False,
        "n_bacterias": 2,   # 1 a 6 — cada uma gerada proceduralmente e diferente das outras

        # ---- NOVO: vírus com genoma NÃO pré-definido (digitado pelo usuário) ----
        "virus_customizado": None,
        # formato: {
        #   "nome": "MeuVirusX", "tipo_acido": "RNA" | "DNA",
        #   "sequencia_acido": "AUGCAUGC...",   # texto cru digitado pelo usuário
        #   "raio_envelope": 6.0, "n_espiculas": 30,   # opcionais, tem defaults
        # }
      }

    progresso: função opcional callback(str) para relatar etapas (usada pelo
    servidor web para atualizar o status do job).
    """
    def log(msg):
        print(msg)
        if progresso:
            progresso(msg)

    cfg = {
        # --- proteína nova ---
        "comprimento_proteina": 30,      # nº total de resíduos da cadeia
        "aminoacidos_pool": [],          # lista de AAs permitidos ([] = todos os 20)
        "composicao_fixa": {},           # {"ALA":4,"LEU":6,...} ignora comprimento/pool
        "n_proteinas": 1,                # quantas proteínas gerar
        # --- ambiente e velocidade ---
        "ambiente": "citoplasma",
        "modo_rapido": True,
        "limite_energetico": None,
        "tipo_resposta": "nenhuma",
        "seed": None,
        # --- ligantes reais ---
        "elementos": ["Paracetamol", "Cafeina", "Cloroquina"],
        "incluir_compostos_reais": False,
        # --- alvos ---
        "gerar_virus": True,
        "gerar_bacterias": False,
        "n_bacterias": 2,
        "virus_customizado": None,
        # --- NOVO: previsão de mutações nos vírus ---
        "prever_mutacao": False,   # gera variante(s) mutante(s) dos vírus fictícios
        "n_mutacoes": 3,           # nº de mutações pontuais sorteadas por vírus (1-10)
    }
    cfg.update({k: v for k, v in (config or {}).items() if v is not None})

    comprimento_prot = max(5, min(int(cfg["comprimento_proteina"]), 80))
    aa_pool   = [aa for aa in cfg.get("aminoacidos_pool", []) if aa in AMINOACIDOS]
    comp_fixa = {aa.upper(): int(q)
                 for aa, q in (cfg.get("composicao_fixa") or {}).items()
                 if aa.upper() in AMINOACIDOS}
    n_proteinas = max(1, min(int(cfg.get("n_proteinas", 1)), 5))
    ambiente_key = cfg["ambiente"] if cfg["ambiente"] in ENVIRONMENTS else "citoplasma"
    k_ambiente = ENVIRONMENTS[ambiente_key]
    modo_rapido = bool(cfg["modo_rapido"])

    limite_energetico = cfg.get("limite_energetico")
    try:
        limite_energetico = float(limite_energetico) if limite_energetico is not None else None
        if limite_energetico is not None and not np.isfinite(limite_energetico):
            limite_energetico = None
    except (TypeError, ValueError):
        limite_energetico = None
    tipo_resposta = str(cfg.get("tipo_resposta") or "nenhuma")
    if tipo_resposta not in TIPOS_RESPOSTA:
        tipo_resposta = "nenhuma"

    elementos_permitidos = [e for e in cfg["elementos"] if e in COMPOSTOS_REAIS]
    if not elementos_permitidos:
        elementos_permitidos = list(COMPOSTOS_REAIS.keys())

    run_id = int(time.time() * 1000) % 10_000_000
    dir_proteinas = os.path.join(output_dir, "proteinas")
    dir_virus = os.path.join(output_dir, "virus")
    dir_docking = os.path.join(output_dir, "docking")
    dir_relatorios = os.path.join(output_dir, "relatorios")
    for d in (output_dir, dir_proteinas, dir_virus, dir_docking, dir_relatorios):
        os.makedirs(d, exist_ok=True)

    prever_mutacao = bool(cfg.get("prever_mutacao", False))
    n_mutacoes = max(1, min(int(cfg.get("n_mutacoes", 3) or 3), 10))

    resultado = {"run_id": run_id, "ambiente": ambiente_key,
                 "ambiente_label": k_ambiente["label"], "arquivos": {},
                 "limite_energetico": limite_energetico,
                 "tipo_resposta": tipo_resposta,
                 "proteinas": [], "virus": [], "bacterias": [], "docking": [], "ranking": [],
                 "mutacoes": [], "mutacao_top": None, "resposta": None}

    # ---- 1) compostos reais (biblioteca de LIGANTES p/ triagem) -------
    compostos_carregados = {}
    if cfg["incluir_compostos_reais"]:
        if not RDKIT_DISPONIVEL:
            raise RuntimeError("rdkit não está instalado — necessário para incluir_compostos_reais=True.")
        log("Carregando geometria 3D dos compostos selecionados (ligantes)...")
        for nome_c in elementos_permitidos:
            pts_c, els_c, ligs_c, cargas_c, raios_c, hidro_c, desc_c = carregar_composto_real(nome_c)
            compostos_carregados[nome_c] = {
                "pontos": pts_c, "elementos": els_c, "ligacoes": ligs_c,
                "cargas": cargas_c, "raios": raios_c, "hidro": hidro_c, "descritores": desc_c,
            }

    # ---- 2) proteína(s) nova(s) — gerada por combinação de aminoácidos
    n_iter_dobra = 150 if modo_rapido else 400
    modo_desc = (f"composição fixa ({sum(comp_fixa.values())} AAs)"
                 if comp_fixa else
                 f"{comprimento_prot} resíduos" + (f" | pool: {', '.join(aa_pool)}" if aa_pool else ""))
    log(f"Gerando {n_proteinas} proteína(s) ({modo_desc}) "
        f"no ambiente '{k_ambiente['label']}'...")

    proteinas_info = {}
    tentativas_energia = 1 if limite_energetico is None else (3 if modo_rapido else 6)
    for p_idx in range(n_proteinas):
        melhor_prot = None
        melhor_energia = float("inf")
        tentativas_feitas = 0
        for tentativa in range(tentativas_energia):
            if cfg["seed"] is None:
                seed_p = None
            else:
                seed_p = int(cfg["seed"]) + p_idx * tentativas_energia + tentativa
            candidato = gerar_proteina_nova(
                comprimento=comprimento_prot,
                aminoacidos_pool=aa_pool or None,
                composicao_fixa=comp_fixa or None,
                k_ambiente=k_ambiente,
                n_iter_dobra=n_iter_dobra,
                seed=seed_p,
            )
            energia_candidato = float(candidato.get("energia_final", 1e12))
            if not np.isfinite(energia_candidato):
                # Evita enviar NaN/Infinity ao navegador, que não são JSON válido.
                energia_candidato = 1e12
                candidato["energia_final"] = energia_candidato
            tentativas_feitas += 1
            if melhor_prot is None or energia_candidato < melhor_energia:
                melhor_prot = candidato
                melhor_energia = energia_candidato
            if limite_energetico is not None and energia_candidato <= limite_energetico:
                break

        prot = melhor_prot
        prot["atingiu_limite_energetico"] = (
            limite_energetico is None or float(prot["energia_final"]) <= limite_energetico
        )
        prot["tentativas_energeticas"] = tentativas_feitas
        stats = calcular_estatisticas_proteina(prot["sequencia"])
        nome_prot = f"Proteina_{p_idx + 1}({stats['n_residuos']}aa)"
        proteinas_info[nome_prot] = prot

        scene_prot = trimesh.Scene()
        for geo in criar_geometria_proteina(prot["pontos"], prot["sequencia"]):
            scene_prot.add_geometry(geo)
        scene_prot.apply_translation(-scene_prot.centroid)
        glb = exportar_glb(scene_prot, os.path.join(dir_proteinas, f"{slug(nome_prot)}_{run_id}.glb"))

        resultado["proteinas"].append({
            "nome":           nome_prot,
            "sequencia_fasta": sequencia_para_fasta(prot["sequencia"]),
            "estatisticas":   stats,
            "energia_formacao": round(float(prot["energia_final"]), 4),
            "limite_energetico": limite_energetico,
            "atingiu_limite_energetico": bool(prot["atingiu_limite_energetico"]),
            "tentativas_energeticas": tentativas_feitas,
            "glb":            os.path.relpath(glb, output_dir) if glb else None,
        })
        status_energia = (
            "atingiu o limite" if prot["atingiu_limite_energetico"] else "melhor candidato após as tentativas"
        )
        limite_log = f" | limite={limite_energetico:.4f}" if limite_energetico is not None else ""
        log(f"  {nome_prot}: {stats['n_residuos']} resíduos, "
            f"energia={prot['energia_final']:.4f}{limite_log} ({status_energia}), "
            f"massa~{stats['massa_total_da']:.1f} Da, "
            f"carga liq={stats['carga_liquida']:+.2f}, "
            f"{stats['n_hidrofobo']}🔴 hidrofóbicos / {stats['n_hidrofilico']}🟢 hidrofílicos")

    # a triagem usa as proteínas geradas (+ compostos reais, se ativos)
    alvos_docking_ligantes = dict(proteinas_info)
    if cfg["incluir_compostos_reais"]:
        alvos_docking_ligantes.update(compostos_carregados)

    # ---- 3) vírus fictícios --------------------------------------------
    virus_runtime = []
    if cfg["gerar_virus"]:
        n_res_scale = 0.6 if modo_rapido else 1.0
        log("Gerando os 5 modelos virais fictícios...")
        for vb in VIRUS_MODELOS_BASE:
            v = dict(vb)
            v["raio_nucleo"] = round(v["raio_envelope"] * 0.62, 2)
            n_res = max(6, int(v["n_res_spike"] * n_res_scale))
            v["sequencia"] = gerar_sequencia(n_res, seed=v["seed"])
            v["carga"] = torch.tensor([AMINOACIDOS[a]["carga"] for a in v["sequencia"]])
            v["hidro"] = torch.tensor([AMINOACIDOS[a]["hidrofobicidade"] for a in v["sequencia"]])
            v["raio"] = torch.tensor([AMINOACIDOS[a]["raio_vdw"] for a in v["sequencia"]])
            v["cor_env"] = PALETA_ENVELOPES[v["cor_env_idx"]]
            cores_internas = PALETA_INTERNA[v["cor_env_idx"]]
            v["estilo"] = {
                "raio_envelope": v["raio_envelope"], "raio_nucleo": v["raio_nucleo"],
                "cor_envelope": v["cor_env"], "num_bases_rna": v["n_bases_rna"],
                "n_espiculas": v["n_espiculas"] if not modo_rapido else max(0, v["n_espiculas"] // 2),
                "raio_espicula": v["raio_espicula"], "altura_espicula": v["altura_espicula"],
                "raio_cabeca_espicula": v["raio_cabeca_espicula"],
                "n_glicanos_por_espicula": v["n_glicanos_por_espicula"],
                "amplitude_organica": v["amplitude_organica"],
                "forma_alongada": v["forma_alongada"], "fator_alongamento": v["fator_alongamento"],
                "cor_spike": cores_internas["spike"], "cor_cabeca_espicula_visual": cores_internas["cabeca"],
                "cor_glicano": cores_internas["glicano"], "cor_nucleo_interno": cores_internas["nucleo"],
            }
            geos, pontos_spike = montar_virus_3d(v["sequencia"], v["seed"], v["estilo"])
            scene_virus = trimesh.Scene()
            for g in geos:
                scene_virus.add_geometry(g)
            scene_virus.apply_translation(-scene_virus.centroid)
            nome_arq = slug(v["nome"])
            glb = exportar_glb(scene_virus, os.path.join(dir_virus, f"virus_{nome_arq}_{run_id}.glb"))
            v["pontos_alvo"] = pontos_spike
            v["glb"] = glb
            v["tipo_alvo"] = "virus"
            virus_runtime.append(v)
            resultado["virus"].append({
                "nome": v["nome"],
                "glb": os.path.relpath(glb, output_dir) if glb else None,
            })
            log(f"  {v['nome']} pronto.")

        # ---- 3B) vírus CUSTOMIZADO (genoma digitado pelo usuário, não
        #          um dos 5 pré-definidos) -------------------------------
        vc = cfg.get("virus_customizado")
        if vc and (vc.get("sequencia_acido") or vc.get("preset")):
            # se um preset real foi escolhido (covid19/hiv/influenza) e
            # nenhuma sequência própria foi digitada, usa o trecho real
            # daquele vírus em vez de exigir que o usuário digite algo
            preset_key = vc.get("preset")
            if preset_key and preset_key in VIRUS_REAIS_PRESETS and not vc.get("sequencia_acido"):
                preset = VIRUS_REAIS_PRESETS[preset_key]
                vc = {**vc, "nome": vc.get("nome") or preset["nome"],
                      "tipo_acido": preset["tipo_acido"],
                      "sequencia_acido": preset["sequencia_acido"]}
            tipo_acido = vc.get("tipo_acido", "RNA")
            bases_customizadas = parse_sequencia_acido_nucleico(vc["sequencia_acido"], tipo_acido)
            nome_vc = vc.get("nome") or "Virus-Customizado"
            seed_vc = int(vc.get("seed", 909))
            n_res_spike_vc = int(vc.get("n_res_spike", 15) * (0.6 if modo_rapido else 1.0))
            raio_envelope_vc = float(vc.get("raio_envelope", 6.0))
            v = {
                "nome": nome_vc, "seed": seed_vc,
                "raio_envelope": raio_envelope_vc, "raio_nucleo": round(raio_envelope_vc * 0.62, 2),
            }
            v["sequencia"] = gerar_sequencia(max(6, n_res_spike_vc), seed=seed_vc)
            v["carga"] = torch.tensor([AMINOACIDOS[a]["carga"] for a in v["sequencia"]])
            v["hidro"] = torch.tensor([AMINOACIDOS[a]["hidrofobicidade"] for a in v["sequencia"]])
            v["raio"] = torch.tensor([AMINOACIDOS[a]["raio_vdw"] for a in v["sequencia"]])
            cores_internas = PALETA_INTERNA[seed_vc % len(PALETA_INTERNA)]
            v["estilo"] = {
                "raio_envelope": v["raio_envelope"], "raio_nucleo": v["raio_nucleo"],
                "cor_envelope": PALETA_ENVELOPES[seed_vc % len(PALETA_ENVELOPES)],
                "num_bases_rna": len(bases_customizadas),
                "n_espiculas": int(vc.get("n_espiculas", 30)) if not modo_rapido else int(vc.get("n_espiculas", 30)) // 2,
                "raio_espicula": float(vc.get("raio_espicula", 0.32)),
                "altura_espicula": float(vc.get("altura_espicula", 1.0)),
                "raio_cabeca_espicula": float(vc.get("raio_cabeca_espicula", 0.36)),
                "n_glicanos_por_espicula": int(vc.get("n_glicanos_por_espicula", 3)),
                "amplitude_organica": float(vc.get("amplitude_organica", 0.07)),
                "forma_alongada": bool(vc.get("forma_alongada", False)),
                "fator_alongamento": float(vc.get("fator_alongamento", 1.0)),
                "cor_spike": cores_internas["spike"], "cor_cabeca_espicula_visual": cores_internas["cabeca"],
                "cor_glicano": cores_internas["glicano"], "cor_nucleo_interno": cores_internas["nucleo"],
            }
            genoma_customizado = {"tipo": tipo_acido, "bases": bases_customizadas}
            v["_genoma_customizado"] = genoma_customizado

            # ---- imagem real de referência (opcional) ----------------
            imagem_ref_info = None
            if vc.get("usar_imagem_referencia"):
                termo_img = vc.get("nome_organismo_imagem") or nome_vc
                ref = obter_referencia_visual(termo_img, tipo="virus", log=log)
                if ref:
                    cor_media = np.array(ref["cor_media"], dtype=float)
                    v["estilo"] = aplicar_referencia_no_virus(v["estilo"], cor_media, ref["proporcao"])
                    imagem_ref_info = {
                        "url": ref["fonte_url"], "titulo": ref["titulo"], "thumb_url": ref["thumb_url"],
                        "tipo": ref.get("tipo", "virus"), "termo_consultado": ref.get("termo_consultado", termo_img),
                        "consulta": ref.get("consulta", ""), "autor": ref.get("autor", ""),
                        "licenca": ref.get("licenca", ""),
                    }
                    log(f"  🖼️ Cor/proporção de '{nome_vc}' ajustadas com imagem real "
                        f"({ref['fonte_url']}).")

            geos, pontos_spike = montar_virus_3d(v["sequencia"], v["seed"], v["estilo"],
                                                  genoma_customizado=genoma_customizado)
            scene_virus = trimesh.Scene()
            for g in geos:
                scene_virus.add_geometry(g)
            scene_virus.apply_translation(-scene_virus.centroid)
            nome_arq = slug(v["nome"])
            glb = exportar_glb(scene_virus, os.path.join(dir_virus, f"virus_{nome_arq}_{run_id}.glb"))
            v["pontos_alvo"] = pontos_spike
            v["glb"] = glb
            v["tipo_alvo"] = "virus"
            virus_runtime.append(v)
            resultado["virus"].append({
                "nome": v["nome"], "customizado": True,
                "genoma": f"{tipo_acido}: {''.join(bases_customizadas)}",
                "glb": os.path.relpath(glb, output_dir) if glb else None,
                "imagem_referencia": imagem_ref_info,
            })
            log(f"  {nome_vc} (genoma customizado, {tipo_acido}, {len(bases_customizadas)} bases) pronto.")

    # ---- 3B-bis) previsão de mutações nos vírus ---------------------------
    # Gera, para cada vírus fictício/customizado já criado acima, uma
    # variante MUTANTE (mutações pontuais sorteadas na proteína de
    # espícula) e um índice heurístico de impacto estrutural — um
    # "preditor" simplificado, físico-químico, sem validade biológica
    # real, mas útil para ilustrar como mutações podem mudar a
    # superfície de ligação de um vírus (o mesmo princípio, em espírito,
    # por trás do acompanhamento de variantes de vírus reais).
    if prever_mutacao and virus_runtime:
        log(f"Prevendo variantes mutantes ({n_mutacoes} mutação(ões) por vírus)...")
        for idx_v, v in enumerate(virus_runtime):
            seed_mut = (v["seed"] * 31 + idx_v * 7 + (cfg.get("seed") or 0)) % (2 ** 31 - 1)
            pred = prever_mutacao_virus(v, n_mutacoes, seed_mut, dir_virus, run_id, modo_rapido=modo_rapido)
            v["_predicao_mutacao"] = pred
            resultado["mutacoes"].append({
                "virus": pred["virus"],
                "mutacoes": [{k: m[k] for k in
                              ("posicao", "de", "para", "de_1letra", "para_1letra", "indice_impacto", "classe")}
                             for m in pred["mutacoes"]],
                "indice_impacto_medio": pred["indice_impacto_medio"],
                "classificacao": pred["classificacao"],
                "glb": os.path.relpath(pred["glb"], output_dir) if pred["glb"] else None,
            })
            resumo_mut = ", ".join(f"{m['de_1letra']}{m['posicao']}{m['para_1letra']}" for m in pred["mutacoes"])
            log(f"  {v['nome']}: mutações {resumo_mut} — impacto médio "
                f"{pred['indice_impacto_medio']}/100 ({pred['classificacao']})")

    # ---- 3C) bactérias fictícias -----------------------------------------
    bacteria_runtime = []
    if cfg["gerar_bacterias"]:
        n_bact = max(1, min(int(cfg["n_bacterias"]), 6))
        seed_base_bact = int(cfg.get("seed") or 600)
        log(f"Gerando {n_bact} bactéria(s) fictícia(s) — cada uma com forma, "
            f"tamanho, flagelos e cores diferentes...")

        # genoma real opcional (via NCBI, ver buscar_acido_nucleico_ncbi) —
        # se fornecido, substitui o nucleoide gerado aleatoriamente da
        # PRIMEIRA bactéria pela sequência de DNA/RNA real digitada/buscada
        bact_genoma_cfg = cfg.get("bacteria_genoma_customizado")
        genoma_bact_parseado = None
        if bact_genoma_cfg and bact_genoma_cfg.get("sequencia_acido"):
            tipo_g = bact_genoma_cfg.get("tipo_acido", "DNA")
            bases_g = parse_sequencia_acido_nucleico(bact_genoma_cfg["sequencia_acido"], tipo_g)
            if bases_g:
                genoma_bact_parseado = {"tipo": tipo_g, "bases": bases_g}

        # imagem real de referência (opcional) — mesma ideia do vírus
        # acima: só ajusta cor/proporção da PRIMEIRA bactéria, e só se
        # `usar_imagem_referencia` foi marcado E deu pra achar um termo
        # (nome do organismo escolhido na busca/detalhe do NCBI).
        imagem_ref_bact = None
        if bact_genoma_cfg and bact_genoma_cfg.get("usar_imagem_referencia"):
            termo_img_bact = bact_genoma_cfg.get("nome_organismo_imagem")
            if termo_img_bact:
                ref_b = obter_referencia_visual(termo_img_bact, tipo="bacteria", log=log)
                if ref_b:
                    imagem_ref_bact = {
                        "cor_media": np.array(ref_b["cor_media"], dtype=float),
                        "proporcao": ref_b["proporcao"],
                        "forma_sugerida": ref_b.get("forma_sugerida"),
                        "info": {
                            "url": ref_b["fonte_url"], "titulo": ref_b["titulo"],
                            "thumb_url": ref_b["thumb_url"], "tipo": ref_b.get("tipo", "bacteria"),
                            "termo_consultado": ref_b.get("termo_consultado", termo_img_bact),
                            "consulta": ref_b.get("consulta", ""), "autor": ref_b.get("autor", ""),
                            "licenca": ref_b.get("licenca", ""),
                        },
                    }
                    log(f"  🖼️ Cor/proporção da 1ª bactéria serão ajustadas com imagem real "
                        f"de '{termo_img_bact}' ({ref_b['fonte_url']}).")
            else:
                log("  ⚠ 'usar imagem de referência' marcado, mas nenhum organismo foi "
                    "identificado na busca — pulando o ajuste por imagem.")

        for i_b in range(n_bact):
            b = gerar_estilo_bacteria(i_b, seed_base=seed_base_bact)
            paleta = PALETA_BACTERIAS[b["paleta_idx"] % len(PALETA_BACTERIAS)]
            b["cor_parede"] = paleta["parede"]
            b["cor_citoplasma"] = paleta["citoplasma"]
            b["cor_proteina_superficie"] = paleta["proteina_superficie"]
            n_res_prot = max(6, int(b["n_res_proteina_superficie"] * (0.6 if modo_rapido else 1.0)))
            b["sequencia_proteina_superficie"] = gerar_sequencia(n_res_prot, seed=b["seed"])
            b["carga"] = torch.tensor([AMINOACIDOS[a]["carga"] for a in b["sequencia_proteina_superficie"]])
            b["hidro"] = torch.tensor([AMINOACIDOS[a]["hidrofobicidade"] for a in b["sequencia_proteina_superficie"]])
            b["raio"] = torch.tensor([AMINOACIDOS[a]["raio_vdw"] for a in b["sequencia_proteina_superficie"]])
            genoma_desta = genoma_bact_parseado if i_b == 0 else None
            b["_genoma_customizado"] = genoma_desta
            imagem_ref_info_b = None
            if i_b == 0 and imagem_ref_bact:
                b = aplicar_referencia_na_bacteria(
                    b, imagem_ref_bact["cor_media"], imagem_ref_bact["proporcao"],
                    forma_sugerida=imagem_ref_bact.get("forma_sugerida")
                )
                imagem_ref_info_b = imagem_ref_bact["info"]
            geos, pontos_prot = montar_bacteria_3d(
                b["sequencia_proteina_superficie"], b["seed"], b, genoma_customizado=genoma_desta)
            scene_bact = trimesh.Scene()
            for g in geos:
                scene_bact.add_geometry(g)
            scene_bact.apply_translation(-scene_bact.centroid)
            nome_arq = slug(b["nome"])
            glb = exportar_glb(scene_bact, os.path.join(dir_virus, f"bacteria_{nome_arq}_{run_id}.glb"))
            b["pontos_alvo"] = pontos_prot
            b["glb"] = glb
            b["tipo_alvo"] = "bacteria"
            bacteria_runtime.append(b)
            resultado["bacterias"].append({
                "nome": b["nome"],
                "glb": os.path.relpath(glb, output_dir) if glb else None,
                "estatisticas": calcular_estatisticas_bacteria(b),
                "imagem_referencia": imagem_ref_info_b,
            })
            log(f"  {b['nome']} pronta ({b['forma']}, {b['n_flagelos']} flagelo(s), "
                f"{b['n_plasmideos']} plasmídeo(s)).")

    alvos_runtime = virus_runtime + bacteria_runtime

    # ---- 4) triagem / docking (contra vírus E bactérias) ------------------
    resultados_triagem = []
    if alvos_runtime and alvos_docking_ligantes:
        iter_dock = 80 if modo_rapido else 250
        n_starts = 1 if modo_rapido else 3
        log("Rodando docking (triagem virtual) contra os vírus/bactérias fictícios...")
        for v in alvos_runtime:
            p_vir = torch.tensor(v["pontos_alvo"], dtype=torch.float32)
            for nome_c, comp in alvos_docking_ligantes.items():
                seed_par = hash((nome_c, v["nome"])) % (2 ** 31 - 1)
                p_final, energia = dockar_molecula_corpo_rigido(
                    comp["pontos"], comp["cargas"], comp["raios"], comp["hidro"],
                    p_vir, v["carga"], v["raio"], v["hidro"],
                    n_iter=iter_dock, seed=seed_par, n_starts=n_starts,
                )
                tipo = "proteina" if "sequencia" in comp else "composto_real"
                resultados_triagem.append({
                    "composto": nome_c, "virus": v["nome"], "energia": energia,
                    "p_virus": v["pontos_alvo"], "p_composto": p_final,
                    "tipo": tipo, "tipo_alvo": v["tipo_alvo"],
                    "sequencia": comp.get("sequencia"),
                    "elementos": comp.get("elementos"), "ligacoes": comp["ligacoes"],
                    "alvo_runtime": v,
                })

        resultados_triagem.sort(key=lambda r: r["energia"])
        n_res = len(resultados_triagem)
        for i, r in enumerate(resultados_triagem):
            r["score"] = 100.0 * (n_res - 1 - i) / (n_res - 1) if n_res > 1 else 50.0
        resultados_triagem.sort(key=lambda r: r["score"], reverse=True)

        for r in resultados_triagem:
            resultado["ranking"].append({
                "composto": r["composto"], "virus": r["virus"],
                "energia": round(float(r["energia"]), 4), "score": round(float(r["score"]), 1),
            })

        if resultados_triagem:
            melhor = resultados_triagem[0]
            scene_dock = trimesh.Scene()
            alvo_r = melhor["alvo_runtime"]
            if melhor["tipo_alvo"] == "bacteria":
                geos_alvo, _ = montar_bacteria_3d(alvo_r["sequencia_proteina_superficie"], alvo_r["seed"], alvo_r,
                                                   genoma_customizado=alvo_r.get("_genoma_customizado"))
            else:
                geos_alvo, _ = montar_virus_3d(alvo_r["sequencia"], alvo_r["seed"], alvo_r["estilo"])
            for g in geos_alvo:
                scene_dock.add_geometry(g)
            if melhor["tipo"] == "proteina":
                geos_comp = criar_geometria_proteina(melhor["p_composto"], melhor["sequencia"])
            else:
                geos_comp = criar_geometria_composto_real(melhor["p_composto"], melhor["elementos"], melhor["ligacoes"])
            for g in geos_comp:
                scene_dock.add_geometry(g)
            scene_dock.apply_translation(-scene_dock.centroid)
            glb_dock = exportar_glb(scene_dock, os.path.join(
                dir_docking, f"docking_{slug(melhor['composto'])}_x_{slug(melhor['virus'])}_{run_id}.glb"))
            resultado["docking"].append({
                "composto": melhor["composto"], "virus": melhor["virus"],
                "score": round(float(melhor["score"]), 1),
                "glb": os.path.relpath(glb_dock, output_dir) if glb_dock else None,
            })
            log(f"Melhor encaixe: {melhor['composto']} x {melhor['virus']} (score {melhor['score']:.1f}/100)")

            # ---- 4B) previsão de mutação aplicada ao par vencedor --------
            # Redocka o MESMO ligante contra a variante mutante do vírus
            # vencedor, pra comparar a energia de ligação antes/depois da
            # mutação (em vez de só o índice heurístico rápido calculado
            # acima) — essa é a única mutação que passa pelo docking
            # completo de novo, pra manter o custo computacional baixo.
            if prever_mutacao and melhor["tipo_alvo"] == "virus":
                pred = alvo_r.get("_predicao_mutacao")
                comp = alvos_docking_ligantes.get(melhor["composto"])
                if pred and comp is not None:
                    seq_mut = pred["_seq_mutada"]
                    p_spike_mut = torch.tensor(pred["_pontos_spike_mutados"], dtype=torch.float32)
                    carga_mut = torch.tensor([AMINOACIDOS[a]["carga"] for a in seq_mut])
                    hidro_mut = torch.tensor([AMINOACIDOS[a]["hidrofobicidade"] for a in seq_mut])
                    raio_mut = torch.tensor([AMINOACIDOS[a]["raio_vdw"] for a in seq_mut])
                    seed_par_m = hash((melhor["composto"], alvo_r["nome"], "mutante")) % (2 ** 31 - 1)
                    p_final_mut, energia_mut = dockar_molecula_corpo_rigido(
                        comp["pontos"], comp["cargas"], comp["raios"], comp["hidro"],
                        p_spike_mut, carga_mut, raio_mut, hidro_mut,
                        n_iter=iter_dock, seed=seed_par_m, n_starts=n_starts,
                    )
                    delta_energia = float(energia_mut) - float(melhor["energia"])
                    if delta_energia > 0.5:
                        interpretacao = "ligação prevista MAIS FRACA no mutante (possível fuga/resistência)"
                    elif delta_energia < -0.5:
                        interpretacao = "ligação prevista MAIS FORTE no mutante"
                    else:
                        interpretacao = "impacto previsto pequeno na ligação"

                    scene_mut = trimesh.Scene()
                    geos_mut_alvo, _ = montar_virus_3d(
                        seq_mut, alvo_r["seed"], alvo_r["estilo"],
                        genoma_customizado=alvo_r.get("_genoma_customizado"),
                        indices_mutados=[m["posicao"] - 1 for m in pred["mutacoes"]])
                    for g in geos_mut_alvo:
                        scene_mut.add_geometry(g)
                    if melhor["tipo"] == "proteina":
                        geos_comp_mut = criar_geometria_proteina(p_final_mut, melhor["sequencia"])
                    else:
                        geos_comp_mut = criar_geometria_composto_real(
                            p_final_mut, melhor["elementos"], melhor["ligacoes"])
                    for g in geos_comp_mut:
                        scene_mut.add_geometry(g)
                    scene_mut.apply_translation(-scene_mut.centroid)
                    glb_mut = exportar_glb(scene_mut, os.path.join(
                        dir_docking,
                        f"docking_mutante_{slug(melhor['composto'])}_x_{slug(melhor['virus'])}_{run_id}.glb"))

                    resultado["mutacao_top"] = {
                        "virus": melhor["virus"], "composto": melhor["composto"],
                        "mutacoes": [{k: m[k] for k in
                                      ("posicao", "de", "para", "de_1letra", "para_1letra",
                                       "indice_impacto", "classe")}
                                     for m in pred["mutacoes"]],
                        "energia_original": round(float(melhor["energia"]), 4),
                        "energia_mutante": round(float(energia_mut), 4),
                        "delta_energia": round(delta_energia, 4),
                        "interpretacao": interpretacao,
                        "glb": os.path.relpath(glb_mut, output_dir) if glb_mut else None,
                    }
                    log(f"Previsão de mutação no vencedor ({melhor['virus']}): {interpretacao} "
                        f"(Δenergia={delta_energia:+.3f})")

    # ---- 4C) interpretação da estratégia de resposta ----------------------
    resultado["resposta"] = construir_resposta_simulada(
        tipo_resposta,
        alvos_runtime,
        resultados_triagem,
        limite_energetico=limite_energetico,
    )
    log(f"Estratégia de resposta: {resultado['resposta']['label']}. "
        f"{resultado['resposta']['resumo']}")

    # ---- 5) relatório + gráfico ------------------------------------------
    rel_path = os.path.join(dir_relatorios, f"relatorio_{run_id}.txt")
    with open(rel_path, "w", encoding="utf-8") as f:
        f.write("RELATÓRIO DE TRIAGEM VIRTUAL (modelo de brinquedo, sem validade biológica)\n")
        f.write(f"Ambiente de dobramento: {k_ambiente['label']}\n")
        f.write(f"Limite energético: {limite_energetico if limite_energetico is not None else 'sem limite'}\n")
        f.write(f"Estratégia de resposta: {resultado['resposta']['label']}\n")
        f.write(f"Proteínas geradas: {n_proteinas} | comprimento: {comprimento_prot} resíduos\n")
        if comp_fixa:
            f.write(f"Composição fixa: {comp_fixa}\n")
        elif aa_pool:
            f.write(f"Pool de aminoácidos: {', '.join(aa_pool)}\n")
        if cfg["incluir_compostos_reais"]:
            f.write(f"Ligantes reais incluídos: {', '.join(elementos_permitidos)}\n")
        f.write("\n")
        for p in resultado["proteinas"]:
            s = p["estatisticas"]
            f.write(f"- {p['nome']}  |  {s['n_residuos']} resíduos  "
                    f"|  energia={p.get('energia_formacao')}  "
                    f"|  massa~{s['massa_total_da']:.1f} Da  |  carga liq={s['carga_liquida']:+.2f}  "
                    f"|  {s['n_hidrofobo']} hidrofóbicos / {s['n_hidrofilico']} hidrofílicos\n")
        f.write("\nRanking (melhor primeiro):\n")
        for r in resultado["ranking"]:
            f.write(f"  {r['composto']} x {r['virus']}: score={r['score']}/100 (energia={r['energia']})\n")
        f.write(f"\nResposta: {resultado['resposta']['resumo']}\n")
        for alvo_resp in resultado["resposta"].get("alvos", []):
            f.write(f"  {alvo_resp['nome']}: {alvo_resp['status']} "
                    f"(energia={alvo_resp.get('energia_melhor')}, "
                    f"candidato={alvo_resp.get('ligante_melhor')})\n")
        if resultado["mutacoes"]:
            f.write("\nPrevisão de mutações (modelo heurístico, sem validade biológica real):\n")
            for m in resultado["mutacoes"]:
                resumo = ", ".join(f"{mm['de_1letra']}{mm['posicao']}{mm['para_1letra']}" for mm in m["mutacoes"])
                f.write(f"  {m['virus']}: {resumo} -> impacto médio {m['indice_impacto_medio']}/100 "
                        f"({m['classificacao']})\n")
        if resultado["mutacao_top"]:
            mt = resultado["mutacao_top"]
            f.write(f"\nPar vencedor re-docado com a variante mutante:\n")
            f.write(f"  {mt['composto']} x {mt['virus']} (mutado): energia {mt['energia_original']} -> "
                    f"{mt['energia_mutante']} (Δ={mt['delta_energia']:+.3f}) — {mt['interpretacao']}\n")
    resultado["arquivos"]["relatorio"] = os.path.relpath(rel_path, output_dir)

    if MATPLOTLIB_DISPONIVEL and resultado["ranking"]:
        try:
            top = resultado["ranking"][:15]
            nomes = [f"{r['composto']} x {r['virus']}" for r in top]
            scores = [r["score"] for r in top]
            plt.figure(figsize=(9, max(4, len(top) * 0.4)))
            bars = plt.barh(nomes, scores, color="#3d8bfd")
            plt.gca().invert_yaxis()
            plt.xlabel("Score ilustrativo (0-100)")
            plt.title(f"Ranking — ambiente: {k_ambiente['label']}")
            plt.xlim(0, 105)
            for bar, sc in zip(bars, scores):
                plt.text(bar.get_width() + 1, bar.get_y() + bar.get_height() / 2, f"{sc:.1f}", va="center", fontsize=9)
            plt.tight_layout()
            graf_path = os.path.join(dir_relatorios, f"grafico_{run_id}.png")
            plt.savefig(graf_path, dpi=140)
            plt.close()
            resultado["arquivos"]["grafico"] = os.path.relpath(graf_path, output_dir)
        except Exception as e:
            log(f"[AVISO] Gráfico não gerado: {e}")

    log("Simulação concluída.")
    return resultado