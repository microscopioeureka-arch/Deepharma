# -*- coding: utf-8 -*-
"""
rodar_local.py — roda a simulação direto no VS Code, sem precisar do site.
Edite CONFIG abaixo e rode com:
    python rodar_local.py
"""
import os
import simulacao

CONFIG = {
    # ── Proteína nova por combinação de aminoácidos ──────────────────
    "comprimento_proteina": 30,   # nº de resíduos da cadeia (5–80)
    "n_proteinas": 2,             # quantas proteínas gerar (1–5)

    # Modo A — comprimento livre, pool restrito (deixe [] para usar todos os 20)
    "aminoacidos_pool": ["ALA", "LEU", "GLY", "LYS", "GLU", "PHE", "SER"],

    # Modo B — composição exata (ignorará comprimento_proteina e aminoacidos_pool)
    # Descomente o bloco abaixo e comente o aminoacidos_pool acima:
    # "composicao_fixa": {"ALA": 5, "LEU": 8, "LYS": 3, "GLU": 4, "GLY": 2},

    # ── Ambiente ────────────────────────────────────────────────────
    "ambiente": "citoplasma",   # citoplasma | membrana | acido | alcalino | salino
    "modo_rapido": False,

    # ── Alvos de docking ────────────────────────────────────────────
    "gerar_virus": True,
    "gerar_bacterias": True,
    "n_bacterias": 3,   # 1 a 6 — cada bactéria gerada é diferente das outras

    # Genoma real opcional para a 1ª bactéria (busque com
    # simulacao.buscar_acido_nucleico_ncbi("nome do organismo") e cole aqui):
    # "bacteria_genoma_customizado": {
    #     "tipo_acido": "DNA",
    #     "sequencia_acido": "ATGCATGCATGC...",
    # },

    # Vírus real pré-definido (preset):
    "virus_customizado": {"preset": "covid19"},  # "covid19" | "hiv" | "influenza"

    # OU vírus com genoma digitado (descomente e comente a linha acima):
    # "virus_customizado": {
    #     "nome": "Meu-Virus",
    #     "tipo_acido": "RNA",
    #     "sequencia_acido": "AUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGC",
    # },

    # ── Ligantes reais (fármacos) ────────────────────────────────────
    "incluir_compostos_reais": True,
    "elementos": ["Paracetamol", "Cafeina", "Cloroquina", "Favipiravir"],

    # ── Previsão de mutações nos vírus ───────────────────────────────
    # Gera, para cada vírus fictício/customizado, uma variante mutante
    # (mutações pontuais sorteadas na proteína de espícula) + um índice
    # heurístico de impacto (0-100). Para o par vencedor do docking, além
    # disso, refaz o docking com a variante mutante pra comparar a
    # energia de ligação antes/depois. Modelo ilustrativo, sem validade
    # biológica real.
    "prever_mutacao": True,
    "n_mutacoes": 3,   # nº de mutações pontuais por vírus (1-10)
}

if __name__ == "__main__":
    saida = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs", "execucao_local")
    print(f"Ambientes: {list(simulacao.ENVIRONMENTS.keys())}")
    print(f"AAs disponíveis: {list(simulacao.AMINOACIDOS.keys())}")
    print(f"Vírus reais: {list(simulacao.VIRUS_REAIS_PRESETS.keys())}")
    resultado = simulacao.executar_simulacao(CONFIG, saida, progresso=print)
    print("\n─── RESUMO ───")
    print(f"Pasta de saída: {os.path.abspath(saida)}")
    print(f"Proteínas geradas: {[p['nome'] for p in resultado['proteinas']]}")
    print(f"Bactérias: {[b['nome'] for b in resultado['bacterias']]}")
    print(f"Vírus: {[v['nome'] for v in resultado['virus']]}")
    print(f"Melhor docking: {resultado['docking']}")
    if resultado.get("mutacoes"):
        print(f"Previsões de mutação: {[(m['virus'], m['classificacao']) for m in resultado['mutacoes']]}")
    if resultado.get("mutacao_top"):
        print(f"Mutação do vencedor: {resultado['mutacao_top']['interpretacao']}")
