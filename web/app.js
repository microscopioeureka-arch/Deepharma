// app.js — Simulador Molecular
// Troque API por "https://SEU-BACKEND.onrender.com" se hospedar separado.
const API = "";

// Todas as chamadas são para a mesma origem. Se a sessão expirar,
// manda o usuário de volta para a tela de login.
async function fetchAutenticado(url, options) {
  const resposta = await fetch(url, options);
  if (resposta.status === 401) {
    window.location.href = "/login";
    throw new Error("Sessão expirada.");
  }
  return resposta;
}

// Lê JSON somente depois de verificar o status. Se o proxy devolver uma
// página HTML (por exemplo, o 502 padrão do Render), evita o erro genérico
// "Unexpected token '<'" e mostra a causa real.
async function lerJson(resposta, contexto = "API") {
  const texto = await resposta.text();
  if (!resposta.ok) {
    const resumo = texto.replace(/<[^>]*>/g, " ").replaceAll("\n", " ").replaceAll("\r", " ").trim().slice(0, 240);
    throw new Error(`${contexto}: HTTP ${resposta.status}${resumo ? ` — ${resumo}` : ""}`);
  }
  try {
    return JSON.parse(texto);
  } catch (_) {
    const tipo = resposta.headers.get("content-type") || "desconhecido";
    throw new Error(`${contexto}: o servidor devolveu ${tipo}, não JSON. Verifique a URL da API e os logs do backend.`);
  }
}


// ── Referências DOM ──────────────────────────────────────────────────
const elComprimento       = document.getElementById("comprimento-proteina");
const elNProteinas        = document.getElementById("n-proteinas");
const elUsarComposicao    = document.getElementById("usar-composicao-fixa");
const elBlocoComposicao   = document.getElementById("bloco-composicao");
const elListaComposicao   = document.getElementById("lista-composicao");
const elUsarPool          = document.getElementById("usar-pool");
const elBlocoPool         = document.getElementById("bloco-pool");
const elListaPool         = document.getElementById("lista-pool");
const elSelectAmbiente    = document.getElementById("select-ambiente");
const elUsarLimiteEnergetico = document.getElementById("usar-limite-energetico");
const elLimiteEnergetico  = document.getElementById("limite-energetico");
const elTipoResposta      = document.getElementById("tipo-resposta");
const elDescricaoResposta = document.getElementById("descricao-resposta");
const elGerarVirus        = document.getElementById("gerar-virus");
const elSelectPreset      = document.getElementById("select-virus-preset");
const elFontePreset       = document.getElementById("fonte-preset");
const elAtivarVc          = document.getElementById("ativar-virus-customizado");
const elBlocoVc           = document.getElementById("bloco-virus-customizado");
const elVcNome            = document.getElementById("vc-nome");
const elVcTipo            = document.getElementById("vc-tipo");
const elVcSeq             = document.getElementById("vc-sequencia");
const elGerarBacterias    = document.getElementById("gerar-bacterias");
const elCampoNBact        = document.getElementById("campo-n-bacterias");
const elAtivarBc          = document.getElementById("ativar-bact-customizada");
const elBlocoBc           = document.getElementById("bloco-bact-customizada");
const elBcTipo            = document.getElementById("bc-tipo");
const elBcSeq             = document.getElementById("bc-sequencia");
const elBuscaVirusTermo   = document.getElementById("busca-virus-termo");
const elBtnBuscarVirus    = document.getElementById("btn-buscar-virus");
const elBuscaVirusStatus  = document.getElementById("busca-virus-status");
const elResultadosVirus   = document.getElementById("resultados-virus");
const elBuscaBactTermo    = document.getElementById("busca-bacteria-termo");
const elBtnBuscarBact     = document.getElementById("btn-buscar-bacteria");
const elBuscaBactStatus   = document.getElementById("busca-bacteria-status");
const elResultadosBact    = document.getElementById("resultados-bacteria");
const elBuscaProtTermo    = document.getElementById("busca-proteina-termo");
const elBtnBuscarProt     = document.getElementById("btn-buscar-proteina");
const elBuscaProtStatus   = document.getElementById("busca-proteina-status");
const elResultadosProt    = document.getElementById("resultados-proteina");
const elPreverMutacao     = document.getElementById("prever-mutacao");
const elCampoNMutacoes    = document.getElementById("campo-n-mutacoes");
const elNMutacoes         = document.getElementById("n-mutacoes");
const elIncluirReais      = document.getElementById("incluir-reais");
const elBlocoLigantes     = document.getElementById("bloco-ligantes");
const elListaElementos    = document.getElementById("lista-elementos");
const elBtnSimular        = document.getElementById("btn-simular");
const elStatusTexto       = document.getElementById("status-texto");
const elStatusLog         = document.getElementById("status-log");
const elVisor             = document.getElementById("visor-3d");
const elVisorLegenda      = document.getElementById("visor-legenda");
const elSecaoResultados   = document.getElementById("secao-resultados");
const elListaResposta      = document.getElementById("lista-resposta");
const elBtnSalvar         = document.getElementById("btn-salvar-seq");
const elBtnVerSalvas      = document.getElementById("btn-ver-salvas");
const elSalvas            = document.getElementById("sequencias-salvas");
const elUsarImagemRef     = document.getElementById("usar-imagem-referencia");

// guarda o nome do organismo do último resultado ESCOLHIDO na busca
// (vem do campo "organismo" da lista do NCBI) — usado como termo de
// busca da imagem real de referência, se a opção acima for marcada.
let ultimoOrganismoVirus = "";
let ultimoOrganismoBact = "";

function atualizarDescricaoResposta() {
  const opcao = elTipoResposta?.selectedOptions?.[0];
  if (elDescricaoResposta) {
    elDescricaoResposta.textContent = opcao?.dataset?.descricao || "";
  }
}

function atualizarLimiteEnergetico() {
  if (!elUsarLimiteEnergetico || !elLimiteEnergetico) return;
  elLimiteEnergetico.disabled = !elUsarLimiteEnergetico.checked;
}

// ── Carregar opções da API ───────────────────────────────────────────
async function carregarOpcoes() {
  try {
    const resp = await fetchAutenticado(`${API}/api/opcoes`);
    const d = await lerJson(resp, "Carregar opções");

    // Composição exata — um campo numérico por AA
    elListaComposicao.innerHTML = "";
    (d.aminoacidos || []).forEach((aa) => {
      const div = document.createElement("div");
      div.className = "campo-aa";
      div.innerHTML = `<label>${aa}</label>
        <input type="number" class="aa-qtd" data-aa="${aa}"
          min="0" max="30" value="0" />`;
      elListaComposicao.appendChild(div);
    });

    // Pool de AAs (checkboxes)
    elListaPool.innerHTML = "";
    (d.aminoacidos || []).forEach((aa) => {
      const label = document.createElement("label");
      label.className = "linha-checkbox";
      label.innerHTML = `<input type="checkbox" class="pool-check" value="${aa}" checked/> ${aa}`;
      elListaPool.appendChild(label);
    });

    // Ambiente
    elSelectAmbiente.innerHTML = "";
    (d.ambientes || []).forEach((a) => {
      const o = document.createElement("option");
      o.value = a.chave; o.textContent = a.label;
      elSelectAmbiente.appendChild(o);
    });

    // Estratégias de resposta contra vírus e bactérias
    if (elTipoResposta) {
      elTipoResposta.innerHTML = "";
      (d.tipos_resposta || []).forEach((r) => {
        const o = document.createElement("option");
        o.value = r.chave;
        o.textContent = r.nome;
        o.dataset.descricao = r.descricao || "";
        elTipoResposta.appendChild(o);
      });
      atualizarDescricaoResposta();
    }

    // Ligantes reais
    elListaElementos.innerHTML = "";
    (d.compostos || []).forEach((nome, i) => {
      const lbl = document.createElement("label");
      lbl.className = "linha-checkbox";
      lbl.innerHTML = `<input type="checkbox" value="${nome}" ${i < 3 ? "checked" : ""}/> ${nome}`;
      elListaElementos.appendChild(lbl);
    });

    // Presets de vírus reais
    const presetsInfo = {};
    (d.virus_presets || []).forEach((p) => {
      const o = document.createElement("option");
      o.value = p.chave; o.textContent = p.nome;
      elSelectPreset.appendChild(o);
      presetsInfo[p.chave] = p;
    });
    elSelectPreset._info = presetsInfo;

  } catch (e) { console.error("Erro ao carregar opções:", e); }
}

// ── Eventos ──────────────────────────────────────────────────────────
elUsarLimiteEnergetico?.addEventListener("change", atualizarLimiteEnergetico);
elTipoResposta?.addEventListener("change", atualizarDescricaoResposta);
atualizarLimiteEnergetico();

elUsarComposicao.addEventListener("change", () => {
  elBlocoComposicao.style.display = elUsarComposicao.checked ? "block" : "none";
  if (elUsarComposicao.checked) {
    elUsarPool.checked = false;
    elBlocoPool.style.display = "none";
  }
});

elUsarPool.addEventListener("change", () => {
  elBlocoPool.style.display = elUsarPool.checked ? "block" : "none";
  if (elUsarPool.checked) {
    elUsarComposicao.checked = false;
    elBlocoComposicao.style.display = "none";
  }
});

elSelectPreset.addEventListener("change", () => {
  const info = (elSelectPreset._info || {})[elSelectPreset.value];
  elFontePreset.textContent = info ? (info.descricao || info.fonte || "") : "";
  // Para um preset real, a comparação científica fica ligada por padrão;
  // o usuário ainda pode desmarcá-la nas opções.
  if (info) elUsarImagemRef.checked = true;
});

elAtivarVc.addEventListener("change", () => {
  elBlocoVc.style.display = elAtivarVc.checked ? "block" : "none";
});

elGerarBacterias.addEventListener("change", () => {
  elCampoNBact.style.display = elGerarBacterias.checked ? "flex" : "none";
});

elAtivarBc.addEventListener("change", () => {
  elBlocoBc.style.display = elAtivarBc.checked ? "block" : "none";
});

elPreverMutacao.addEventListener("change", () => {
  elCampoNMutacoes.style.display = elPreverMutacao.checked ? "flex" : "none";
});

elIncluirReais.addEventListener("change", () => {
  elBlocoLigantes.style.display = elIncluirReais.checked ? "block" : "none";
});

// ── Salvar / carregar sequências (localStorage) ───────────────────────
const LS_KEY = "sim_sequencias";

elBtnSalvar.addEventListener("click", () => {
  const nome = elVcNome.value.trim() || "Minha Sequência";
  const seq  = elVcSeq.value.trim();
  if (seq.length < 3) { alert("Digite ao menos 3 bases antes de salvar."); return; }
  const salvas = JSON.parse(localStorage.getItem(LS_KEY) || "[]");
  const idx = salvas.findIndex(s => s.nome === nome);
  const entrada = { nome, tipo: elVcTipo.value, seq, data: new Date().toLocaleDateString("pt-BR") };
  if (idx >= 0) {
    if (!confirm(`Substituir a sequência "${nome}" já salva?`)) return;
    salvas[idx] = entrada;
  } else { salvas.push(entrada); }
  localStorage.setItem(LS_KEY, JSON.stringify(salvas));
  renderSalvas();
  alert(`✅ "${nome}" salva!`);
});

elBtnVerSalvas.addEventListener("click", () => {
  const vis = elSalvas.style.display !== "none";
  elSalvas.style.display = vis ? "none" : "block";
  if (!vis) renderSalvas();
});

function renderSalvas() {
  const salvas = JSON.parse(localStorage.getItem(LS_KEY) || "[]");
  elSalvas.innerHTML = "";
  if (!salvas.length) {
    elSalvas.innerHTML = "<p class='ajuda' style='margin:4px 0;'>Nenhuma sequência salva.</p>";
    return;
  }
  salvas.forEach((s, i) => {
    const d = document.createElement("div");
    d.className = "item-seq-salva";
    d.innerHTML = `<span class="seq-info"><b>${s.nome}</b><small>${s.tipo} · ${s.data}</small></span>
      <span class="seq-acoes">
        <button onclick="carregarSalva(${i})">↩ Usar</button>
        <button onclick="excluirSalva(${i})">✕</button>
      </span>`;
    elSalvas.appendChild(d);
  });
}
window.carregarSalva = (i) => {
  const s = JSON.parse(localStorage.getItem(LS_KEY) || "[]")[i];
  if (!s) return;
  elVcNome.value = s.nome; elVcTipo.value = s.tipo; elVcSeq.value = s.seq;
  elAtivarVc.checked = true; elBlocoVc.style.display = "block";
};
window.excluirSalva = (i) => {
  if (!confirm("Excluir esta sequência salva?")) return;
  const sal = JSON.parse(localStorage.getItem(LS_KEY) || "[]");
  sal.splice(i, 1);
  localStorage.setItem(LS_KEY, JSON.stringify(sal));
  renderSalvas();
};

// ── Buscar dados reais na internet (NCBI / UniProt) ────────────────────
// Fluxo em DUAS etapas: 1) lista candidatos (rápido, não baixa a
// sequência inteira) 2) o usuário clica no que quer e SÓ ENTÃO a
// sequência completa daquele item é baixada e aplicada nos campos.

async function buscarComStatus(btn, elStatus, url, aoSucesso) {
  const original = btn.textContent;
  btn.disabled = true;
  btn.textContent = "Buscando...";
  elStatus.textContent = "";
  elStatus.style.color = "";
  try {
    const r = await fetchAutenticado(url);
    const d = await lerJson(r, "Buscar dados");
    aoSucesso(d);
  } catch (e) {
    elStatus.textContent = `⚠ ${e.message || e}`;
    elStatus.style.color = "#e05c5c";
  } finally {
    btn.disabled = false;
    btn.textContent = original;
  }
}

// Igual a buscarComStatus, mas devolve o array de resultados em vez de
// chamar um callback (usado para a etapa de LISTAGEM).
async function buscarListaComStatus(btn, elStatus, url) {
  const original = btn.textContent;
  btn.disabled = true;
  btn.textContent = "Buscando...";
  elStatus.textContent = "";
  elStatus.style.color = "";
  try {
    const r = await fetchAutenticado(url);
    const d = await lerJson(r, "Listar resultados");
    return d.resultados || [];
  } catch (e) {
    elStatus.textContent = `⚠ ${e.message || e}`;
    elStatus.style.color = "#e05c5c";
    return null;
  } finally {
    btn.disabled = false;
    btn.textContent = original;
  }
}

function esconderResultados(container) {
  container.style.display = "none";
  container.innerHTML = "";
}

// Renderiza a lista de candidatos como itens clicáveis. `montarLinha`
// monta o HTML interno de cada item; `aoEscolher` roda quando o item é
// clicado.
function renderResultados(container, itens, montarLinha, aoEscolher) {
  container.innerHTML = "";
  if (!itens || !itens.length) { container.style.display = "none"; return; }
  const cab = document.createElement("div");
  cab.className = "resultados-busca-cabecalho";
  cab.textContent = `${itens.length} resultado(s) — clique no que você quer usar:`;
  container.appendChild(cab);
  itens.forEach((item) => {
    const div = document.createElement("div");
    div.className = "resultado-item";
    div.innerHTML = montarLinha(item);
    div.addEventListener("click", () => aoEscolher(item));
    container.appendChild(div);
  });
  container.style.display = "block";
}

elBtnBuscarVirus.addEventListener("click", async () => {
  const termo = elBuscaVirusTermo.value.trim();
  if (!termo) { elBuscaVirusStatus.textContent = "Digite um termo de busca."; return; }
  esconderResultados(elResultadosVirus);
  const itens = await buscarListaComStatus(elBtnBuscarVirus, elBuscaVirusStatus,
    `${API}/api/buscar/acido_nucleico/lista?termo=${encodeURIComponent(termo)}&max_resultados=8`);
  if (!itens) return;
  elBuscaVirusStatus.style.color = "#9aa0ac";
  elBuscaVirusStatus.textContent = `${itens.length} resultado(s) no NCBI — escolha um abaixo.`;
  renderResultados(elResultadosVirus, itens, (item) => `
      <b>${item.titulo}</b>
      <small>${item.organismo ? item.organismo + " · " : ""}${item.accession}${
        item.tamanho_bases ? " · " + item.tamanho_bases.toLocaleString("pt-BR") + " bases" : ""}</small>
    `, (item) => {
    esconderResultados(elResultadosVirus);
    ultimoOrganismoVirus = item.organismo || item.titulo || "";
    elUsarImagemRef.checked = true;
    buscarComStatus(elBtnBuscarVirus, elBuscaVirusStatus,
      `${API}/api/buscar/acido_nucleico/detalhe?accession=${encodeURIComponent(item.accession)}&tipo=RNA`,
      (d) => {
        elVcNome.value = d.nome;
        elVcTipo.value = d.tipo_acido;
        elVcSeq.value = d.sequencia_acido;
        elAtivarVc.checked = true;
        elBlocoVc.style.display = "block";
        elBuscaVirusStatus.style.color = "#50c878";
        elBuscaVirusStatus.textContent = `✅ ${d.fonte}`;
      });
  });
});

elBtnBuscarBact.addEventListener("click", async () => {
  const termo = elBuscaBactTermo.value.trim();
  if (!termo) { elBuscaBactStatus.textContent = "Digite um termo de busca."; return; }
  esconderResultados(elResultadosBact);
  const itens = await buscarListaComStatus(elBtnBuscarBact, elBuscaBactStatus,
    `${API}/api/buscar/acido_nucleico/lista?termo=${encodeURIComponent(termo)}&max_resultados=8`);
  if (!itens) return;
  elBuscaBactStatus.style.color = "#9aa0ac";
  elBuscaBactStatus.textContent = `${itens.length} resultado(s) no NCBI — escolha um abaixo.`;
  renderResultados(elResultadosBact, itens, (item) => `
      <b>${item.titulo}</b>
      <small>${item.organismo ? item.organismo + " · " : ""}${item.accession}${
        item.tamanho_bases ? " · " + item.tamanho_bases.toLocaleString("pt-BR") + " bases" : ""}</small>
    `, (item) => {
    esconderResultados(elResultadosBact);
    ultimoOrganismoBact = item.organismo || item.titulo || "";
    elUsarImagemRef.checked = true;
    buscarComStatus(elBtnBuscarBact, elBuscaBactStatus,
      `${API}/api/buscar/acido_nucleico/detalhe?accession=${encodeURIComponent(item.accession)}&tipo=DNA`,
      (d) => {
        elBcTipo.value = d.tipo_acido;
        elBcSeq.value = d.sequencia_acido;
        elBuscaBactStatus.style.color = "#50c878";
        elBuscaBactStatus.textContent = `✅ ${d.fonte}`;
      });
  });
});

elBtnBuscarProt.addEventListener("click", async () => {
  const termo = elBuscaProtTermo.value.trim();
  if (!termo) { elBuscaProtStatus.textContent = "Digite um termo de busca."; return; }
  esconderResultados(elResultadosProt);
  const itens = await buscarListaComStatus(elBtnBuscarProt, elBuscaProtStatus,
    `${API}/api/buscar/proteina/lista?termo=${encodeURIComponent(termo)}&max_resultados=8`);
  if (!itens) return;
  elBuscaProtStatus.style.color = "#9aa0ac";
  elBuscaProtStatus.textContent = `${itens.length} resultado(s) no UniProt — escolha um abaixo.`;
  renderResultados(elResultadosProt, itens, (item) => `
      <b>${item.nome}</b>
      <small>${item.organismo ? item.organismo + " · " : ""}${item.accession}${
        item.tamanho_residuos ? " · " + item.tamanho_residuos + " aa" : ""}</small>
    `, (item) => {
    esconderResultados(elResultadosProt);
    buscarComStatus(elBtnBuscarProt, elBuscaProtStatus,
      `${API}/api/buscar/proteina/detalhe?accession=${encodeURIComponent(item.accession)}`,
      (d) => {
        elUsarComposicao.checked = true;
        elBlocoComposicao.style.display = "block";
        elUsarPool.checked = false;
        elBlocoPool.style.display = "none";
        document.querySelectorAll(".aa-qtd").forEach(inp => {
          inp.value = d.composicao_fixa[inp.dataset.aa] || 0;
        });
        elBuscaProtStatus.style.color = "#50c878";
        elBuscaProtStatus.textContent =
          `✅ ${d.nome}${d.organismo ? " (" + d.organismo + ")" : ""} — ${d.fonte}`;
      });
  });
});

// ── Helpers ──────────────────────────────────────────────────────────
function getComposicaoFixa() {
  const c = {};
  document.querySelectorAll(".aa-qtd").forEach(inp => {
    const v = parseInt(inp.value, 10);
    if (v > 0) c[inp.dataset.aa] = v;
  });
  return c;
}
function getPool() {
  return Array.from(document.querySelectorAll(".pool-check:checked")).map(c => c.value);
}
function getElementos() {
  return Array.from(elListaElementos.querySelectorAll("input:checked")).map(c => c.value);
}
function setStatus(t) { elStatusTexto.textContent = t; }
function apendarLog(ls) {
  elStatusLog.textContent = (ls || []).join("\n");
  elStatusLog.scrollTop = elStatusLog.scrollHeight;
}

// ── Simulação ────────────────────────────────────────────────────────
async function iniciarSimulacao() {
  const vcAtivo = elAtivarVc.checked;
  if (vcAtivo && elVcSeq.value.trim().length < 3) {
    alert("Digite uma sequência com pelo menos 3 bases, ou desmarque a opção de vírus personalizado.");
    return;
  }
  const incluirReais = elIncluirReais.checked;
  if (incluirReais && getElementos().length === 0) {
    alert("Selecione pelo menos um composto ou desmarque a opção de fármacos.");
    return;
  }

  const bactCustAtivo = elAtivarBc.checked;
  if (bactCustAtivo && elBcSeq.value.trim().length < 3) {
    alert("Digite uma sequência com pelo menos 3 bases para o genoma real da bactéria, ou desmarque a opção.");
    return;
  }

  let virusCustomizado = null;
  const preset = elSelectPreset.value;
  const usarImagemRef = elUsarImagemRef.checked;
  if (preset) {
    virusCustomizado = { preset, usar_imagem_referencia: usarImagemRef };
  } else if (vcAtivo) {
    virusCustomizado = {
      nome: elVcNome.value || "Meu-Virus",
      tipo_acido: elVcTipo.value,
      sequencia_acido: elVcSeq.value.trim(),
      usar_imagem_referencia: usarImagemRef,
      nome_organismo_imagem: ultimoOrganismoVirus,
    };
  }

  const tipoResposta = elTipoResposta?.value || "nenhuma";
  const limiteDigitado = Number(elLimiteEnergetico?.value);
  const limiteEnergetico = elUsarLimiteEnergetico?.checked && Number.isFinite(limiteDigitado)
    ? limiteDigitado : null;

  const corpo = {
    comprimento_proteina: parseInt(elComprimento.value, 10),
    n_proteinas:          parseInt(elNProteinas.value, 10),
    composicao_fixa:      elUsarComposicao.checked ? getComposicaoFixa() : {},
    aminoacidos_pool:     elUsarPool.checked ? getPool() : [],
    ambiente:             elSelectAmbiente.value,
    modo_rapido:          document.getElementById("modo-rapido").checked,
    limite_energetico:    limiteEnergetico,
    tipo_resposta:        tipoResposta,
    gerar_virus:          elGerarVirus.checked,
    virus_customizado:    virusCustomizado,
    gerar_bacterias:      elGerarBacterias.checked,
    n_bacterias:          parseInt(document.getElementById("n-bacterias").value, 10),
    bacteria_genoma_customizado: bactCustAtivo ? {
      tipo_acido: elBcTipo.value,
      sequencia_acido: elBcSeq.value.trim(),
      usar_imagem_referencia: usarImagemRef,
      nome_organismo_imagem: ultimoOrganismoBact,
    } : null,
    incluir_compostos_reais: incluirReais,
    elementos:            getElementos(),
    prever_mutacao:       elPreverMutacao.checked,
    n_mutacoes:           parseInt(elNMutacoes.value, 10),
  };

  elBtnSimular.disabled = true;
  elBtnSimular.textContent = "Simulando...";
  setStatus("Enviando simulação...");
  elSecaoResultados.style.display = "none";

  try {
    const r = await fetchAutenticado(`${API}/api/simular`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(corpo),
    });
    const { job_id } = await lerJson(r, "Iniciar simulação");
    await acompanharJob(job_id);
  } catch (e) {
    setStatus("Erro ao iniciar simulação.");
    apendarLog([String(e)]);
  } finally {
    elBtnSimular.disabled = false;
    elBtnSimular.textContent = "▶ Rodar simulação";
  }
}

async function acompanharJob(id) {
  setStatus("Simulando… (pode levar de segundos a poucos minutos)");
  while (true) {
    const r = await fetchAutenticado(`${API}/api/status/${id}`);
    const d = await lerJson(r, "Consultar status da simulação");
    apendarLog(d.log || []);
    if (d.estado === "concluido") { setStatus("Concluído ✅"); mostrarResultado(d.resultado); return; }
    if (d.estado === "erro")     { setStatus("Erro ❌"); apendarLog([...(d.log||[]), d.erro||""]); return; }
    await new Promise(res => setTimeout(res, 1500));
  }
}

// ── Exibir resultados ────────────────────────────────────────────────
function mostrarResultado(res) {
  elSecaoResultados.style.display = "block";

  // Estratégia de resposta selecionada
  if (elListaResposta) {
    const resp = res.resposta;
    if (!resp) {
      elListaResposta.innerHTML = "";
    } else {
      const linhasAlvos = (resp.alvos || []).map((a) => {
        const tipoAlvo = a.tipo === "bacteria" ? "Bactéria" : "Vírus";
        const energia = a.energia_melhor == null ? "—" : a.energia_melhor;
        const score = a.score_local == null ? "—" : a.score_local;
        const limite = a.dentro_limite === null || a.dentro_limite === undefined
          ? "" : (a.dentro_limite ? " · dentro ✅" : " · fora ⚠️");
        return `<tr><td>${tipoAlvo}: ${a.nome}</td><td>${a.ligante_melhor || "—"}</td><td>${energia}</td><td>${score}${limite}</td></tr>`;
      }).join("");
      const cobertura = resp.cobertura_percentual == null ? "sem corte" : `${resp.cobertura_percentual}%`;
      elListaResposta.innerHTML = `<div class="item-resultado">
        <h3>🎯 ${resp.label}</h3>
        <p>${resp.descricao}</p>
        <small><b>Critério:</b> ${resp.regra} · <b>Cobertura:</b> ${cobertura}</small>
        <p>${resp.resumo}</p>
        ${linhasAlvos ? `<table class="ranking"><tr><th>Alvo</th><th>Melhor candidato</th><th>Energia</th><th>Score local</th></tr>${linhasAlvos}</table>` : "<p>Nenhum alvo foi selecionado.</p>"}
        <small>⚠️ ${resp.observacao}</small>
      </div>`;
    }
  }

  // Proteínas
  const listaP = document.getElementById("lista-proteinas");
  listaP.innerHTML = "<h3>Proteínas geradas</h3>";
  (res.proteinas || []).forEach(p => {
    const s = p.estatisticas;
    const div = document.createElement("div");
    div.className = "item-resultado";
    div.innerHTML = `<b>${p.nome}</b>
      <br/><small>${s.n_residuos} resíduos &nbsp;·&nbsp;
        massa ≈ ${s.massa_total_da.toFixed(1)} Da &nbsp;·&nbsp;
        carga ${s.carga_liquida >= 0 ? "+" : ""}${s.carga_liquida.toFixed(2)}</small>
      ${p.energia_formacao != null ? `<br/><small>⚡ Energia de formação: <b>${p.energia_formacao}</b>${p.limite_energetico != null ? ` / limite ${p.limite_energetico} — ${p.atingiu_limite_energetico ? "atingido ✅" : "não atingido; melhor candidato ⚠️"}` : ""}</small>` : ""}
      <br/><small>
        <span style="color:#f05050">🔴 ${s.n_hidrofobo} hidrofóbicos</span> &nbsp;
        <span style="color:#50c878">🟢 ${s.n_hidrofilico} hidrofílicos</span>
      </small>
      ${p.sequencia_fasta ? `<br/><small class="fasta">${p.sequencia_fasta}</small>` : ""}
      ${p.glb ? `<br/><button data-glb="${p.glb}" class="btn-ver">Ver em 3D</button>` : ""}`;
    listaP.appendChild(div);
  });

  // Vírus
  const listaV = document.getElementById("lista-virus");
  listaV.innerHTML = "<h3>Vírus</h3>";
  (res.virus || []).forEach(v => {
    const div = document.createElement("div");
    div.className = "item-resultado";
    div.innerHTML = `${v.nome}
      ${v.customizado ? `<small>(genoma: ${v.genoma})</small>` : ""}
      ${v.imagem_referencia ? `
        <br/>
        <div class="ref-imagem">
          <div class="ref-titulo">📷 Comparação com referência científica</div>
          ${v.imagem_referencia.thumb_url ? `<img src="${v.imagem_referencia.thumb_url}" alt="Micrografia científica de referência" class="thumb-ref"/>` : ""}
          <small>Paleta/proporção do modelo ajustadas com <a href="${v.imagem_referencia.url}" target="_blank">micrografia filtrada (${v.imagem_referencia.titulo || "Wikimedia Commons"})</a>${v.imagem_referencia.licenca ? ` · ${v.imagem_referencia.licenca}` : ""}</small>
        </div>` : ""}
      ${v.glb ? `<button data-glb="${v.glb}" class="btn-ver">Ver em 3D</button>` : ""}`;
    listaV.appendChild(div);
  });

  // Bactérias
  const listaB = document.getElementById("lista-bacterias");
  listaB.innerHTML = "<h3>Bactérias</h3>";
  (res.bacterias || []).forEach(b => {
    const s = b.estatisticas;
    const div = document.createElement("div");
    div.className = "item-resultado";
    let detalhes = "";
    if (s) {
      const forma = s.forma === "coco" ? "esférica (coco)" : "bastonete (bacilo)";
      const tamanho = s.forma === "coco"
        ? `raio ≈ ${s.raio_corpo}`
        : `comprimento ≈ ${s.comprimento} · raio ≈ ${s.raio_corpo}`;
      detalhes = `<br/><small>forma ${forma} &nbsp;·&nbsp; ${tamanho}</small>
        <br/><small>
          🧫 ${s.n_flagelos} flagelo(s) &nbsp;·&nbsp;
          🔵 ${s.n_plasmideos} plasmídeo(s) &nbsp;·&nbsp;
          ⚪ ${s.n_ribossomos} ribossomos
        </small>
        <br/><small>proteína de superfície: ${s.n_residuos_proteina_superficie} resíduos
          ${s.genoma_customizado ? " &nbsp;·&nbsp; 🌐 genoma real" : ""}</small>
        ${b.imagem_referencia ? `
          <br/>
          <div class="ref-imagem">
            <div class="ref-titulo">📷 Comparação com referência científica</div>
            ${b.imagem_referencia.thumb_url ? `<img src="${b.imagem_referencia.thumb_url}" alt="Micrografia científica de referência" class="thumb-ref"/>` : ""}
            <small>Forma/paleta ajustadas com <a href="${b.imagem_referencia.url}" target="_blank">micrografia filtrada (${b.imagem_referencia.titulo || "Wikimedia Commons"})</a>${b.imagem_referencia.licenca ? ` · ${b.imagem_referencia.licenca}` : ""}</small>
          </div>` : ""}`;
    }
    div.innerHTML = `<b>${b.nome}</b>${detalhes}
      ${b.glb ? `<br/><button data-glb="${b.glb}" class="btn-ver">Ver em 3D</button>` : ""}`;
    listaB.appendChild(div);
  });

  // Docking
  const listaDock = document.getElementById("lista-docking");
  listaDock.innerHTML = "<h3>Melhor encaixe (docking)</h3>";
  (res.docking || []).forEach(d => {
    const div = document.createElement("div");
    div.className = "item-resultado";
    div.innerHTML = `<b>${d.composto}</b> × <b>${d.virus}</b> — score ${d.score}/100
      ${d.glb ? `<button data-glb="${d.glb}" class="btn-ver">Ver em 3D</button>` : ""}`;
    listaDock.appendChild(div);
  });

  // Previsão de mutações
  const listaMut = document.getElementById("lista-mutacoes");
  listaMut.innerHTML = "";
  if ((res.mutacoes && res.mutacoes.length) || res.mutacao_top) {
    listaMut.innerHTML = "<h3>🧬 Previsão de mutações (modelo ilustrativo)</h3>";

    (res.mutacoes || []).forEach(m => {
      const cor = m.classificacao.startsWith("radical") ? "#f05050"
                : m.classificacao.startsWith("moderada") ? "#e6b800" : "#50c878";
      const listaMuts = m.mutacoes.map(x => `${x.de_1letra}${x.posicao}${x.para_1letra}`).join(", ");
      const div = document.createElement("div");
      div.className = "item-resultado";
      div.innerHTML = `<b>${m.virus}</b>
        <br/><small class="fasta">${listaMuts}</small>
        <br/><small>Impacto médio: <span style="color:${cor}">${m.indice_impacto_medio}/100 — ${m.classificacao}</span></small>
        ${m.glb ? `<br/><button data-glb="${m.glb}" class="btn-ver">Ver variante mutante em 3D</button>` : ""}`;
      listaMut.appendChild(div);
    });

    if (res.mutacao_top) {
      const mt = res.mutacao_top;
      const div = document.createElement("div");
      div.className = "item-resultado";
      div.innerHTML = `<b>Impacto no melhor encaixe:</b> ${mt.composto} × ${mt.virus} (mutado)
        <br/><small>Energia de ligação: ${mt.energia_original} → ${mt.energia_mutante}
          (Δ${mt.delta_energia >= 0 ? "+" : ""}${mt.delta_energia})</small>
        <br/><small><b>${mt.interpretacao}</b></small>
        ${mt.glb ? `<br/><button data-glb="${mt.glb}" class="btn-ver">Ver docking mutante em 3D</button>` : ""}`;
      listaMut.appendChild(div);
    }
  }

  // Ranking
  const elRanking = document.getElementById("ranking");
  if (res.ranking && res.ranking.length) {
    let html = "<h3>Ranking</h3><table class='ranking'><tr><th>Proteína/Ligante</th><th>Alvo</th><th>Score</th></tr>";
    res.ranking.slice(0, 20).forEach(r => {
      html += `<tr><td>${r.composto}</td><td>${r.virus}</td><td>${r.score}</td></tr>`;
    });
    elRanking.innerHTML = html + "</table>";
  } else { elRanking.innerHTML = ""; }

  // Arquivos
  const elArq = document.getElementById("links-arquivos");
  elArq.innerHTML = "<h3>Arquivos</h3>";
  Object.entries(res.arquivos || {}).forEach(([nome, url]) => {
    if (!url) return;
    const a = document.createElement("a");
    a.href = url; a.target = "_blank"; a.className = "link-arquivo";
    a.textContent = `⬇ ${nome}`;
    elArq.appendChild(a);
  });

  // Botões "Ver em 3D"
  document.querySelectorAll(".btn-ver").forEach(btn =>
    btn.addEventListener("click", () => {
      elVisor.setAttribute("src", btn.dataset.glb);
      elVisorLegenda.style.display = "none";
    })
  );
  // Abre automaticamente a primeira proteína
  const primeiraGlb = (res.proteinas || []).find(p => p.glb)?.glb;
  if (primeiraGlb) { elVisor.setAttribute("src", primeiraGlb); elVisorLegenda.style.display = "none"; }
}

// ── Arranque ─────────────────────────────────────────────────────────
elBtnSimular.addEventListener("click", iniciarSimulacao);
carregarOpcoes();
