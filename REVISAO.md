# Revisão do Deepharma — 14/10/2026

## Resultado e limites

Esta entrega corrige problemas de software e acrescenta uma auditoria numérica e
de procedência. **Não converte o núcleo atual em simulador biologicamente validado.**
Nenhuma alteração foi publicada automaticamente no GitHub ou na VPS.

O código público de `main` foi baixado nesta data. A revisão é direcionada, não
uma auditoria exaustiva de segurança, de conformidade ou de todas as rotinas físicas.

## Achados confirmados

1. `simulacao.py` constrói geometrias e sequências ilustrativas; o modo lento não
   implementa dinâmica molecular atomística validada. As espículas não são
   reconstruídas da sequência genômica importada.
2. O ranking reescala posições para um score de 0 a 100. Esse número depende dos
   candidatos presentes e não é probabilidade, afinidade ou eficácia.
3. A busca UniProt é usada como composição de aminoácidos: a ordem original não é
   preservada na geração. Portanto não simula necessariamente a proteína buscada.
4. A geração de variantes é heurística; não demonstra evolução ou efeito biológico.
   Não foi aprimorada nesta revisão e permanece sem validade preditiva.
5. Uma energia NaN/Infinity era substituída por 1e12, ocultando falha numérica.
   Agora a execução falha explicitamente nesse ponto.
6. Os pedidos não limitavam composição total nem várias dimensões de entrada.
   Agora tamanhos, composição e opções principais são validados antes da execução.
7. Cada pedido criava uma thread de simulação sem limite. Agora há uma execução
   por processo e no máximo 100 registros de jobs em memória. Erros ao criar a
   thread liberam a vaga. Os arquivos em disco não são apagados automaticamente.
8. Exceções expunham traceback ao cliente; os detalhes agora ficam no log do
   servidor. O usuário recebe mensagem genérica com identificação do job.
9. Comparação de senhas Unicode podia levantar erro. A comparação agora usa UTF-8.
10. Nomes de usuários, resultados externos e rótulos eram interpolados como HTML.
    Os principais campos dinâmicos foram escapados e URLs de visualização filtradas
    para HTTP/HTTPS. Isso não constitui certificação completa contra XSS.

## Auditoria disponível em cada execução web

`qualidade.json` registra configuração, horário UTC, versões, plataforma e SHA-256
do código. A resposta da API inclui `qualidade_cientifica` e link para download.
Valores finais não finitos impedem a publicação do resultado como concluído.
O relatório continua dizendo `validade_biologica: nao_demonstrada`, mesmo quando
todos os valores são finitos. A interface mostra isso antes e depois da execução.

O relatório não registra todo dado remoto recuperado nem garante reprodução bit
a bit. A checagem não substitui inspeção geométrica ou validação de trajetórias.

## O que falta para pesquisa de alta qualidade

É necessário separar o demonstrador atual de um futuro motor científico validado.
Para benchmarks iniciais, usar sistemas de referência não patogênicos, com dados
experimentais e condições documentadas. Não extrapolar esses benchmarks para vírus.

| Verificação necessária | Evidência que deve existir | Estado atual |
|---|---|---|
| Identidade da estrutura | Fonte, cadeia, sequência, cobertura, resíduos ausentes, montagem biológica e relatório wwPDB | Não implementada |
| Integridade numérica | Rejeição de NaN/Infinity nos resultados | Implementada parcialmente; intermediários não são todos inspecionados |
| Geometria | Ligações, ângulos, colisões e qualidade estrutural com ferramentas reconhecidas | Não implementada |
| Modelo físico | Campo de força, unidades, solvente e parâmetros registrados e compatíveis | Não implementado |
| Reprodutibilidade | Entradas imutáveis, versões fixadas, sementes controladas e execuções repetidas | Apenas procedência de código/configuração adicionada |
| Estabilidade e convergência | Trajetórias, controles e análise de incerteza de réplicas independentes | Não implementadas |
| Concordância experimental | Métrica e tolerância definidas antes do teste, dados independentes e relatório de erros | Não demonstrada |

Um benchmark precisa poder reprovar o modelo. Não basta exibir RMSD ou uma energia
baixa; os observáveis precisam responder à pergunta científica e ser comparados
com referências apropriadas. Critérios de aceitação dependem do método e do sistema.

Referências primárias para essa etapa:

- [wwPDB — relatórios de validação estrutural](https://www.wwpdb.org/validation/validation-reports)
- [GROMACS — documentação, validação e reprodutibilidade](https://manual.gromacs.org/current/user-guide/)

## Recursos necessários

Não é possível dimensionar a máquina científica sem número de átomos, método,
duração das trajetórias e volume de tarefas. GPU por si só não confere validade.

- Manter a VPS como interface/API e usar worker separado para cálculos longos.
- Para um motor reconhecido como GROMACS ou OpenMM, verificar compatibilidade de
  hardware/software e medir um benchmark representativo antes de contratar GPU.
- Medir tempo por execução, pico de RAM, memória da GPU e tamanho dos resultados.
  Usar essas medidas para limitar fila, concorrência e armazenamento.
- Guardar entradas, estruturas, versões e resultados em armazenamento persistente,
  com política explícita de retenção. O sistema atual usa jobs em memória.
- Fixar dependências num ambiente validado. Nesta entrega apenas Pydantic ganhou
  faixa de versão exigida pelo código; não foi criado lockfile do motor completo.

## Pendências de produto e operação

- O repositório continua usando uma senha compartilhada. Não há cadastro individual,
  recuperação de senha ou isolamento por usuário para jobs/arquivos.
- Faltam limitação de tentativas de login, fila persistente e quotas de disco.
- As buscas externas ainda precisam de limites, cache e tratamento de indisponibilidade
  uniformes. Esta revisão não certificou a seleção taxonômica dos resultados.
- O serviço deve continuar com **um worker**: a trava adicionada é local ao processo.
- Para HTTPS, configurar `COOKIE_SECURE=1`, conforme o guia original.

## Testes e implantação

Testes HTTP usam um substituto do núcleo molecular. Eles verificam contrato,
autenticação, limites, concorrência e publicação do relatório; não medem ciência.
Não foi executada a simulação completa com Torch/RDKit nem validada a VPS real.

Em ambiente de teste com as dependências instaladas:

```text
python -m pip install httpx
python -m unittest test_quality test_server_contract -v
node --check web/app.js
```

O arquivo `ALTERACOES.patch` permite revisar as diferenças. O ZIP contém apenas
arquivos novos/modificados, com caminhos relativos à raiz do repositório. Não
substitua o diretório inteiro por esse pacote: aplique os arquivos sobre uma cópia
do projeto, revise as diferenças, teste num ambiente de homologação e faça backup
antes de atualizar o serviço existente. `requirements.txt` deve ser atualizado.

Valide login, abertura do simulador, uma execução e download de `qualidade.json` em
homologação antes de promover. Guarde a versão anterior para reversão. Os comandos
de reinício dependem do serviço real da VPS, que não foi identificado nesta revisão.
