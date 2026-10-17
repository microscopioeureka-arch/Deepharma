# Deepharma — contas, cotas e pagamentos (versão de trabalho)

Esta branch é uma implementação de teste. **A cobrança começa desativada e não deve ser ativada na hospedagem atual até que o banco tenha persistência real.** A landing page só libera o Checkout se os preços mensais do Stripe puderem ser validados.

## O que esta versão implementa

- Cadastro por e-mail com hash de senha PBKDF2-HMAC-SHA256; as senhas não são armazenadas em texto puro.
- Cota controlada pelo servidor e separada por conta.
- Plano Gratuito: 5 simulações por mês.
- Deeplatinum: cota provisória de 40 simulações nas últimas 5 horas.
- DeepharmaX: sem cota de quantidade. O backend atual ainda permite apenas uma simulação em execução no servidor por vez.
- Checkout hospedado pelo Stripe, portal de gerenciamento e verificação de assinatura de webhook com processamento idempotente.
- Preços que a interface mostra são lidos dos preços ativos e recorrentes no Stripe quando a integração está configurada.

Os limites e preços dos planos são uma **proposta de produto**, ainda precisam ser aprovados antes de qualquer lançamento. O modelo científico continua exploratório: as fitas proteicas e morfologias virais geradas são esquemáticas, não estruturas atômicas experimentalmente validadas.

## Configuração local / homologação

Mantenha o pagamento desligado enquanto desenvolve:

```text
ENABLE_BILLING=false
```

As variáveis existentes do servidor continuam necessárias:

- `APP_PASSWORD`: senha temporária do acesso legado de desenvolvimento.
- `SESSION_SECRET`: segredo aleatório longo para assinar os cookies de sessão.
- `COOKIE_SECURE=1` em HTTPS.

Para um ambiente Stripe de teste, configure no ambiente do servidor — nunca no JavaScript, HTML ou Git:

- `ENABLE_BILLING=true`
- `STRIPE_SECRET_KEY`: chave de teste do servidor.
- `STRIPE_WEBHOOK_SECRET`: segredo de assinatura do endpoint webhook.
- `STRIPE_PRICE_DEEPLATINUM`: ID de um preço ativo, recorrente, mensal (intervalo de 1 mês) para Deeplatinum.
- `STRIPE_PRICE_DEEPHARMAX`: ID de um preço ativo, recorrente, mensal (intervalo de 1 mês) para DeepharmaX.
- `DEEPHARMA_PUBLIC_URL`: origem pública completa, por exemplo `https://seu-dominio.example`, sem barra final.
- `DEEPHARMA_DB_PATH`: caminho absoluto para o arquivo SQLite em **armazenamento persistente**.

Os preços devem ser fixos e usar USD, BRL ou EUR neste protótipo. A página lê o valor dos objetos Price no Stripe, de modo que o botão só fica habilitado quando ambos os preços puderem ser verificados. O checkout também revalida o preço antes de criar a sessão.

## Webhook do Stripe

Cadastre no painel do provedor o endpoint:

```text
https://SEU-DOMINIO/api/billing/webhook
```

Inscreva o endpoint, no mínimo, nestes eventos:

- `checkout.session.completed`
- `customer.subscription.created`
- `customer.subscription.updated`
- `customer.subscription.deleted`
- `invoice.payment_failed`

Copie o segredo do endpoint para `STRIPE_WEBHOOK_SECRET`. Não use um webhook falso nem ative o plano só porque a pessoa voltou para a página de sucesso: o código depende do evento assinado pelo Stripe. Teste primeiro com chaves e preços de teste e confira cadastro, checkout cancelado, assinatura ativa, falha de pagamento, cancelamento e reentrega do mesmo evento.

## Persistência é um bloqueio para lançamento

O serviço deste repositório está configurado como Render Free em `render.yaml`. O filesystem de serviços Render é efêmero por padrão, e serviços gratuitos não aceitam discos persistentes. Se SQLite ficar no filesystem efêmero, contas e cotas podem desaparecer em reinícios/deploys; **não habilite cobrança com essa configuração**. Antes do lançamento, conecte um banco persistente suportado ou mude para um serviço com disco persistente e configure `DEEPHARMA_DB_PATH` dentro do ponto de montagem. A documentação do Render confirma essas limitações: https://render.com/docs/free e https://render.com/docs/disks.

Também falta implementar e testar verificação de e-mail, recuperação de senha, proteção contra tentativas repetidas de login, política de privacidade/termos e configuração de backup/retention. O acesso legado por senha compartilhada é apenas transitório e deve ser removido antes de abrir o serviço ao público.

A conta de cobrança do Stripe e os contratos/preços devem ser configurados pelo titular autorizado, com o apoio de um responsável adulto quando aplicável. Nunca compartilhe chaves secretas em commits, prints ou mensagens.

## Checklist antes de merge/deploy

1. CI verde para compilação e testes.
2. Testes manuais de criação de conta, login/logout, cota gratuita e cota Deeplatinum no banco de teste.
3. Verificar que duas contas diferentes não conseguem consultar o job uma da outra.
4. Validar preços reais (valor, moeda, recorrência mensal) na página antes de habilitar Checkout.
5. Testar todos os eventos do webhook no modo de teste e confirmar que assinatura só muda após evento válido.
6. Resolver o armazenamento persistente e remover o login legado antes de habilitar a cobrança real.
