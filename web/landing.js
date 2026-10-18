(() => {
  const toggle = document.querySelector('.menu-toggle');
  const nav = document.querySelector('.nav-links');

  const closeMenu = () => {
    if (!toggle || !nav) return;
    nav.classList.remove('open');
    toggle.setAttribute('aria-expanded', 'false');
    toggle.setAttribute('aria-label', 'Abrir menu');
    toggle.textContent = '☰';
  };

  if (toggle && nav) {
    toggle.addEventListener('click', () => {
      const open = nav.classList.toggle('open');
      toggle.setAttribute('aria-expanded', String(open));
      toggle.setAttribute('aria-label', open ? 'Fechar menu' : 'Abrir menu');
      toggle.textContent = open ? '×' : '☰';
    });

    nav.querySelectorAll('a').forEach(link => link.addEventListener('click', closeMenu));

    document.addEventListener('keydown', event => {
      if (event.key === 'Escape') closeMenu();
    });

    window.addEventListener('resize', () => {
      if (window.innerWidth > 800) closeMenu();
    });
  }

  const year = document.getElementById('year');
  if (year) year.textContent = new Date().getFullYear();

  const billingStatus = document.getElementById('billing-status');
  const checkoutButtons = [...document.querySelectorAll('[data-checkout-plan]')];
  let billingEnabled = false;

  async function readJson(response, context) {
    let data;
    try {
      data = await response.json();
    } catch (_) {
      throw new Error(context + ': resposta inválida do servidor.');
    }
    if (!response.ok) {
      const detail = typeof data.detail === 'string'
        ? data.detail
        : (data.detail && data.detail.message) || 'A solicitação não pôde ser concluída.';
      throw new Error(detail);
    }
    return data;
  }

  async function startCheckout(plan, button) {
    if (!['deeplatinum', 'deepharmax'].includes(plan)) return;
    if (!billingEnabled) {
      if (billingStatus) billingStatus.textContent = 'Os pagamentos ainda não estão configurados.';
      return;
    }

    if (button) {
      button.disabled = true;
      button.textContent = 'Conectando ao checkout...';
    }
    if (billingStatus) billingStatus.textContent = 'Preparando o checkout seguro...';

    try {
      const csrfResponse = await fetch('/api/csrf');
      const csrfData = await readJson(csrfResponse, 'Token de segurança');
      const response = await fetch('/api/billing/checkout', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'X-CSRF-Token': csrfData.csrf_token
        },
        body: JSON.stringify({ plan })
      });

      if (response.status === 401) {
        window.location.href = '/login?checkout=' + encodeURIComponent(plan);
        return;
      }

      const data = await readJson(response, 'Criar checkout');
      if (!data.url || !data.url.startsWith('https://checkout.stripe.com/')) {
        throw new Error('O provedor não retornou um link seguro de checkout.');
      }
      window.location.assign(data.url);
    } catch (error) {
      if (billingStatus) billingStatus.textContent = error.message || 'Não foi possível iniciar o checkout.';
    } finally {
      if (button && window.location.pathname === '/') {
        button.disabled = !billingEnabled;
        button.textContent = 'Assinar quando disponível';
      }
    }
  }

  async function loadBillingStatus() {
    if (!checkoutButtons.length) return;
    try {
      const response = await fetch('/api/billing/status');
      const data = await readJson(response, 'Verificar pagamentos');
      billingEnabled = data.enabled === true;
      if (billingEnabled && data.plans) {
        Object.entries(data.plans).forEach(([plan, price]) => {
          const priceElement = document.querySelector(`[data-plan-price="${plan}"]`);
          if (!priceElement || !Number.isInteger(price.amount) || !price.currency) return;
          priceElement.textContent = new Intl.NumberFormat('pt-BR', {
            style: 'currency',
            currency: String(price.currency).toUpperCase()
          }).format(price.amount / 100);
          const interval = document.createElement('small');
          interval.textContent = price.interval === 'month' ? ' / mês' : ' / ' + price.interval;
          priceElement.append(' ', interval);
        });
      }
    } catch (_) {
      billingEnabled = false;
    }

    checkoutButtons.forEach(button => {
      button.disabled = !billingEnabled;
      button.textContent = billingEnabled ? 'Assinar' : 'Pagamento indisponível';
      button.addEventListener('click', () => {
        startCheckout(button.dataset.checkoutPlan, button);
      });
    });

    if (billingStatus) {
      billingStatus.textContent = billingEnabled
        ? 'Checkout seguro configurado. Entre ou crie uma conta individual para assinar.'
        : 'Planos pagos desativados: o provedor, os preços e o webhook ainda precisam ser configurados.';
    }

    const query = new URLSearchParams(window.location.search);
    const pendingPlan = query.get('checkout');
    if (billingEnabled && ['deeplatinum', 'deepharmax'].includes(pendingPlan)) {
      window.history.replaceState(null, '', window.location.pathname + window.location.hash);
      await startCheckout(pendingPlan, checkoutButtons.find(button => button.dataset.checkoutPlan === pendingPlan));
    }
  }

  void loadBillingStatus();
})();
