// fgls site: the theme switch of emiliogcs.com in the footer, and the
// contents as a closed "contents +" above the text on narrow screens, as on
// emiliogcs.com/code/. Runs after pkgdown.js has built the contents (#toc).

window.jQuery(() => {
  // ---------------------------------------------------------------- theme
  // auto follows the system; light or dark is remembered in localStorage
  // under the same key as on the rest of emiliogcs.com, and applied before
  // first paint by the inline script in <head> (_pkgdown.yml).
  const row = document.querySelector('footer .footer-row');
  if (row) {
    const root = document.documentElement;
    const metas = document.querySelectorAll('meta[name="theme-color"]');
    const paper = { light: '#F5F3EE', dark: '#0F0F0E' };

    const switcher = document.createElement('div');
    switcher.className = 'theme-switch';
    switcher.setAttribute('role', 'group');
    switcher.setAttribute('aria-labelledby', 'theme-label');
    const label = document.createElement('span');
    label.className = 'theme-label';
    label.id = 'theme-label';
    label.textContent = 'theme';
    switcher.append(label);
    for (const mode of ['auto', 'light', 'dark']) {
      const b = document.createElement('button');
      b.type = 'button';
      b.value = mode;
      b.textContent = mode;
      switcher.append(b);
    }
    row.append(switcher);

    const apply = (mode) => {
      if (mode === 'light' || mode === 'dark') root.dataset.theme = mode;
      else delete root.dataset.theme;
      for (const m of metas) {
        const own = (m.getAttribute('media') || '').includes('dark') ? 'dark' : 'light';
        m.setAttribute('content', paper[root.dataset.theme || own]);
      }
      for (const b of switcher.querySelectorAll('button')) {
        b.setAttribute('aria-pressed', String(b.value === (root.dataset.theme || 'auto')));
      }
    };

    apply(root.dataset.theme || 'auto');
    switcher.addEventListener('click', (e) => {
      const b = e.target.closest('button');
      if (!b) return;
      apply(b.value);
      try {
        if (b.value === 'auto') localStorage.removeItem('theme');
        else localStorage.setItem('theme', b.value);
      } catch {}
    });
  }

  // ------------------------------------------------------------- contents
  // Below 1200px the contents column is hidden (extra.css); a copy of it
  // opens above the text instead, after the page's title.
  const toc = document.querySelector('#toc > .nav');
  const header = document.querySelector('main > .page-header');
  if (toc && header && toc.querySelector('a')) {
    const details = document.createElement('details');
    details.className = 'toc-top';
    const summary = document.createElement('summary');
    summary.innerHTML = 'contents <span class="sign" aria-hidden="true"></span>';
    const nav = document.createElement('nav');
    nav.setAttribute('aria-label', 'Contents of this page');
    const list = toc.cloneNode(true);
    for (const a of list.querySelectorAll('.active')) a.classList.remove('active');
    nav.append(list);
    details.append(summary, nav);
    header.after(details);
  }
});
