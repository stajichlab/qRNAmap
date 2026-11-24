// Copied from qmutfind 4404e1c93f47ea527287b93d6b496070563b193a, site/theme.js. Changes: see site/COPIED.md.
// Theme buttons: set the theme, save the choice, and mark the active button.
(function () {
  var root = document.documentElement;
  function mark() {
    document.querySelectorAll('[data-set-theme]').forEach(function (b) {
      b.setAttribute('aria-pressed', String(b.getAttribute('data-set-theme') === root.getAttribute('data-theme')));
    });
  }
  document.querySelectorAll('[data-set-theme]').forEach(function (button) {
    button.addEventListener('click', function () {
      var theme = button.getAttribute('data-set-theme');
      root.setAttribute('data-theme', theme);
      try { localStorage.setItem('qrnamap-theme', theme); } catch (e) {}
      mark();
      document.dispatchEvent(new Event('qrnamap-theme'));
    });
  });
  mark();

  // Layout buttons: standard, wide (the whole window) and phone (the stacked skin). A click is a choice: it is saved, and the
  // layout then stops following the width of the window. Without a choice the layout follows the window (phone when narrow).
  function markLayout() {
    document.querySelectorAll('[data-set-layout]').forEach(function (b) {
      b.setAttribute('aria-pressed', String(b.getAttribute('data-set-layout') === root.getAttribute('data-layout')));
    });
  }
  document.querySelectorAll('[data-set-layout]').forEach(function (button) {
    button.addEventListener('click', function () {
      var layout = button.getAttribute('data-set-layout');
      root.setAttribute('data-layout', layout);
      root.setAttribute('data-layout-auto', 'no');
      try { localStorage.setItem('qrnamap-layout', layout); } catch (e) {}
      markLayout();
    });
  });
  var narrow = window.matchMedia('(max-width: 720px)');
  function follow() {
    if (root.getAttribute('data-layout-auto') === 'yes') { root.setAttribute('data-layout', narrow.matches ? 'phone' : 'standard'); markLayout(); }
  }
  if (narrow.addEventListener) narrow.addEventListener('change', follow); else narrow.addListener(follow);
  markLayout();
})();
