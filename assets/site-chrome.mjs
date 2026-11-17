// The viewer is published at sahmed73.github.io/molecule-viewer/ from its own repository.
// To match the rest of the website, it borrows the site's current header and footer from
// /assets/site-chrome.json, which the website builds from its own templates. If that file
// is missing (for example in a local copy), the simple header and footer in index.html stay.

const chromeUrl = '/assets/site-chrome.json';

function replaceWith(selector, html) {
  const current = document.querySelector(selector);
  if (!current || !html) return;
  const template = document.createElement('template');
  template.innerHTML = html.trim();
  const next = template.content.firstElementChild;
  if (next?.matches(selector)) current.replaceWith(next);
}

try {
  const response = await fetch(chromeUrl, {cache: 'no-cache'});
  if (response.ok) {
    const {header, footer} = await response.json();
    replaceWith('header.site-header', header);
    replaceWith('footer.site-footer', footer);
    const here = location.pathname.replace(/index\.html$/, '');
    for (const link of document.querySelectorAll('.site-nav a')) {
      if (new URL(link.href, location.href).pathname === here) link.setAttribute('aria-current', 'page');
      else link.removeAttribute('aria-current');
    }
  }
} catch {
  // Keep the built-in header and footer.
}
