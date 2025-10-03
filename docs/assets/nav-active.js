(() => {
  // Normalize a path: remove trailing /index.html and trailing slash
  const norm = p => p.replace(/\/index\.html$/i,'').replace(/\/$/,'') || '/';

  const here = norm(new URL(location.href).pathname);

  let matched = null;
  document.querySelectorAll('header nav a').forEach(a => {
    const target = norm(new URL(a.getAttribute('href'), location.origin).pathname);
    // Exact match OR "we are inside that section" (e.g. /library vs /library/foo)
    const isActive = (here === target) || (target !== '/' && here.startsWith(target));
    if (isActive && !matched) matched = a; // prefer the first best match

    a.classList.remove('text-indigo-300','underline','underline-offset-4');
  });

  if (matched) {
    matched.classList.add('text-indigo-300','underline','underline-offset-4');
  }
})();


