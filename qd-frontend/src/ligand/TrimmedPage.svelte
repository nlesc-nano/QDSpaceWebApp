<script>
  // A pipeline report page (props/*.html, same origin) showing only some of its <h2> sections.
  // Pages written before ground_state.html and stability.html were split hold every section;
  // keep(title) picks the ones to show. autosize: the frame grows to its content (no inner scroll).
  let { src, keep = () => true, header = true, autosize = false, title = "Report page", cls = "" } = $props();

  let frame = $state(null);
  let height = $state(600);
  let observer = null;

  function trim() {
    const doc = frame?.contentDocument;
    if (!doc) return;
    const wrap = doc.querySelector(".wrap");
    if (wrap) {
      let section = null;
      for (const el of wrap.children) {
        if (el.tagName === "H2") section = el.textContent.trim();
        const show = section === null ? header : keep(section);
        el.style.display = show ? "" : "none";
      }
      frame.contentWindow.dispatchEvent(new Event("resize"));    // Plotly panels that were hidden
    }
    if (autosize) {
      const fit = () => (height = Math.max(200, doc.documentElement.scrollHeight + 4));
      fit();
      observer?.disconnect();
      observer = new ResizeObserver(fit);
      observer.observe(doc.body);
    }
  }
  $effect(() => () => observer?.disconnect());
</script>

<iframe bind:this={frame} {title} {src} onload={trim}
        class="w-full rounded-[1rem] border border-slate-200 bg-white {cls}"
        style={autosize ? `height: ${height}px` : ""}
        sandbox="allow-scripts allow-same-origin allow-popups allow-modals" referrerpolicy="no-referrer"></iframe>
