<script>
  // Optical properties of one dot from its QDEX webapp.json (dft/qdex/):
  // absorption spectrum, band-edge excitons and the gap summary, spin free and SOC.
  let { url } = $props();

  let data = $state(null);
  let status = $state("idle");
  let kind = $state("sf");        // "sf" | "soc": spectrum, markers and table
  let compare = $state(false);    // overlay the other kind's spectrum (dashed)
  let sigma = $state("");         // broadening key, e.g. "sigma_0.030"
  let edgeOnly = $state(false);   // x range: band edge only (lowest exciton + 1.5 eV) or the whole grid
  let hoverE = $state(null);

  const KINDS = { sf: { label: "Spin free", summary: "spin_free", color: "#4f46e5" },
                  soc: { label: "SOC", summary: "soc", color: "#e11d48" } };
  const other = (k) => (k === "sf" ? "soc" : "sf");

  $effect(() => {
    const u = url;
    if (!u) return;
    status = "loading";
    data = null;
    fetch(u)
      .then((res) => { if (!res.ok) throw new Error(`HTTP ${res.status}`); return res.json(); })
      .then((d) => {
        if (u !== url) return;
        data = d;
        const keys = sigmaKeys(d, "sf");
        sigma = keys.includes("sigma_0.030") ? "sigma_0.030" : keys[0] || "";
        status = "ready";
      })
      .catch(() => { if (u === url) status = "error"; });
  });

  const sigmaKeys = (d, k) => Object.keys(d?.spectrum?.[k] || {}).filter((s) => s.startsWith("sigma_")).sort();
  const sigmaLabel = (s) => `${Number(s.slice(6)).toFixed(2)} eV`;
  const fmt = (v, n = 3) => (v === undefined || v === null ? "–" : Number(v).toFixed(n));

  // --- Plot geometry (viewBox units) ---
  const W = 720, H = 300, M = { l: 52, r: 16, t: 26, b: 40 };
  const PW = W - M.l - M.r, PH = H - M.t - M.b;

  // [energy, intensity] pairs of one spectrum.
  function series(k) {
    const sp = data?.spectrum?.[k];
    if (!sp || !sp[sigma]) return [];
    const [e0, , de] = sp.energy_ev;
    return sp[sigma].map((v, i) => [e0 + i * de, v]);
  }
  let shown = $derived(data ? (compare ? [kind, other(kind)] : [kind]).filter((k) => data.spectrum?.[k]) : []);
  // x range: from just below the lowest exciton to emax_ev (or 1.5 eV above it). QDEX keeps the
  // excitations up to emax_ev only, so the convolved spectrum above it is incomplete.
  let xRange = $derived.by(() => {
    if (!data) return [0, 1];
    const [e0, grid1] = data.spectrum?.[kind]?.energy_ev || [0, 1];
    const emax = shown.map((k) => data.summary?.[KINDS[k].summary]?.emax_ev).filter((v) => v !== undefined);
    const e1 = Math.min(grid1, ...emax);
    const lows = shown.map((k) => data.summary?.[KINDS[k].summary]?.lowest_ev).filter((v) => v !== undefined);
    if (!lows.length) return [e0, e1];
    const lo = Math.max(e0, Math.floor((Math.min(...lows) - 0.4) * 2) / 2);
    return [lo, edgeOnly ? Math.min(e1, Math.ceil((Math.min(...lows) + 1.5) * 2) / 2) : e1];
  });
  let curves = $derived(shown.map((k) => ({ k, pts: series(k).filter(([e]) => e >= xRange[0] && e <= xRange[1]) })));
  let yMax = $derived(Math.max(1e-9, ...curves.flatMap((c) => c.pts.map(([, v]) => v))) * 1.05);
  const sx = (e) => M.l + ((e - xRange[0]) / (xRange[1] - xRange[0])) * PW;
  const sy = (v) => M.t + PH - (v / yMax) * PH;
  const pathOf = (pts) => pts.map(([e, v], i) => `${i ? "L" : "M"}${sx(e).toFixed(1)},${sy(v).toFixed(1)}`).join("");

  let xTicks = $derived.by(() => {
    const out = [];
    for (let e = Math.ceil(xRange[0] * 2) / 2; e <= xRange[1] + 1e-9; e += 0.5) out.push(e);
    return out;
  });
  let yTicks = $derived.by(() => {
    const raw = yMax / 4, mag = 10 ** Math.floor(Math.log10(raw));
    const step = [1, 2, 5, 10].map((m) => m * mag).find((s) => s >= raw);
    const out = [];
    for (let v = 0; v <= yMax; v += step) out.push(v);
    return out;
  });

  let summary = $derived(data?.summary?.[KINDS[kind].summary] || {});
  let excitons = $derived(data?.excitons?.[kind] || []);
  let brightF = $derived(data?.summary?.bright_f ?? 0.05);
  // Oscillator-strength sticks of the band-edge excitons, up to a third of the plot height.
  let fMax = $derived(Math.max(1e-9, ...excitons.map((x) => x.f_osc || 0)));
  let markers = $derived.by(() => {
    const s = summary, out = [];
    if (s.lowest_ev !== undefined && s.first_bright_ev !== undefined && Math.abs(s.lowest_ev - s.first_bright_ev) < 1e-4) {
      out.push({ e: s.lowest_ev, label: "lowest = first bright", dash: "" });
    } else {
      if (s.lowest_ev !== undefined) out.push({ e: s.lowest_ev, label: "lowest", dash: "3 3" });
      if (s.first_bright_ev !== undefined) out.push({ e: s.first_bright_ev, label: "first bright", dash: "" });
    }
    return out.filter((m) => m.e >= xRange[0] && m.e <= xRange[1]);
  });

  function onMove(ev) {
    const r = ev.currentTarget.getBoundingClientRect();
    const x = ((ev.clientX - r.left) / r.width) * W;
    hoverE = x < M.l || x > W - M.r ? null : xRange[0] + ((x - M.l) / PW) * (xRange[1] - xRange[0]);
  }
  let hover = $derived.by(() => {
    if (hoverE === null) return null;
    const rows = curves.map((c) => {
      let best = c.pts[0];
      for (const p of c.pts) if (Math.abs(p[0] - hoverE) < Math.abs(best[0] - hoverE)) best = p;
      return { k: c.k, e: best?.[0], v: best?.[1] };
    }).filter((r) => r.e !== undefined);
    return rows.length ? { e: rows[0].e, rows } : null;
  });
</script>

{#if status === "loading"}
  <div class="flex items-center justify-center h-full text-brand-600 text-sm font-bold">
    <div class="w-6 h-6 border-4 border-brand-500 border-t-transparent rounded-full animate-spin mr-3"></div>Loading optical data…
  </div>
{:else if status === "error"}
  <div class="flex items-center justify-center h-full text-slate-500 text-sm">Could not load {url}.</div>
{:else if data}
  <div class="flex flex-col gap-5 pb-4">
    <!-- Controls -->
    <div class="flex flex-wrap items-center gap-3 text-xs">
      <div class="flex gap-1 bg-slate-100 p-1 rounded-xl border border-slate-200">
        {#each Object.entries(KINDS) as [k, meta]}
          {#if data.spectrum?.[k] || data.excitons?.[k]}
            <button class="px-3 py-1.5 rounded-lg font-bold transition-all {kind === k ? 'bg-white text-slate-900 shadow-sm' : 'text-slate-500 hover:text-slate-800'}"
                    onclick={() => (kind = k)}>{meta.label}</button>
          {/if}
        {/each}
      </div>
      <label class="flex items-center gap-1.5 cursor-pointer text-slate-600 font-medium">
        <input type="checkbox" bind:checked={compare} class="accent-accent-600" /> overlay {KINDS[other(kind)].label}
      </label>
      <span class="font-bold text-slate-400 uppercase tracking-widest ml-2">Range</span>
      <div class="flex gap-1 bg-slate-100 p-1 rounded-xl border border-slate-200">
        {#each [[false, "Full"], [true, "Band edge"]] as [v, label]}
          <button class="px-3 py-1.5 rounded-lg font-bold transition-all {edgeOnly === v ? 'bg-white text-slate-900 shadow-sm' : 'text-slate-500 hover:text-slate-800'}"
                  onclick={() => (edgeOnly = v)}>{label}</button>
        {/each}
      </div>
      <span class="font-bold text-slate-400 uppercase tracking-widest ml-2">Broadening</span>
      <div class="flex gap-1 bg-slate-100 p-1 rounded-xl border border-slate-200">
        {#each sigmaKeys(data, kind) as s}
          <button class="px-3 py-1.5 rounded-lg font-bold transition-all {sigma === s ? 'bg-white text-slate-900 shadow-sm' : 'text-slate-500 hover:text-slate-800'}"
                  onclick={() => (sigma = s)}>{sigmaLabel(s)}</button>
        {/each}
      </div>
    </div>

    <div class="grid grid-cols-1 xl:grid-cols-3 gap-5">
      <!-- Spectrum -->
      <div class="xl:col-span-2 relative">
        <h3 class="text-xs font-bold text-slate-400 uppercase tracking-widest mb-1">Absorption spectrum</h3>
        <svg viewBox="0 0 {W} {H}" class="w-full h-auto select-none" role="img"
             aria-label="Absorption spectrum, {KINDS[kind].label}" onpointermove={onMove} onpointerleave={() => (hoverE = null)}>
          {#each yTicks as v}
            <line x1={M.l} x2={W - M.r} y1={sy(v)} y2={sy(v)} stroke="#e2e8f0" stroke-width="1" />
            <text x={M.l - 6} y={sy(v) + 3} text-anchor="end" font-size="10" fill="#94a3b8">{v.toPrecision(3).replace(/\.?0+$/, "")}</text>
          {/each}
          {#each xTicks as e}
            <line x1={sx(e)} x2={sx(e)} y1={M.t + PH} y2={M.t + PH + 4} stroke="#94a3b8" />
            <text x={sx(e)} y={M.t + PH + 16} text-anchor="middle" font-size="10" fill="#64748b">{e.toFixed(1)}</text>
          {/each}
          <line x1={M.l} x2={W - M.r} y1={M.t + PH} y2={M.t + PH} stroke="#94a3b8" />
          <text x={M.l + PW / 2} y={H - 4} text-anchor="middle" font-size="11" fill="#475569">Energy (eV)</text>
          <text transform="translate(12 {M.t + PH / 2}) rotate(-90)" text-anchor="middle" font-size="11" fill="#475569">Absorption (arb. u.)</text>

          {#each excitons as x}
            {#if x.energy_ev >= xRange[0] && x.energy_ev <= xRange[1] && x.f_osc > 0}
              <line x1={sx(x.energy_ev)} x2={sx(x.energy_ev)} y1={M.t + PH} y2={M.t + PH - (x.f_osc / fMax) * PH * 0.33}
                    stroke="#94a3b8" stroke-width="1.5"><title>state {x.state} · {fmt(x.energy_ev)} eV · f = {fmt(x.f_osc, 4)} · {x.transition}</title></line>
            {/if}
          {/each}
          {#each markers as m}
            <line x1={sx(m.e)} x2={sx(m.e)} y1={M.t} y2={M.t + PH} stroke="#334155" stroke-width="1" stroke-dasharray={m.dash} />
            <text x={sx(m.e) + 4} y={M.t - 8} font-size="10" fill="#334155">{m.label} · {fmt(m.e)} eV</text>
          {/each}
          {#each curves as c}
            <path d={pathOf(c.pts)} fill="none" stroke={KINDS[c.k].color} stroke-width="2"
                  stroke-dasharray={c.k === kind ? "" : "5 4"} stroke-linejoin="round" />
          {/each}
          {#if hover}
            <line x1={sx(hover.e)} x2={sx(hover.e)} y1={M.t} y2={M.t + PH} stroke="#cbd5e1" stroke-width="1" />
            {#each hover.rows as r}
              <circle cx={sx(r.e)} cy={sy(r.v)} r="4" fill={KINDS[r.k].color} stroke="white" stroke-width="2" />
            {/each}
          {/if}
        </svg>
        {#if hover}
          <div class="absolute pointer-events-none bg-white/95 border border-slate-200 rounded-lg shadow-soft px-2.5 py-1.5 text-[11px] font-mono text-slate-700"
               style="left: {(sx(hover.e) / W) * 100}%; top: 1.75rem; transform: translateX({sx(hover.e) > W * 0.7 ? 'calc(-100% - 8px)' : '8px'});">
            <div class="font-bold text-slate-900">{fmt(hover.e)} eV</div>
            {#each hover.rows as r}
              <div class="flex items-center gap-1.5"><span class="w-2.5 h-0.5 inline-block" style="background: {KINDS[r.k].color}"></span>{KINDS[r.k].label} {fmt(r.v, 2)}</div>
            {/each}
          </div>
        {/if}
        <div class="flex flex-wrap gap-x-4 gap-y-1 text-[11px] text-slate-500 mt-1 px-1">
          {#each curves as c}
            <span class="flex items-center gap-1.5"><svg width="18" height="6"><line x1="0" x2="18" y1="3" y2="3" stroke={KINDS[c.k].color} stroke-width="2" stroke-dasharray={c.k === kind ? "" : "5 4"} /></svg>{KINDS[c.k].label}, σ = {sigmaLabel(sigma)}</span>
          {/each}
          <span class="flex items-center gap-1.5"><span class="w-0.5 h-2.5 bg-slate-400 inline-block"></span>band-edge excitons, height ∝ f</span>
        </div>
      </div>

      <!-- Summary -->
      <div>
        <h3 class="text-xs font-bold text-slate-400 uppercase tracking-widest mb-1">Summary</h3>
        <table class="w-full text-xs text-slate-700">
          <thead><tr class="text-slate-400 text-left"><th class="font-bold pb-1"></th><th class="font-bold pb-1">Spin free</th><th class="font-bold pb-1">SOC</th></tr></thead>
          <tbody>
            {#each [["DFT gap", "dft_gap_ev"], ["QP gap", "qp_gap_ev"], ["Lowest exciton", "lowest_ev", "lowest_f"],
                    ["First bright exciton", "first_bright_ev", "first_bright_f"]] as [label, eKey, fKey]}
              <tr class="border-t border-slate-100">
                <td class="py-1.5 font-bold text-slate-900">{label}</td>
                {#each ["spin_free", "soc"] as s}
                  {@const v = data.summary?.[s] || {}}
                  <td class="py-1.5 font-mono whitespace-nowrap">{fmt(v[eKey])} eV{#if fKey && v[fKey] !== undefined}<div class="text-slate-400">f = {fmt(v[fKey])}</div>{/if}</td>
                {/each}
              </tr>
            {/each}
          </tbody>
        </table>
        <p class="text-[11px] text-slate-500 mt-3 leading-relaxed">
          QP: gap between the quasiparticle-corrected HOMO and LUMO. Exciton energies include
          the electron–hole interaction. Bright: f ≥ {brightF}.
          {#if data.summary?.qdex_commit}<br />QDEX <span class="font-mono">{data.summary.qdex_commit}</span>{/if}
        </p>
      </div>
    </div>

    <!-- Band-edge excitons -->
    <div>
      <h3 class="text-xs font-bold text-slate-400 uppercase tracking-widest mb-1">
        Band-edge excitons · {KINDS[kind].label}
        <span class="normal-case font-medium text-slate-400">· {excitons.length} states within {data.edge_window_ev ?? "?"} eV of the lowest</span>
      </h3>
      <div class="max-h-[340px] overflow-y-auto border border-slate-100 rounded-xl">
        <table class="w-full text-xs text-slate-700">
          <thead class="sticky top-0 bg-slate-50">
            <tr class="text-slate-500 text-left">
              <th class="font-bold px-3 py-1.5">State</th><th class="font-bold px-3 py-1.5">E (eV)</th><th class="font-bold px-3 py-1.5">f</th>
              <th class="font-bold px-3 py-1.5">Transition</th><th class="font-bold px-3 py-1.5">Type</th>
              <th class="font-bold px-3 py-1.5" title="Mean electron–hole distance">d<sub>eh</sub> (Å)</th>
              <th class="font-bold px-3 py-1.5" title="Charge-transfer character (0: local, 1: full CT)">CT character</th>
            </tr>
          </thead>
          <tbody>
            {#each excitons as x}
              <tr class="border-t border-slate-100 {x.f_osc >= brightF ? 'font-bold text-slate-900' : ''}">
                <td class="px-3 py-1 font-mono">{x.state}</td>
                <td class="px-3 py-1 font-mono">{fmt(x.energy_ev)}</td>
                <td class="px-3 py-1 font-mono">{fmt(x.f_osc, 4)}</td>
                <td class="px-3 py-1 font-mono">{x.transition}</td>
                <td class="px-3 py-1">{x.type}</td>
                <td class="px-3 py-1 font-mono">{fmt(x.d_eh_ang, 2)}</td>
                <td class="px-3 py-1 font-mono">{fmt(x.ct_character, 2)}</td>
              </tr>
            {/each}
          </tbody>
        </table>
      </div>
    </div>
  </div>
{/if}
