<script>
  // Small SVG chart for the ligand-shell tab: lines, bands, bars, reference lines and a hover readout.
  //   layers: [{type: "line", pts: [[x, y]], color, dash, width}
  //            {type: "band", lo: [[x, y]], hi: [[x, y]], color, opacity}
  //            {type: "bars", bars: [{x, y0, y1, w, color, opacity}]}      (w in x units)
  //            {type: "ticks", ticks: [{x, y, w, color}]}                     (short horizontal marks)
  //            {type: "points", pts: [[x, y]], color, r, hollow}
  //            {type: "hline", y, color, dash, label} / {type: "vline", x, color, dash, label}
  //            {type: "cells", cells: [{x0, x1, y0, y1, color}]}]
  //   hover(x) -> {title, rows: [{color, dash, text}], x?} or null
  let { layers = [], xDomain, yDomain, xLabel = "", yLabel = "", xTicks = null, yTicks = null,
        height = 260, width = 720, hover = null, legend = [], ariaLabel = "chart" } = $props();

  // viewBox width = the rendered width (text stays at its pixel size); `width` until it is measured
  let cw = $state(0);
  let W = $derived(cw > 0 ? Math.max(420, cw) : width);
  const M = { l: 58, r: 14, t: 14, b: 40 };
  const clipId = "clip-" + Math.random().toString(36).slice(2);
  let H = $derived(height);
  let PW = $derived(W - M.l - M.r), PH = $derived(H - M.t - M.b);
  const sx = (x) => M.l + ((x - xDomain[0]) / (xDomain[1] - xDomain[0])) * PW;
  const sy = (y) => M.t + PH - ((Math.min(Math.max(y, yDomain[0] - 1e6), yDomain[1] + 1e6) - yDomain[0]) / (yDomain[1] - yDomain[0])) * PH;
  const path = (pts) => pts.filter(([x, y]) => isFinite(x) && isFinite(y))
    .map(([x, y], i) => `${i ? "L" : "M"}${sx(x).toFixed(1)},${sy(y).toFixed(1)}`).join("");
  const bandPath = (lo, hi) => path(lo) + hi.slice().reverse()
    .filter(([x, y]) => isFinite(x) && isFinite(y)).map(([x, y]) => `L${sx(x).toFixed(1)},${sy(y).toFixed(1)}`).join("") + "Z";

  function nice(a, b, n = 5) {
    const raw = (b - a) / n, mag = 10 ** Math.floor(Math.log10(raw));
    const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => s >= raw);
    const out = [];
    for (let v = Math.ceil(a / step - 1e-9) * step; v <= b + 1e-9 * step; v += step) out.push(+v.toPrecision(10));
    return out.map((v) => ({ v, label: String(+v.toFixed(6)) }));
  }
  let xt = $derived(xTicks || nice(xDomain[0], xDomain[1], 8));
  let yt = $derived(yTicks || nice(yDomain[0], yDomain[1], 5));

  let hx = $state(null);
  function onMove(ev) {
    if (!hover) return;
    const r = ev.currentTarget.getBoundingClientRect();
    const x = ((ev.clientX - r.left) / r.width) * W;
    hx = x < M.l || x > W - M.r ? null : xDomain[0] + ((x - M.l) / PW) * (xDomain[1] - xDomain[0]);
  }
  let info = $derived(hover && hx !== null ? hover(hx) : null);
  let infoX = $derived(info ? (info.x ?? hx) : null);
</script>

<div class="relative" bind:clientWidth={cw}>
  <svg viewBox="0 0 {W} {H}" class="w-full h-auto select-none" role="img" aria-label={ariaLabel}
       onpointermove={onMove} onpointerleave={() => (hx = null)}>
    <defs><clipPath id={clipId}><rect x={M.l} y={M.t} width={PW} height={PH} /></clipPath></defs>
    {#each yt as t}
      <line x1={M.l} x2={W - M.r} y1={sy(t.v)} y2={sy(t.v)} stroke="#eef2f6" stroke-width="1" />
      <text x={M.l - 6} y={sy(t.v) + 3} text-anchor="end" font-size="11" fill="#94a3b8">{t.label}</text>
    {/each}
    {#each xt as t}
      <line x1={sx(t.v)} x2={sx(t.v)} y1={M.t + PH} y2={M.t + PH + 4} stroke="#94a3b8" />
      <text x={sx(t.v)} y={M.t + PH + 15} text-anchor="middle" font-size="11" fill="#64748b">{t.label}</text>
    {/each}
    <g clip-path="url(#{clipId})">
      {#each layers as ly}
        {#if ly.type === "cells"}
          {#each ly.cells as c}
            <rect x={sx(c.x0)} y={sy(c.y1)} width={Math.max(0.5, sx(c.x1) - sx(c.x0) + 0.4)} height={Math.max(0.5, sy(c.y0) - sy(c.y1) + 0.4)} fill={c.color} />
          {/each}
        {:else if ly.type === "band"}
          <path d={bandPath(ly.lo, ly.hi)} fill={ly.color} opacity={ly.opacity ?? 0.15} />
        {:else if ly.type === "bars"}
          {#each ly.bars as b}
            {@const ya = sy(Math.max(b.y0, b.y1))}
            {@const yb = sy(Math.min(b.y0, b.y1))}
            <rect x={sx(b.x - b.w / 2)} y={ya} width={Math.max(1, sx(b.x + b.w / 2) - sx(b.x - b.w / 2))} height={Math.max(1, yb - ya)}
                  fill={b.color} opacity={b.opacity ?? 1} rx="2" />
          {/each}
        {:else if ly.type === "ticks"}
          {#each ly.ticks as k}
            <line x1={sx(k.x - k.w / 2)} x2={sx(k.x + k.w / 2)} y1={sy(k.y)} y2={sy(k.y)} stroke={k.color} stroke-width="3" stroke-linecap="round" />
          {/each}
        {:else if ly.type === "line"}
          <path d={path(ly.pts)} fill="none" stroke={ly.color} stroke-width={ly.width ?? 2} stroke-dasharray={ly.dash || ""}
                stroke-linejoin="round" stroke-linecap="round" />
        {:else if ly.type === "points"}
          {#each ly.pts as [x, y]}
            {#if isFinite(y)}
              <circle cx={sx(x)} cy={sy(y)} r={ly.r ?? 3.5} fill={ly.hollow ? "white" : ly.color} stroke={ly.color} stroke-width="1.5" />
            {/if}
          {/each}
        {/if}
      {/each}
    </g>
    {#each layers as ly}
      {#if ly.type === "hline" && ly.y >= yDomain[0] && ly.y <= yDomain[1]}
        <line x1={M.l} x2={W - M.r} y1={sy(ly.y)} y2={sy(ly.y)} stroke={ly.color || "#334155"} stroke-width="1.2" stroke-dasharray={ly.dash ?? "5 4"} />
        {#if ly.label}<text x={W - M.r - 4} y={sy(ly.y) - 5} text-anchor="end" font-size="11" fill={ly.color || "#334155"}>{ly.label}</text>{/if}
      {:else if ly.type === "vline" && ly.x >= xDomain[0] && ly.x <= xDomain[1]}
        <line x1={sx(ly.x)} x2={sx(ly.x)} y1={M.t} y2={M.t + PH} stroke={ly.color || "#94a3b8"} stroke-width="1" stroke-dasharray={ly.dash ?? "2 3"} />
        {#if ly.label}<text x={sx(ly.x) + 4} y={M.t + 10} font-size="11" fill={ly.color || "#64748b"}>{ly.label}</text>{/if}
      {/if}
    {/each}
    <line x1={M.l} x2={W - M.r} y1={M.t + PH} y2={M.t + PH} stroke="#94a3b8" />
    <text x={M.l + PW / 2} y={H - 4} text-anchor="middle" font-size="12" fill="#475569">{xLabel}</text>
    <text transform="translate(13 {M.t + PH / 2}) rotate(-90)" text-anchor="middle" font-size="12" fill="#475569">{yLabel}</text>
    {#if info && infoX !== null}
      <line x1={sx(infoX)} x2={sx(infoX)} y1={M.t} y2={M.t + PH} stroke="#cbd5e1" stroke-width="1" />
    {/if}
  </svg>
  {#if info}
    <div class="absolute pointer-events-none bg-white/95 border border-slate-200 rounded-lg shadow-sm px-2.5 py-1.5 text-[11px] text-slate-700 z-10"
         style="left: {(sx(infoX) / W) * 100}%; top: 0.5rem; transform: translateX({sx(infoX) > W * 0.62 ? 'calc(-100% - 8px)' : '8px'});">
      <div class="font-bold text-slate-900">{@html info.title}</div>
      {#each info.rows as r}
        <div class="flex items-center gap-1.5 whitespace-nowrap font-mono">
          {#if r.color}<svg width="14" height="6"><line x1="0" x2="14" y1="3" y2="3" stroke={r.color} stroke-width="2.5" stroke-dasharray={r.dash || ""} /></svg>{/if}
          {@html r.text}
        </div>
      {/each}
    </div>
  {/if}
  {#if legend.length}
    <div class="flex flex-wrap gap-x-4 gap-y-1 text-[11px] text-slate-500 mt-1 px-1">
      {#each legend as g}
        <span class="flex items-center gap-1.5">
          <svg width="18" height="10">
            {#if g.kind === "box"}<rect x="4" y="0" width="10" height="10" rx="2" fill={g.color} opacity={g.opacity ?? 1} />
            {:else if g.kind === "circle"}<circle cx="9" cy="5" r="3.5" fill="white" stroke={g.color} stroke-width="1.5" />
            {:else if g.kind === "tick"}<line x1="3" x2="15" y1="5" y2="5" stroke={g.color} stroke-width="3" stroke-linecap="round" />
            {:else}<line x1="0" x2="18" y1="5" y2="5" stroke={g.color} stroke-width="2" stroke-dasharray={g.dash || ""} />{/if}
          </svg>{@html g.label}
        </span>
      {/each}
    </div>
  {/if}
</div>
