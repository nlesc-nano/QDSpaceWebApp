<script>
  // Excited states tab: QDEX excitons of the published DFT chain, spin free and with spin-orbit coupling,
  // with the surface-trap filter off or on. Data: dft/qdex/webapp.json (and webapp_traps_filtered.json):
  // the summary, the band-edge excitons with their descriptors, and the absorption spectrum of every state.
  import Plot from "./Plot.svelte";
  import { rebroaden } from "./broaden.js";

  let { url, dashboards = null, functional = "PBE" } = $props();   // url: .../webapp.json; dashboards: dftPages group

  let runs = $state({});                   // {included, filtered}
  let status = $state("loading");
  $effect(() => {
    const u = url;
    status = "loading";
    const get = (v) => fetch(v).then((r) => (r.ok && (r.headers.get("content-type") || "").includes("json") ? r.json() : null)).catch(() => null);
    Promise.all([get(u), get(u.replace(/webapp\.json$/, "webapp_traps_filtered.json"))]).then(([a, b]) => {
      if (u !== url) return;
      runs = { included: a, filtered: b };
      status = a ? "ready" : "missing";
    });
  });

  // ---------------- choices ----------------
  let spin = $state("sf");                 // "sf" | "soc"
  let traps = $state("included");          // "included" | "filtered"
  let sigma = $state(0.03);                // Gaussian broadening of the spectrum (eV)
  let unit = $state("ev");                 // energy axis: "ev" | "nm"
  let range = $state("edge");              // spectrum energy range: "edge" | "all"
  let sel = $state(null);                  // highlighted state number
  let run = $derived(runs[traps] || runs.included);
  let hasSoc = $derived(!!runs.included?.excitons?.soc?.length);
  let states = $derived(run?.excitons?.[spin] || []);
  let sum = $derived(run?.summary?.[spin === "sf" ? "spin_free" : "soc"] || {});
  let sumOther = $derived(runs[traps === "filtered" ? "included" : "filtered"]?.summary?.[spin === "sf" ? "spin_free" : "soc"] || null);
  // the same band edge with and without the filter (energies and strengths; degenerate pairs may swap labels)
  const edgeKey = (list) => (list || []).map((s) => `${s.energy_ev}:${s.f_osc}`).sort().join("|");
  let sameAsUnfiltered = $derived(traps === "filtered" && !!runs.included && edgeKey(runs.included.excitons?.[spin]) === edgeKey(states));
  $effect(() => { if (!hasSoc && spin === "soc") spin = "sf"; });
  $effect(() => { spin; traps; sel = null; });

  const seg = (on) => `px-3 py-1.5 rounded-lg text-xs font-bold transition-all ${on ? "bg-white text-slate-900 shadow-sm" : "text-slate-500 hover:text-slate-800"}`;
  const fmt = (v, n = 2) => (v === null || v === undefined || !isFinite(v) ? "–" : Number(v).toFixed(n));
  const nm = (e) => (e ? 1239.84193 / e : null);
  const arrow = (t) => String(t || "").replace("->", " → ").replace(/sp(H|L)/g, "$1");

  // ---------------- energy axes ----------------
  // descriptors exist for the band-edge states; every panel shares their window so a state lines up across panels
  let eLo = $derived(states.length ? states[0].energy_ev : 0);
  let eHi = $derived(states.length ? states.at(-1).energy_ev : 1);
  let edgeDomain = $derived([Math.floor((eLo - 0.25) * 20) / 20, Math.ceil((eHi + 0.2) * 20) / 20]);
  let specDomain = $derived(range === "edge" ? edgeDomain : [edgeDomain[0], Math.ceil(((sum.emax_ev || eHi) + 0.3) * 10) / 10]);

  // ---------------- spectrum ----------------
  let fmax = $derived(Math.max(1e-9, ...states.map((s) => s.f_osc)));
  let spec = $derived.by(() => {
    const sp = run?.spectrum?.[spin];
    if (!sp) return null;
    // S(E) = Σ f_i g_σ(E − E_i) with unit-area Gaussians; times σ√(2π) an isolated state's peak is as tall as its f.
    // σ >= 0.03 eV: the stored 0.03 eV curve re-broadened (exact, and it includes the states QDEX solved above emax);
    // narrower: summed from the sticks (every stored state), complete up to emax − 3σ
    const [e0, , de] = sp.energy_ev, c = sigma * Math.sqrt(2 * Math.PI);
    let pts, cut = null;
    if (sigma >= 0.03 - 1e-9 || !sp.sticks) {
      pts = rebroaden(e0, de, sp["sigma_0.030"] || [], 0.03, Math.max(sigma, 0.03)).map(([x, v]) => [x, v * c]);
    } else {
      const E = sp.sticks.energy_ev, F = sp.sticks.f_osc, step = Math.min(0.002, sigma / 4);
      cut = (sum.emax_ev ?? Infinity) - 3 * sigma;
      const hi = Math.min(specDomain[1], cut);
      pts = [];
      for (let x = specDomain[0]; x <= hi + 1e-9; x += step) {
        let v = 0;
        for (let i = 0; i < E.length; i++) { const d = (x - E[i]) / sigma; if (d > -6 && d < 6) v += F[i] * Math.exp(-0.5 * d * d); }
        pts.push([x, v]);
      }
    }
    pts = pts.filter(([x]) => x >= specDomain[0] && x <= specDomain[1]);
    const inWin = states.filter((s) => s.energy_ev >= specDomain[0] && s.energy_ev <= specDomain[1]).map((s) => s.f_osc);
    // sticks of the states above the band edge (no descriptors), in grey
    const eTop = states.at(-1)?.energy_ev ?? 0;
    const others = sp.sticks ? sp.sticks.energy_ev.map((e, i) => [e, sp.sticks.f_osc[i]]).filter(([e]) => e > eTop + 1e-6 && e >= specDomain[0] && e <= specDomain[1]) : [];
    return { pts, cut: cut !== null && cut < specDomain[1] ? cut : null, others, ymax: Math.max(1e-9, ...pts.map((p) => p[1]), ...inWin, ...others.map((o) => o[1])) * 1.08 };
  });
  let specLayers = $derived.by(() => {
    if (!spec) return [];
    const w = (specDomain[1] - specDomain[0]) / 450;
    return [
      { type: "band", lo: spec.pts.map(([x]) => [x, 0]), hi: spec.pts, color: "#7c3aed", opacity: 0.12 },
      { type: "line", pts: spec.pts, color: "#7c3aed", width: 2 },
      { type: "bars", bars: spec.others.map(([e, f]) => ({ x: e, y0: 0, y1: f, w, color: "#94a3b8", opacity: 0.7 })) },
      ...(spec.cut ? [{ type: "vline", x: spec.cut, color: "#94a3b8", dash: "3 3", label: unit === "nm" ? "complete above" : "complete below" }] : []),
      { type: "bars", bars: states.map((s) => ({ x: s.energy_ev, y0: 0, y1: s.f_osc, w, color: TYPE[s.type]?.color || "#334155", opacity: sel === s.state ? 1 : 0.8 })) },
      ...(sum.qp_gap_ev ? [{ type: "vline", x: sum.qp_gap_ev, color: "#0f172a", dash: "5 4", label: "QP gap" }] : []),
      ...selLine,
    ];
  });

  // ---------------- descriptors ----------------
  const TYPE = { Wannier: { color: "#2563eb", label: "Wannier (delocalised)" }, CT: { color: "#ea580c", label: "charge transfer" },
                 Frenkel: { color: "#16a34a", label: "Frenkel (local)" } };
  let types = $derived([...new Set(states.map((s) => s.type))]);
  const rad = (s) => 2.5 + 5.5 * Math.sqrt(Math.max(0, s.f_osc) / fmax);
  // one points layer per exciton type; key(s) gives y
  const byType = (key, extra = {}) => types.map((t) => ({ type: "points", color: TYPE[t]?.color || "#334155", opacity: 0.55, ...extra,
    pts: states.filter((s) => s.type === t && isFinite(key(s))).map((s) => [s.energy_ev, key(s), rad(s)]) }));
  let selState = $derived(states.find((s) => s.state === sel) || null);
  let selLine = $derived(selState ? [{ type: "vline", x: selState.energy_ev, color: "#db2777", dash: "", label: `#${selState.state}` }] : []);
  const yMax = (vals, pad = 1.12) => Math.max(1e-9, ...vals.filter(isFinite)) * pad;

  const nearest = (x) => (states.length ? states.reduce((b, s) => (Math.abs(s.energy_ev - x) < Math.abs(b.energy_ev - x) ? s : b), states[0]) : null);
  const stateRows = (s) => [
    { text: `f = ${fmt(s.f_osc, 3)} · ${arrow(s.transition)} · ${s.type}` },
    { text: `d<sub>eh</sub> ${fmt(s.d_eh_ang, 1)} Å · d<sub>CT</sub> ${fmt(s.d_ct_ang, 1)} Å · CT ${fmt(s.ct_character, 2)}` },
    { text: `σ<sub>h</sub> ${fmt(s.sigma_h_ang, 1)} Å · σ<sub>e</sub> ${fmt(s.sigma_e_ang, 1)} Å` },
    ...(s.singlet_fraction !== undefined ? [{ text: `singlet ${fmt(100 * s.singlet_fraction, 0)}%` }] : []),
  ];
  const stateHover = (x) => {
    const s = nearest(x);
    return s ? { x: s.energy_ev, title: `#${s.state} · ${fmt(s.energy_ev, 3)} eV · ${fmt(nm(s.energy_ev), 0)} nm`, rows: stateRows(s) } : null;
  };
  const specHover = (x) => {
    const h = stateHover(x);
    if (!spec) return h;
    const p = spec.pts.reduce((b, q) => (Math.abs(q[0] - x) < Math.abs(b[0] - x) ? q : b), spec.pts[0]);
    const near = h && Math.abs(h.x - x) < (specDomain[1] - specDomain[0]) / 60;
    return near ? h : { x: p[0], title: `${fmt(p[0], 3)} eV · ${fmt(nm(p[0]), 0)} nm`, rows: [{ text: `broadened f ${fmt(p[1], 3)}` }] };
  };

  let panels = $derived.by(() => {
    if (!states.length) return [];
    const out = [
      { title: "Exciton size and charge-transfer distance", caption: "d<sub>eh</sub>: root-mean-square electron–hole distance; d<sub>CT</sub>: distance between the hole and electron centroids (open circles).",
        layers: [...byType((s) => s.d_eh_ang), ...byType((s) => s.d_ct_ang, { hollow: true, opacity: 1 })], y: [0, yMax(states.flatMap((s) => [s.d_eh_ang, s.d_ct_ang]))], yLabel: "distance (Å)",
        legend: [{ kind: "dot", color: "#475569", label: "d<sub>eh</sub>" }, { kind: "circle", color: "#475569", label: "d<sub>CT</sub>" }] },
      { title: "Hole and electron spread", caption: "σ<sub>h</sub>, σ<sub>e</sub>: the radius of the hole and electron clouds; against the dot's radius, how far each carrier spreads.",
        layers: [{ type: "points", color: "#dc2626", opacity: 0.55, pts: states.map((s) => [s.energy_ev, s.sigma_h_ang, rad(s)]) },
                 { type: "points", color: "#0284c7", opacity: 0.55, pts: states.map((s) => [s.energy_ev, s.sigma_e_ang, rad(s)]) }],
        y: [0, yMax(states.flatMap((s) => [s.sigma_h_ang, s.sigma_e_ang]))], yLabel: "spread (Å)",
        legend: [{ kind: "dot", color: "#dc2626", label: "hole σ<sub>h</sub>" }, { kind: "dot", color: "#0284c7", label: "electron σ<sub>e</sub>" }] },
      { title: "Charge-transfer character", caption: "d<sub>CT</sub> / d<sub>eh</sub>: 0 when the hole and electron share a centre, 1 when they are fully separated. Above 0.6 the state counts as charge transfer.",
        layers: [...byType((s) => s.ct_character), { type: "hline", y: 0.6, color: "#ea580c", label: "CT above 0.6", labelAt: "start" }], y: [0, 1], yLabel: "CT character" },
      { title: "Electron–hole interaction", caption: "E = D<sub>QP</sub> + K<sub>x</sub> − K<sub>d</sub>: the quasiparticle transition energy, the exchange (repulsive, open circles) and the screened direct attraction (filled).",
        layers: [...byType((s) => s.minus_kd_ev), ...byType((s) => s.kx_ev, { hollow: true, opacity: 1 }), { type: "hline", y: 0, color: "#94a3b8", dash: "" }],
        y: [Math.min(0, ...states.map((s) => s.minus_kd_ev)) * 1.12, Math.max(0.05, ...states.map((s) => s.kx_ev)) * 1.3], yLabel: "energy (eV)",
        legend: [{ kind: "dot", color: "#475569", label: "−K<sub>d</sub> (direct)" }, { kind: "circle", color: "#475569", label: "K<sub>x</sub> (exchange)" }] },
    ];
    if (spin === "soc" && states.some((s) => s.singlet_fraction !== undefined))
      out.push({ title: "Singlet character", caption: "Weight of the spin-singlet component of each spin-orbit state; mostly triplet states are dark.",
                 layers: byType((s) => s.singlet_fraction), y: [0, 1], yLabel: "singlet fraction" });
    return out.map((p) => ({ ...p, layers: [...p.layers, ...selLine] }));
  });

  // ---------------- energy or wavelength axis ----------------
  // everything is built in eV; in nm the x values are mapped at the plot (hc / E), the broadening stays in energy
  const HC = 1239.84193;
  const X = (e) => (unit === "nm" ? HC / e : e);
  const dom = ([a, b]) => (unit === "nm" ? [Math.floor(HC / b / 5) * 5, Math.ceil(HC / a / 5) * 5] : [a, b]);
  let xLabel = $derived(unit === "nm" ? "wavelength (nm)" : "excitation energy (eV)");
  const toX = (layers) => unit === "ev" ? layers : layers.map((l) => {
    if (l.type === "line" || l.type === "points") return { ...l, pts: l.pts.map(([x, ...r]) => [X(x), ...r]) };
    if (l.type === "band") return { ...l, lo: l.lo.map(([x, y]) => [X(x), y]), hi: l.hi.map(([x, y]) => [X(x), y]) };
    if (l.type === "bars") return { ...l, bars: l.bars.map((b) => ({ ...b, x: X(b.x), w: (b.w * HC) / (b.x * b.x) })) };
    if (l.type === "vline") return { ...l, x: X(l.x) };
    return l;
  });
  const hov = (f) => (x) => { const h = f(unit === "nm" ? HC / x : x); return h ? { ...h, x: X(h.x) } : null; };

  // ---------------- key numbers ----------------
  let cards = $derived.by(() => {
    if (!states.length) return [];
    const s0 = states[0], b = sum.first_bright_state ? states[sum.first_bright_state - 1] : null;
    return [
      ["Lowest exciton", `${fmt(s0.energy_ev)} eV`, `${fmt(nm(s0.energy_ev), 0)} nm · f = ${fmt(s0.f_osc)} · ${arrow(s0.transition)}`],
      ["First bright exciton", sum.first_bright_ev ? `${fmt(sum.first_bright_ev)} eV` : "–",
       sum.first_bright_ev ? `state #${sum.first_bright_state} · f = ${fmt(sum.first_bright_f)} (f ≥ ${run?.summary?.bright_f})` : `none with f ≥ ${run?.summary?.bright_f}`],
      ["Quasiparticle gap", `${fmt(sum.qp_gap_ev)} eV`, `DFT ${functional} gap ${fmt(sum.dft_gap_ev)} eV`],
      ["Exciton binding (lowest)", `${fmt(s0.d_qp_ev - s0.energy_ev)} eV`, "D<sub>QP</sub> − E of the lowest state"],
      ["States computed", `${sum.n_states_stored ?? "–"}`, `up to ${fmt(sum.emax_ev, 1)} eV; ${states.length} band-edge states analysed`],
      ["Lowest state's character", s0.type, `d<sub>eh</sub> ${fmt(s0.d_eh_ang, 1)} Å · CT ${fmt(s0.ct_character)}`],
    ];
  });
  let showAll = $state(false);
  let dashUrl = $derived(dashboards ? (traps === "filtered" ? dashboards[`${spin}_trap`] : null) || dashboards[spin] : null);
</script>

{#snippet section(title, subtitle)}
  <div class="flex items-baseline gap-3 flex-wrap border-b border-slate-200 pb-2 mt-2">
    <h3 class="text-base font-bold text-slate-900">{@html title}</h3>
    <span class="text-xs text-slate-500">{@html subtitle}</span>
  </div>
{/snippet}
{#snippet panel(title, caption)}
  <h4 class="text-xs font-bold text-slate-500 uppercase tracking-widest mb-1">{@html title}</h4>
  {#if caption}<p class="text-[11px] text-slate-500 mb-1 leading-relaxed">{@html caption}</p>{/if}
{/snippet}

{#if status === "loading"}
  <div class="flex items-center justify-center py-20 text-brand-600 text-sm font-bold">
    <div class="w-6 h-6 border-4 border-brand-500 border-t-transparent rounded-full animate-spin mr-3"></div>Loading…
  </div>
{:else if status === "missing"}
  <p class="text-sm text-slate-500 py-10 text-center">No QDEX excitons published for this dot.</p>
{:else}
  <div class="flex flex-col gap-6 pb-6">
    {@render section("Excited states", `QDEX on the DFT ${functional} ground state`)}
    <p class="text-sm text-slate-700 max-w-4xl leading-relaxed">
      Each exciton is one hole → electron transition between quasiparticle levels (DFT {functional} orbitals with a scaled-GW correction),
      with the screened electron–hole attraction and the exchange from atom-centred transition charges (diagonal BSE, Resta screening).
      Its descriptors say how big the exciton is, how far the electron sits from the hole and how far each carrier spreads.
      With surface traps filtered out, transitions from or to trap states (Facet or Localized by their angular coverage) are left out.
    </p>
    <div class="flex flex-wrap items-center gap-3 text-xs">
      <div class="flex gap-1 bg-slate-100 p-1 rounded-xl border border-slate-200">
        <button class={seg(spin === "sf")} onclick={() => (spin = "sf")}>spin free</button>
        <button class={seg(spin === "soc")} disabled={!hasSoc} onclick={() => (spin = "soc")}>spin–orbit</button>
      </div>
      {#if runs.filtered}
        <span class="font-bold text-slate-400">Surface traps</span>
        <div class="flex gap-1 bg-slate-100 p-1 rounded-xl border border-slate-200">
          <button class={seg(traps === "included")} onclick={() => (traps = "included")}>included</button>
          <button class={seg(traps === "filtered")} onclick={() => (traps = "filtered")}>filtered out</button>
        </div>
      {/if}
      {#if dashUrl}<a href={dashUrl} target="_blank" rel="noopener" class="ml-auto text-brand-700 font-bold hover:underline">QDEX dashboard ↗</a>{/if}
    </div>
    {@render cardsSnippet()}
    {#if traps === "filtered" && sameAsUnfiltered}
      <p class="text-xs text-amber-800 bg-amber-50 border border-amber-100 rounded-lg px-3 py-2 -mt-3">No band-edge {spin === "sf" ? "spin-free" : "spin–orbit"} state involves a trap level, so the band edge is the same as with traps included;
        the filter removes {(sumOther?.n_states_stored ?? 0) - (sum.n_states_stored ?? 0)} higher states ({sumOther?.n_states_stored} → {sum.n_states_stored} up to {fmt(sum.emax_ev, 1)} eV).</p>
    {:else if sumOther?.lowest_ev && Math.abs(sumOther.lowest_ev - sum.lowest_ev) > 1e-4}
      <p class="text-[11px] text-slate-500 -mt-3">With surface traps {traps === "filtered" ? "included" : "filtered out"} the lowest exciton is at {fmt(sumOther.lowest_ev)} eV (f = {fmt(sumOther.lowest_f)}).</p>
    {/if}

    <!-- spectrum (top left) and the states -->
    <div class="grid grid-cols-1 xl:grid-cols-2 gap-6">
      <div>
        <div class="flex flex-wrap items-center gap-3 mb-1">
          <h4 class="text-xs font-bold text-slate-500 uppercase tracking-widest">Absorption spectrum</h4>
          <div class="flex gap-1 bg-slate-100 p-1 rounded-xl border border-slate-200 text-xs">
            <button class={seg(range === "edge")} onclick={() => (range = "edge")}>band edge</button>
            <button class={seg(range === "all")} onclick={() => (range = "all")}>to {fmt(sum.emax_ev, 1)} eV</button>
          </div>
          <label class="flex items-center gap-2 text-xs">
            <span class="font-bold text-slate-500">σ</span>
            <input type="range" min="0.01" max="0.3" step="0.005" value={sigma} oninput={(e) => (sigma = +e.currentTarget.value)} class="w-32 accent-brand-600" />
            <span class="font-mono text-slate-900 w-16">{sigma.toFixed(3)} eV</span>
          </label>
          <div class="flex gap-1 bg-slate-100 p-1 rounded-xl border border-slate-200 text-xs">
            <button class={seg(unit === "ev")} onclick={() => (unit = "ev")}>eV</button>
            <button class={seg(unit === "nm")} onclick={() => (unit = "nm")}>nm</button>
          </div>
        </div>
        <p class="text-[11px] text-slate-500 mb-1 leading-relaxed">Sticks: oscillator strengths, coloured by character for the band-edge states, grey above. Line: every computed state broadened by a Gaussian of width σ (in energy), scaled so that an isolated state's peak is as tall as its f.{#if spec?.cut}{' '}Below σ = 0.03 eV the line is summed from the stored states, which end at {fmt(sum.emax_ev, 1)} eV: it is complete on the low-energy side of the grey dashed line.{/if} Hover for a state.</p>
        <Plot layers={toX(specLayers)} xDomain={dom(specDomain)} yDomain={[0, spec?.ymax || 1]} xLabel={xLabel} yLabel="oscillator strength f"
              height={320} hover={hov(specHover)} ariaLabel="absorption spectrum"
              legend={[{ color: "#7c3aed", label: `broadened, σ = ${sigma.toFixed(3)} eV` }, ...types.map((t) => ({ kind: "box", color: TYPE[t]?.color || "#334155", label: `f, ${TYPE[t]?.label || t}` }))]} />
      </div>
      <div>
        {@render panel("Band-edge states", `${states.length} states within ${run?.edge_window_ev} eV of the lowest. Click a row to mark it in every plot.`)}
        <div class="max-h-[330px] overflow-y-auto border border-slate-100 rounded-xl">
          <table class="w-full text-xs">
            <thead class="sticky top-0 bg-slate-50 text-slate-500">
              <tr class="text-right">
                <th class="px-2 py-1.5 text-left">#</th><th class="px-2">E (eV)</th><th class="px-2">λ (nm)</th><th class="px-2">f</th>
                <th class="px-2 text-left">transition</th><th class="px-2 text-left">type</th><th class="px-2">d<sub>eh</sub> (Å)</th><th class="px-2">CT</th>
                {#if spin === "soc"}<th class="px-2">singlet</th>{/if}
              </tr>
            </thead>
            <tbody>
              {#each showAll ? states : states.slice(0, 20) as s}
                <tr class="text-right font-mono border-t border-slate-100 cursor-pointer {sel === s.state ? 'bg-pink-50' : 'hover:bg-slate-50'}"
                    onclick={() => (sel = sel === s.state ? null : s.state)}>
                  <td class="px-2 py-1 text-left text-slate-400">{s.state}</td><td class="px-2">{fmt(s.energy_ev, 3)}</td><td class="px-2">{fmt(nm(s.energy_ev), 0)}</td>
                  <td class="px-2 {s.f_osc >= (run?.summary?.bright_f ?? 1) ? 'font-bold text-slate-900' : 'text-slate-500'}">{fmt(s.f_osc, 3)}</td>
                  <td class="px-2 text-left">{arrow(s.transition)}</td>
                  <td class="px-2 text-left font-sans"><span class="inline-block w-2 h-2 rounded-full mr-1" style="background:{TYPE[s.type]?.color || '#334155'}"></span>{s.type}</td>
                  <td class="px-2">{fmt(s.d_eh_ang, 1)}</td><td class="px-2">{fmt(s.ct_character, 2)}</td>
                  {#if spin === "soc"}<td class="px-2">{fmt(100 * (s.singlet_fraction ?? 0), 0)}%</td>{/if}
                </tr>
              {/each}
            </tbody>
          </table>
          {#if states.length > 20}
            <button class="w-full text-xs font-bold text-brand-700 py-1.5 hover:bg-slate-50 border-t border-slate-100" onclick={() => (showAll = !showAll)}>
              {showAll ? "show the first 20" : `show all ${states.length}`}</button>
          {/if}
        </div>
      </div>
    </div>

    {@render section("Exciton descriptors", "the band-edge states, on the same energy axis; marker size: oscillator strength")}
    <div class="flex flex-wrap gap-x-4 gap-y-1 text-[11px] text-slate-500 -mt-3">
      {#each types as t}<span class="flex items-center gap-1.5"><span class="w-2.5 h-2.5 rounded-full inline-block" style="background:{TYPE[t]?.color || '#334155'}"></span>{TYPE[t]?.label || t}</span>{/each}
      <span>Character: charge transfer when CT &gt; 0.6, Frenkel when d<sub>eh</sub> &lt; 3 Å and CT &lt; 0.2, Wannier otherwise.</span>
    </div>
    <div class="grid grid-cols-1 xl:grid-cols-2 gap-6">
      {#each panels as p}
        <div>
          {@render panel(p.title, p.caption)}
          <Plot layers={toX(p.layers)} xDomain={dom(edgeDomain)} yDomain={p.y} xLabel={xLabel} yLabel={p.yLabel} height={260}
                hover={hov(stateHover)} legend={p.legend || []} ariaLabel={p.title} />
        </div>
      {/each}
    </div>
  </div>
{/if}

{#snippet cardsSnippet()}
  <div class="grid grid-cols-2 md:grid-cols-3 xl:grid-cols-6 gap-3">
    {#each cards as [label, value, note]}
      <div class="rounded-xl border p-3 bg-violet-50/40 border-violet-100">
        <div class="text-[11px] text-slate-500">{@html label}</div>
        <div class="text-lg font-bold text-slate-900">{@html value}</div>
        {#if note}<div class="text-[11px] text-slate-500 leading-snug">{@html note}</div>{/if}
      </div>
    {/each}
  </div>
{/snippet}
