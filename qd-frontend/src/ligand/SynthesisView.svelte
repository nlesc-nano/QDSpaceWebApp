<script>
  // Synthesis and post-synthesis tab of a dot [MX]_m[ML_n]_k0, in vacuum, in an implicit solvent, or both.
  //   Formation: the dot as built, from the X precursor (converted to MX monomers) and ML_n, with ML_n
  //   complexed by oleylamine: the supersaturation it needs (S*) and the shell it keeps there.
  //   Post-synthesis: the dots exist; dilution, oleylamine, added ML_n, washing and heating change their shell.
  // Free energies: MACE-MH-1 (+ DFT COSMO solvation); oleylamine complexes: MACE-MH-1. Formulas in model.js.
  import Plot from "./Plot.svelte";
  import { OAM_COMPLEXES } from "./oamComplexes.js";
  import { KB, UNITS, buildModel, solventModel, populations, dmuOf, amineLnZ, amineSpeciation, solveClosed,
           hasFormation, log10Sstar } from "./model.js";

  let { detUrl = null, stabUrl = null, cubesUrl = null } = $props();

  let model = $state(null);
  let sm = $state(null);
  let status = $state("idle");
  $effect(() => {
    const d = detUrl, s = stabUrl, c = cubesUrl;
    if (!d) { model = null; status = "none"; return; }
    status = "loading";
    const get = (u) => (u ? fetch(u).then((r) => (r.ok ? r.json() : null)).catch(() => null) : Promise.resolve(null));
    Promise.all([fetch(d).then((r) => { if (!r.ok) throw new Error(r.status); return r.json(); }), get(s), get(c)])
      .then(([det, stab, cubes]) => {
        if (d !== detUrl) return;
        model = buildModel(det, stab);
        sm = solventModel(cubes);
        medium = sm ? "both" : "gas";
        status = "ready";
      })
      .catch(() => { if (d === detUrl) status = "error"; });
  });

  // ---------------- settings ----------------
  let unit = $state("eV");
  let confMode = $state("both");
  let medium = $state("gas");
  let eps = $state(2.38);
  let delta = $state(0);                  // solution correction per bound OAm (eV)
  const C_NEAT = 3.04;                    // neat oleylamine (M)
  const SOLVENTS = [["hexane", 1.88], ["toluene", 2.38], ["chloroform", 4.81], ["THF", 7.58]];
  const PRECURSOR = { Se: "TOP-Se or Se in oleylamine", S: "S in oleylamine or (TMS)₂S", Te: "TOP-Te", P: "P(NEt₂)₃ or P(TMS)₃",
                      As: "As(NMe₂)₃", Sb: "Sb(NMe₂)₃" };

  const iTof = (t) => Math.max(0, model.T.findIndex((v) => Math.abs(v - t) < 1e-6));
  const nearestT = (x) => model.T.reduce((b, t, j) => (Math.abs(t - x) < Math.abs(model.T[b] - x) ? j : b), 0);
  const uf = () => UNITS[unit][0];
  const ul = () => UNITS[unit][1];
  const fe = (v, d = UNITS[unit][2]) => (isFinite(v) ? (v * uf() >= 0 ? "+" : "−") + Math.abs(v * uf()).toFixed(d) : "–");
  const fk = (v) => (isFinite(v) ? v.toFixed(1) : "–");
  const pow10 = (e) => "10" + String(e).replace("-", "⁻").replace(/\d/g, (d) => "⁰¹²³⁴⁵⁶⁷⁸⁹"[d]);
  const html = (s) => String(s).replace(/(\d+)/g, "<sub>$1</sub>");
  const cTxt = (c) => {
    if (!(c > 0)) return "0";
    const e = Math.floor(Math.log10(c)), m = c / 10 ** e;
    return e >= -2 && e <= 0 ? `${+c.toPrecision(2)} M` : `${m.toFixed(1)}·10<sup>${e}</sup> M`;
  };
  const confs = () => (!model?.hasConf ? [false] : confMode === "both" ? [true, false] : [confMode === "ideal"]);
  const CONF = { true: { color: "#eb6834", label: "ideal S<sub>conf</sub>" }, false: { color: "#2a78d6", label: "S<sub>conf</sub> = 0" } };
  const solvName = () => SOLVENTS.find(([, e]) => Math.abs(e - eps) < 0.005)?.[0];
  let media = $derived.by(() => {
    const gas = { key: "gas", sv: null, dash: "6 4", label: "vacuum" };
    const sol = { key: "sol", sv: sm ? { sm, eps } : null, dash: "", label: `solution, ε = ${eps.toFixed(2)}${solvName() ? ` (${solvName()})` : ""}` };
    return !sm || medium === "gas" ? [gas] : medium === "sol" ? [sol] : [gas, sol];
  });
  const combos = () => media.flatMap((m) => confs().map((c) => ({ m, c })));
  const comboLabel = ({ m, c }) => `${m.label}, ${CONF[c].label}`;
  const legend = (extra = []) => [...confs().map((c) => ({ color: CONF[c].color, label: CONF[c].label })),
                                   ...media.map((m) => ({ color: "#334155", dash: m.dash, label: m.label })), ...extra];

  let complexes = $derived(model ? OAM_COMPLEXES.complexes[model.unit] || null : null);
  const lnZ = (iT, cOAm) => amineLnZ(complexes, iT, model.T[iT], cOAm, delta);

  // ---------------- labels ----------------
  let dotHtml = $derived(model ? `[${html(model.core || "MX")}]<sub>${model.m ?? "m"}</sub>[${html(model.unit)}]<sub>${model.k0 ?? "k"}</sub>` : "");
  let unitHtml = $derived(model ? html(model.unit) : "");
  let coreHtml = $derived(model ? html(model.core || "MX") : "");
  let X = $derived(model?.coreX || "X");
  let precursorHtml = $derived(`${X} precursor`);
  const tempTicks = [100, 200, 300, 400, 500, 600, 700, 800].map((v) => ({ v, label: String(v) }));
  function kDomOf(layers) {
    const ys = layers.filter((l) => l.type === "line").flatMap((l) => l.pts.map((p) => p[1])).filter(isFinite).concat([model.k0 - 2, model.k0 + 2]);
    return [Math.max(model.states[0].k, Math.floor(Math.min(...ys))) - 0.5, Math.min(model.states.at(-1).k, Math.ceil(Math.max(...ys))) + 0.5];
  }
  const lineHover = (layers, fmtX) => (x) => {
    const lines = layers.filter((l) => l.type === "line" && l.name);
    if (!lines.length) return null;
    const rows = lines.map((l) => {
      const p = l.pts.reduce((b, q) => (Math.abs(q[0] - x) < Math.abs(b[0] - x) ? q : b), l.pts[0]);
      return { color: l.color, dash: l.dash, text: `${l.name}: ⟨k⟩ = ${fk(p[1])}`, x: p[0] };
    });
    return { x: rows[0].x, title: fmtX(rows[0].x), rows };
  };

  // =====================================================================================
  // Formation
  // =====================================================================================
  let fT = $state(500);
  let fLogCL = $state(-1);               // total ML_n
  let fLogCX = $state(-1.3);             // total X precursor
  let fXi = $state(0.5);                 // fraction of the X precursor converted
  let fOAm = $state(1.0);                // oleylamine (M), the rest a co-solvent
  const C_DOT_F = 1e-6;
  const F_DOM = [275, 800];
  const fTicks = [300, 400, 500, 600, 700, 800].map((v) => ({ v, label: String(v) }));
  let fI = $derived(model ? iTof(fT) : 0);
  let fLeft = $derived(10 ** fLogCL - fXi * 10 ** fLogCX);           // ML_n not used by the conversion (M)
  const fFree = (i, cOAm = fOAm) => Math.max(fLeft, 1e-30) / Math.exp(lnZ(i, cOAm));
  let fRows = $derived(model && hasFormation(model) && fLeft > 0 ? combos().map((cb) => ({ ...cb,
    S: log10Sstar(model, fI, fFree(fI), C_DOT_F, cb.c, cb.m.sv),
    k: populations(model, fI, dmuOf(fFree(fI), model.T[fI]), cb.c, cb.m.sv).kmean })) : []);
  let sStarLayers = $derived.by(() => {
    if (!model || !hasFormation(model) || fLeft <= 0) return [];
    const L = [];
    for (const [cOAm, op] of [[fOAm, 1], ...(fOAm > 0 && complexes ? [[0, 0.35]] : [])])
      for (const { m, c } of combos())
        L.push({ type: "line", color: CONF[c].color, dash: m.dash, opacity: op, name: op === 1 ? comboLabel({ m, c }) : null,
                 pts: model.T.map((t, i) => [t, t >= F_DOM[0] ? log10Sstar(model, i, fFree(i, cOAm), C_DOT_F, c, m.sv) : NaN]) });
    L.push({ type: "vline", x: model.T[fI], label: `${model.T[fI]} K` });
    return L;
  });
  let sStarDom = $derived.by(() => {
    const ys = sStarLayers.filter((l) => l.pts).flatMap((l) => l.pts.map((p) => p[1])).filter(isFinite);
    return ys.length ? [Math.floor(Math.min(...ys) - 0.2), Math.ceil(Math.max(...ys) + 0.2)] : [0, 1];
  });
  const sStarHover = (x) => {
    const i = nearestT(x);
    return { x: model.T[i], title: `T = ${model.T[i]} K · free ${unitHtml} ${cTxt(fFree(i))}`,
             rows: combos().map((cb) => ({ color: CONF[cb.c].color, dash: cb.m.dash, text: `${comboLabel(cb)}: log<sub>10</sub> S* = ${log10Sstar(model, i, fFree(i), C_DOT_F, cb.c, cb.m.sv).toFixed(2)}` })) };
  };
  let fShellLayers = $derived.by(() => {
    if (!model || fLeft <= 0) return [];
    const L = [];
    for (const [cOAm, op] of [[fOAm, 1], ...(fOAm > 0 && complexes ? [[0, 0.35]] : [])])
      for (const { m, c } of combos())
        L.push({ type: "line", color: CONF[c].color, dash: m.dash, opacity: op, name: op === 1 ? comboLabel({ m, c }) : null,
                 pts: model.T.map((t, i) => [t, t >= F_DOM[0] ? populations(model, i, dmuOf(fFree(i, cOAm), t), c, m.sv).kmean : NaN]) });
    L.push({ type: "hline", y: model.k0, color: "#334155", label: `as built, k = ${model.k0}` });
    L.push({ type: "vline", x: model.T[fI], label: `${model.T[fI]} K` });
    return L;
  });
  let speciesF = $derived(model && complexes ? amineSpeciation(complexes, fI, model.T[fI], fOAm, delta) : null);

  // =====================================================================================
  // Post-synthesis: closed system
  // =====================================================================================
  let pT = $state(300);
  let pLogDot = $state(-6);
  let pOAm = $state(0);
  let pLogAdd = $state(-12);              // -12 = nothing added
  let pWash = $state(0.9);
  let pI = $derived(model ? iTof(pT) : 0);
  let pAdd = $derived(pLogAdd <= -11.99 ? 0 : 10 ** pLogAdd);
  const closed = (i, cdot, cOAm, c, sv, add = pAdd, kRef = null) => solveClosed(model, i, cdot, add, lnZ(i, cOAm), c, kRef, sv);
  let pRows = $derived(model ? combos().map((cb) => ({ ...cb, ...closed(pI, 10 ** pLogDot, pOAm, cb.c, cb.m.sv) })) : []);
  const curves = (xs, f) => combos().map(({ m, c }) => ({ type: "line", color: CONF[c].color, dash: m.dash, name: comboLabel({ m, c }),
                                                         pts: xs.map((x) => [x, f(x, c, m.sv)]) }));
  const asBuilt = () => ({ type: "hline", y: model.k0, color: "#334155", label: `as built, k = ${model.k0}` });
  const LOGDOT = Array.from({ length: 25 }, (_, i) => -9 + 0.25 * i);
  let dilLayers = $derived(model ? [...curves(LOGDOT, (l, c, sv) => closed(pI, 10 ** l, pOAm, c, sv).kmean), asBuilt(), { type: "vline", x: pLogDot }] : []);
  const OAMS = Array.from({ length: 31 }, (_, i) => (C_NEAT * i) / 30);
  let oamLayers = $derived(model && complexes ? [...curves(OAMS, (a, c, sv) => closed(pI, 10 ** pLogDot, a, c, sv).kmean), asBuilt(), { type: "vline", x: pOAm }] : []);
  const ADDS = Array.from({ length: 33 }, (_, i) => -9 + 0.25 * i);
  let addLayers = $derived(model?.hasUp ? [...curves(ADDS, (l, c, sv) => closed(pI, 10 ** pLogDot, pOAm, c, sv, 10 ** l).kmean), asBuilt(),
                                          ...(pAdd ? [{ type: "vline", x: pLogAdd }] : [])] : []);
  let washLayers = $derived.by(() => {
    if (!model) return [];
    const cdot = 10 ** pLogDot, L = [];
    for (const { m, c } of combos()) {
      let r = closed(pI, cdot, pOAm, c, m.sv), pts = [[0, r.kmean]];
      for (let n = 1; n <= 8; n++) {
        const left = (1 - pWash) * Math.exp(lnZ(pI, pOAm)) * r.cFree;
        r = closed(pI, cdot, pOAm, c, m.sv, left, r.kmean);
        pts.push([n, r.kmean]);
      }
      L.push({ type: "line", color: CONF[c].color, dash: m.dash, name: comboLabel({ m, c }),
               pts: pts.flatMap(([n, k], j) => (j ? [[n - 0.5, pts[j - 1][1]], [n - 0.5, k], [n, k]] : [[0, k]])) });
      L.push({ type: "points", color: CONF[c].color, pts, r: 3 });
    }
    L.push(asBuilt());
    return L;
  });
  let annealLayers = $derived(model ? [...combos().map(({ m, c }) => ({ type: "line", color: CONF[c].color, dash: m.dash, name: comboLabel({ m, c }),
                                          pts: model.T.map((t, i) => [t, closed(i, 10 ** pLogDot, pOAm, c, m.sv).kmean]) })),
                                       asBuilt(), { type: "vline", x: model.T[pI], label: `${model.T[pI]} K` }] : []);
  let firstUnit = $derived.by(() => {
    if (!model) return null;
    const u = model.units.find((q) => q.k === model.k0);
    if (!u) return null;
    const kT = KB * model.T[pI];
    return { dG: u.cost[pI], dGc: u.cost[pI] - kT * u.dS, thr: kT * Math.log(1 / 10 ** pLogDot) + kT * lnZ(pI, pOAm) };
  });
  let speciesP = $derived(model && complexes ? amineSpeciation(complexes, pI, model.T[pI], pOAm, delta) : null);
  const seg = (on) => `px-3 py-1.5 rounded-lg font-bold transition-all ${on ? "bg-white text-slate-900 shadow-sm" : "text-slate-500 hover:text-slate-800"}`;
</script>

{#snippet slider(label, value, min, max, step, set, text)}
  <label class="flex items-center gap-2 text-xs text-slate-600">
    <span class="font-bold text-slate-500 whitespace-nowrap">{@html label}</span>
    <input type="range" {min} {max} {step} value={value} oninput={(e) => set(+e.currentTarget.value)} class="w-40 accent-brand-600" />
    <span class="font-mono text-slate-900 min-w-[5.5rem]">{@html text}</span>
  </label>
{/snippet}
{#snippet mediumCtl()}
  <div class="flex flex-wrap items-center gap-3 text-xs">
    <div class="flex gap-1 bg-slate-100 p-1 rounded-xl border border-slate-200">
      {#each [["gas", "vacuum"], ["sol", "solution"], ["both", "both"]] as [k, l]}
        <button class={seg(medium === k)} disabled={!sm && k !== "gas"} onclick={() => (medium = k)}>{l}</button>
      {/each}
    </div>
    {#if sm && medium !== "gas"}
      {@render slider("ε of the solvent mixture", eps, 1.5, 12, 0.01, (v) => (eps = v), `${eps.toFixed(2)}${solvName() ? ` · ${solvName()}` : ""}`)}
      <div class="flex gap-1">
        {#each SOLVENTS as [n, e]}<button class="px-2 py-1 rounded-lg border text-[11px] font-bold {Math.abs(eps - e) < 0.005 ? 'bg-teal-50 border-teal-300 text-teal-800' : 'border-slate-200 text-slate-500 hover:text-slate-800'}" onclick={() => (eps = e)}>{n}</button>{/each}
      </div>
    {/if}
  </div>
{/snippet}
{#snippet section(title, subtitle)}
  <div class="flex items-baseline gap-3 flex-wrap border-b border-slate-200 pb-2 mt-2">
    <h3 class="text-base font-bold text-slate-900">{@html title}</h3>
    <span class="text-xs text-slate-500">{@html subtitle}</span>
  </div>
{/snippet}
{#snippet equations(rows)}
  <div class="bg-slate-50 border border-slate-200 rounded-xl px-4 py-3 flex flex-col gap-2 self-start max-w-full overflow-x-auto">
    {#each rows as [tag, eq, note]}
      <div class="flex flex-wrap items-baseline gap-x-3">
        <span class="text-xs font-bold text-slate-400 w-6">{tag}</span>
        <span class="text-base md:text-lg text-slate-900">{@html eq}</span>
        {#if note}<span class="text-xs text-slate-500">{@html note}</span>{/if}
      </div>
    {/each}
  </div>
{/snippet}
{#snippet panel(title, caption)}
  <h4 class="text-xs font-bold text-slate-500 uppercase tracking-widest mb-1">{@html title}</h4>
  {#if caption}<p class="text-[11px] text-slate-500 mb-1 max-w-4xl leading-relaxed">{@html caption}</p>{/if}
{/snippet}

{#if status === "loading"}
  <div class="flex items-center justify-center py-20 text-brand-600 text-sm font-bold">
    <div class="w-6 h-6 border-4 border-brand-500 border-t-transparent rounded-full animate-spin mr-3"></div>Loading…
  </div>
{:else if status === "error"}
  <div class="py-20 text-center text-slate-500 text-sm">Could not load {detUrl}.</div>
{:else if model}
  <div class="flex flex-col gap-6 pb-6">
    <!-- shared settings -->
    <div class="flex flex-col gap-2">
      <div class="flex flex-wrap items-center gap-4 text-xs">
        <div class="flex gap-1 bg-slate-100 p-1 rounded-xl border border-slate-200">
          {#each Object.entries(UNITS) as [k, u]}<button class={seg(unit === k)} onclick={() => (unit = k)}>{u[1]}</button>{/each}
        </div>
        {#if model.hasConf}
          <span class="font-bold text-slate-400">S<sub>conf</sub></span>
          <div class="flex gap-1 bg-slate-100 p-1 rounded-xl border border-slate-200">
            {#each [["both", "both bounds"], ["ideal", "ideal count"], ["none", "zero"]] as [k, l]}<button class={seg(confMode === k)} onclick={() => (confMode = k)}>{l}</button>{/each}
          </div>
        {/if}
        {#if complexes}{@render slider("OAm binding correction", delta, -0.4, 0.2, 0.01, (v) => (delta = v), `${fe(delta)} ${ul()} per OAm`)}{/if}
      </div>
      {@render mediumCtl()}
      <p class="text-[11px] text-slate-500 leading-relaxed max-w-4xl">
        Free energies: MACE-MH-1{sm ? "; solution: COSMO on the DFT PBE potential (the dot, the free molecules), ε of the solvent mixture" : " (vacuum)"}.
        {#if complexes}Oleylamine (OAm) binds dissolved {@html unitHtml} as {@html complexes.map((c) => html(c.label)).join(", ")} (MACE-MH-1, vacuum); the binding correction
          adds a constant per bound OAm for what vacuum misses (it overestimates the entropy the amine loses).{:else}No oleylamine complexes computed for {@html unitHtml}.{/if}
        Dashed: vacuum; solid: solution; faded: without oleylamine.
      </p>
    </div>

    <!-- formation -->
    {@render section("Formation of the dot as built", "which conditions make this dot, with this shell, stable against its monomers")}
    <p class="text-sm text-slate-700 max-w-4xl leading-relaxed">
      The Thermodynamics tab looks at the dot alone (against bulk {@html coreHtml}, the free molecules, or its own {@html unitHtml} in a closed system).
      Here the dot is in a synthesis: the {@html coreHtml} monomers come from a precursor that is being converted, the {@html unitHtml} is in excess and partly
      bound by oleylamine. The question is how supersaturated the monomers must be for this exact dot to be stable, and whether it keeps its shell there.
      No claim about which dot a synthesis actually makes.
    </p>
    {@render equations([
      ["(1)", `${unitHtml} + ${precursorHtml} → ${coreHtml} + byproducts`, `irreversible; ${X} precursor e.g. ${PRECURSOR[X] || "the anion source"}; ξ = fraction of it converted (0 to 1), a progress variable, not a rate`],
      ["(2)", `${model.m} ${coreHtml} + ${model.k0} ${unitHtml} ⇌ ${dotHtml}`, `fast; ${coreHtml} enters as its supersaturation against bulk ${coreHtml}, S = c / c<sub>sat</sub>`],
      ...(complexes ? [["(3)", `${unitHtml} + p OAm ⇌ ${unitHtml}(OAm)<sub>p</sub>`, "fast; K<sub>p</sub> = exp(−ΔG<sub>p</sub>/k<sub>B</sub>T)"]] : []),
      ["", `c<sup>free</sup><sub>${html(model.unit)}</sub> = (C<sub>${html(model.unit)}</sub> − ξ C<sub>${X}</sub>) / (1 + Σ<sub>p</sub> K<sub>p</sub> c<sub>OAm</sub><sup>p</sup>)`, `mass balance: each converted ${X} uses one ${unitHtml}`],
      ["", `log<sub>10</sub> S* = [G(dot) − ${model.m} G(${coreHtml}, bulk) − ${model.k0} μ(${unitHtml})] / (${model.m} k<sub>B</sub>T ln 10)`, "the supersaturation at which the dot is in equilibrium with its monomers (dots at 1 µM)"]])}
    {#if hasFormation(model)}
      <div class="flex flex-wrap gap-x-6 gap-y-2">
        {@render slider("T", fT, 300, 650, 25, (v) => (fT = v), `${fT} K`)}
        {@render slider(`total ${unitHtml}`, fLogCL, -4, 0.3, 0.05, (v) => (fLogCL = v), cTxt(10 ** fLogCL))}
        {@render slider(`total ${precursorHtml}`, fLogCX, -4, 0.3, 0.05, (v) => (fLogCX = v), cTxt(10 ** fLogCX))}
        {@render slider("ξ, converted", fXi, 0, 1, 0.01, (v) => (fXi = v), `${(100 * fXi).toFixed(0)}%`)}
        {#if complexes}{@render slider("oleylamine", fOAm, 0, C_NEAT, 0.02, (v) => (fOAm = v), fOAm >= C_NEAT - 1e-9 ? "3.04 M (neat)" : `${fOAm.toFixed(2)} M`)}{/if}
      </div>
      {#if fLeft <= 0}
        <div class="text-xs text-amber-800 bg-amber-50 border border-amber-200 rounded-lg px-3 py-2">The conversion uses up all the {@html unitHtml}: ξ C<sub>{X}</sub> ≥ C<sub>{@html unitHtml}</sub>. Lower ξ or the {X} precursor.</div>
      {:else}
        <div class="flex flex-wrap gap-4 text-xs text-slate-700">
          <span>{@html unitHtml} left: <b class="font-mono">{@html cTxt(fLeft)}</b>, free: <b class="font-mono">{@html cTxt(fFree(fI))}</b>
            {#if speciesF}<span class="text-slate-500">({@html speciesF.species.filter((s) => s.frac > 0.005).map((s) => `${s.label === "free" ? "free" : html(s.label)} ${(100 * s.frac).toFixed(0)}%`).join(", ")})</span>{/if}</span>
          {#each fRows as r}<span class="flex items-center gap-1.5"><svg width="16" height="6"><line x1="0" x2="16" y1="3" y2="3" stroke={CONF[r.c].color} stroke-width="2.5" stroke-dasharray={r.m.dash} /></svg>{@html comboLabel(r)}: log<sub>10</sub> S* = <b class="font-mono">{r.S.toFixed(2)}</b>, ⟨k⟩ = <b class="font-mono">{fk(r.k)}</b></span>{/each}
        </div>
        <div>
          {@render panel("Critical supersaturation of the dot as built", "Above S*, the monomers assemble into this dot; below, it dissolves.")}
          <Plot layers={sStarLayers} xDomain={F_DOM} yDomain={sStarDom} xTicks={fTicks} xLabel="T (K)" yLabel={"log_{10} S*"} height={320}
                hover={sStarHover} ariaLabel="critical supersaturation" legend={legend(complexes && fOAm > 0 ? [{ color: "#94a3b8", label: "without oleylamine" }] : [])} />
        </div>
        <div>
          {@render panel("Shell during the synthesis", `⟨k⟩ against the free ${unitHtml} left after conversion and complexation (an open system: the excess is a reservoir).`)}
          <Plot layers={fShellLayers} xDomain={F_DOM} yDomain={kDomOf(fShellLayers)} xTicks={fTicks} xLabel="T (K)" yLabel="⟨k⟩" height={320}
                hover={lineHover(fShellLayers, (x) => `T = ${Math.round(x)} K`)} ariaLabel="shell during synthesis"
                legend={legend(complexes && fOAm > 0 ? [{ color: "#94a3b8", label: "without oleylamine" }] : [])} />
        </div>
      {/if}
    {:else}
      <div class="text-xs text-slate-500">No props/stability.json for this dot: the formation free energies are missing.</div>
    {/if}

    <!-- post-synthesis -->
    {@render section("After the synthesis", `the dots exist: every dissolved ${unitHtml} comes from them, or is added on purpose`)}
    {@render equations([
      ["(1)", `${dotHtml} ⇌ [${coreHtml}]<sub>${model.m}</sub>[${unitHtml}]<sub>k</sub> + (${model.k0} − k) ${unitHtml}`, "the shell exchanges units with the solution"],
      ...(complexes ? [["(2)", `${unitHtml} + p OAm ⇌ ${unitHtml}(OAm)<sub>p</sub>`, "residual or added oleylamine"]] : []),
      ["", `(k<sub>0</sub> − ⟨k⟩) c<sub>dot</sub> + C<sub>added</sub> = c<sup>free</sup><sub>${html(model.unit)}</sub> (1 + Σ<sub>p</sub> K<sub>p</sub> c<sub>OAm</sub><sup>p</sup>)`, `mass balance; C<sub>added</sub>: ${unitHtml} added as a post-treatment`]])}
    <div class="flex flex-wrap gap-x-6 gap-y-2">
      {@render slider("T", pT, 50, 800, 25, (v) => (pT = v), `${pT} K`)}
      {@render slider("dots", pLogDot, -9, -3, 0.25, (v) => (pLogDot = v), cTxt(10 ** pLogDot))}
      {#if complexes}{@render slider("oleylamine", pOAm, 0, C_NEAT, 0.02, (v) => (pOAm = v), pOAm >= C_NEAT - 1e-9 ? "3.04 M (neat)" : pOAm === 0 ? "none" : `${pOAm.toFixed(2)} M`)}{/if}
      {@render slider(`added ${unitHtml}`, pLogAdd, -12, -1, 0.25, (v) => (pLogAdd = v), pAdd ? cTxt(pAdd) : "none")}
      {@render slider("each wash removes", pWash, 0.5, 0.999, 0.001, (v) => (pWash = v), `${(100 * pWash).toFixed(1)}%`)}
    </div>
    <div class="flex flex-wrap gap-4 text-xs text-slate-700">
      {#each pRows as r}<span class="flex items-center gap-1.5"><svg width="16" height="6"><line x1="0" x2="16" y1="3" y2="3" stroke={CONF[r.c].color} stroke-width="2.5" stroke-dasharray={r.m.dash} /></svg>{@html comboLabel(r)}: ⟨k⟩ = <b class="font-mono">{fk(r.kmean)}</b>, free {@html unitHtml} <b class="font-mono">{@html cTxt(r.cFree)}</b></span>{/each}
      {#if speciesP && pOAm > 0}<span class="text-slate-500">(dissolved {@html unitHtml}: {@html speciesP.species.filter((s) => s.frac > 0.005).map((s) => `${s.label === "free" ? "free" : html(s.label)} ${(100 * s.frac).toFixed(0)}%`).join(", ")})</span>{/if}
    </div>
    {#if firstUnit}
      {@const lose0 = firstUnit.dG < firstUnit.thr}
      {@const lose1 = firstUnit.dGc < firstUnit.thr}
      <div class="text-xs text-slate-700 bg-slate-50 border border-slate-200 rounded-lg px-3 py-2 max-w-4xl leading-relaxed">
        <b>First unit, vacuum</b> (nothing added): it leaves if its cost is below k<sub>B</sub>T ln(1 M / c<sub>dot</sub>){#if pOAm > 0 && complexes} + k<sub>B</sub>T ln(1 + Σ<sub>p</sub> K<sub>p</sub> c<sub>OAm</sub><sup>p</sup>){/if} = <b class="font-mono">{fe(firstUnit.thr)} {ul()}</b>.
        Cost ΔG = <b class="font-mono">{fe(firstUnit.dG)}</b>{#if model.hasConf}, ΔG − TΔS<sub>conf</sub> = <b class="font-mono">{fe(firstUnit.dGc)}</b>{/if} {ul()}:
        <span class="font-bold {lose0 || lose1 ? 'text-rose-700' : 'text-emerald-700'}">{lose0 && lose1 ? "the first unit leaves" : !lose0 && !lose1 ? "the first unit stays" : "borderline: it depends on the configurational entropy"}</span>.
        A quick test; the curves solve the whole ladder.
      </div>
    {/if}
    <div>
      {@render panel("Dilution", `⟨k⟩ against the dot concentration: the units that leave stay in solution, free or bound to oleylamine.`)}
      <Plot layers={dilLayers} xDomain={[-9, -3]} yDomain={kDomOf(dilLayers)}
            xTicks={[-9, -8, -7, -6, -5, -4, -3].map((v) => ({ v, label: v === -6 ? "1 µM" : v === -3 ? "1 mM" : v === -9 ? "1 nM" : pow10(v) }))}
            xLabel="dot concentration (M)" yLabel="⟨k⟩" height={300} hover={lineHover(dilLayers, (x) => `dots ${cTxt(10 ** x)}`)} ariaLabel="dilution" legend={legend()} />
    </div>
    {#if complexes}
      <div>
        {@render panel("Oleylamine", `⟨k⟩ against the oleylamine concentration (the rest of the solvent inert): the amine pulls ${unitHtml} off by binding it.`)}
        <Plot layers={oamLayers} xDomain={[0, C_NEAT]} yDomain={kDomOf(oamLayers)} xLabel="oleylamine (M), 3.04 M = neat" yLabel="⟨k⟩" height={300}
              hover={lineHover(oamLayers, (x) => `OAm ${x.toFixed(2)} M`)} ariaLabel="oleylamine" legend={legend()} />
      </div>
    {/if}
    {#if model.hasUp}
      <div>
        {@render panel(`${unitHtml} post-treatment`, `⟨k⟩ against the ${unitHtml} added to the purified dots: the dot takes up units on its free surface sites.`)}
        <Plot layers={addLayers} xDomain={[-9, -1]} yDomain={kDomOf(addLayers)}
              xTicks={[-9, -7, -5, -3, -1].map((v) => ({ v, label: pow10(v) }))}
              xLabel={`added ${model.unit} (M)`} yLabel="⟨k⟩" height={300} hover={lineHover(addLayers, (x) => `added ${cTxt(10 ** x)}`)} ariaLabel="post-treatment" legend={legend()} />
      </div>
    {/if}
    <div>
      {@render panel("Washing cycles", `Each wash keeps the dots (same concentration, same medium) and removes ${(100 * pWash).toFixed(1)}% of the dissolved ${unitHtml}, free and complexed; the dots then re-equilibrate.`)}
      <Plot layers={washLayers} xDomain={[-0.4, 8.4]} yDomain={kDomOf(washLayers)} xTicks={[0, 1, 2, 3, 4, 5, 6, 7, 8].map((v) => ({ v, label: v ? String(v) : "start" }))}
            xLabel="wash" yLabel="⟨k⟩" height={300} hover={lineHover(washLayers, (x) => (Math.round(x) ? `after wash ${Math.round(x)}` : "start"))} ariaLabel="washing" legend={legend()} />
    </div>
    <div>
      {@render panel("Heating the purified dots", `The dots are heated in the same solvent, nothing added or removed: ⟨k⟩ at each T, at equilibrium. The units that leave stay in solution, so desorption limits itself, and they come back on cooling unless they are washed away.`)}
      <Plot layers={annealLayers} xDomain={[50, 800]} yDomain={kDomOf(annealLayers)} xTicks={tempTicks} xLabel="T (K)" yLabel="⟨k⟩" height={300}
            hover={lineHover(annealLayers, (x) => `T = ${Math.round(x)} K`)} ariaLabel="heating" legend={legend()} />
    </div>
  </div>
{/if}
