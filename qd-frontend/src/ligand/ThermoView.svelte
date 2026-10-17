<script>
  // Thermodynamics tab: stability and the Z-type ladder of a dot [MX]_m[ML_n]_k0, in vacuum, in an implicit
  // solvent, or both on the same plots. Data: props/detachment.json, props/stability.json and the DFT
  // solvation cubes (props/solvation_cubes.json: COSMO on the CP2K PBE potential). Formulas in model.js.
  // The ladder is a closed system: purified dots in a solvent, every free ML_n comes from the dots, so the
  // knobs are T, the dot concentration and the solvent; adding ML_n is a post-treatment (Synthesis tab).
  import Plot from "./Plot.svelte";
  import { KB, UNITS, buildModel, solventModel, unitCosts, decompositionPerCore, bindingPerUnit, siteCosts,
           solveClosed, modelAt, siteGroup, closedKmeanGrid } from "./model.js";

  let { detUrl = null, stabUrl = null, cubesUrl = null } = $props();

  let model = $state(null);
  let sm = $state(null);                  // solvent model (null: vacuum only)
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
  let confMode = $state("both");          // "both" | "ideal" | "none"
  let medium = $state("gas");             // "gas" | "sol" | "both"
  let eps = $state(2.38);
  const SOLVENTS = [["hexane", 1.88], ["toluene", 2.38], ["chloroform", 4.81], ["THF", 7.58]];

  const iTof = (t) => Math.max(0, model.T.findIndex((v) => Math.abs(v - t) < 1e-6));
  const nearestT = (x) => model.T.reduce((b, t, j) => (Math.abs(t - x) < Math.abs(model.T[b] - x) ? j : b), 0);
  const uf = () => UNITS[unit][0];
  const ul = () => UNITS[unit][1];
  const fe = (v, d = UNITS[unit][2]) => (isFinite(v) ? (v * uf() >= 0 ? "+" : "−") + Math.abs(v * uf()).toFixed(d) : "–");
  const fk = (v) => (isFinite(v) ? v.toFixed(1) : "–");
  const sub = (s) => String(s).replace(/\d/g, (d) => "₀₁₂₃₄₅₆₇₈₉"[d]);
  const pow10 = (e) => "10" + String(e).replace("-", "⁻").replace(/\d/g, (d) => "⁰¹²³⁴⁵⁶⁷⁸⁹"[d]);
  const html = (s) => String(s).replace(/(\d+)/g, "<sub>$1</sub>");
  const cTxt = (c) => {
    if (!(c > 0)) return "0";
    const e = Math.floor(Math.log10(c)), m = c / 10 ** e;
    return e >= -2 && e <= 0 ? `${+c.toPrecision(2)} M` : `${m.toFixed(1)}·10<sup>${e}</sup> M`;
  };
  const cDotTicks = [-9, -8, -7, -6, -5, -4, -3].map((v) => ({ v, label: v === -6 ? "1 µM" : v === -3 ? "1 mM" : v === -9 ? "1 nM" : pow10(v) }));
  const confs = () => (!model?.hasConf ? [false] : confMode === "both" ? [true, false] : [confMode === "ideal"]);
  const CONF = { true: { color: "#eb6834", label: "ideal S<sub>conf</sub>" }, false: { color: "#2a78d6", label: "S<sub>conf</sub> = 0" } };
  const solvName = () => SOLVENTS.find(([, e]) => Math.abs(e - eps) < 0.005)?.[0];
  // media on the plots: dashed = vacuum, solid = solution; bars and single curves take the medium's colour
  let media = $derived.by(() => {
    const gas = { key: "gas", sv: null, dash: "6 4", color: "#64748b", label: "vacuum" };
    const sol = { key: "sol", sv: sm ? { sm, eps } : null, dash: "", color: "#0f766e",
                  label: `solution, ε = ${eps.toFixed(2)}${solvName() ? ` (${solvName()})` : ""}` };
    return !sm || medium === "gas" ? [gas] : medium === "sol" ? [sol] : [gas, sol];
  });
  const combos = () => media.flatMap((m) => confs().map((c) => ({ m, c })));
  const comboLabel = ({ m, c }) => `${m.label}, ${CONF[c].label}`;
  const legendCombos = () => [
    ...confs().map((c) => ({ color: CONF[c].color, label: CONF[c].label })),
    ...media.map((m) => ({ color: "#334155", dash: m.dash, label: m.label }))];
  const mediaLegend = () => media.map((m) => ({ color: m.color, dash: m.dash, label: m.label }));

  // ---------------- labels ----------------
  let dotHtml = $derived(model ? `[${html(model.core || "MX")}]<sub>${model.m ?? "m"}</sub>[${html(model.unit)}]<sub>${model.k0 ?? "k"}</sub>` : "");
  let unitHtml = $derived(model ? html(model.unit) : "");
  let coreHtml = $derived(model ? html(model.core || "MX") : "");
  let kLabel = $derived(model ? `k, ${sub(model.unit)} units on the dot` : "");
  let kTicks = $derived.by(() => {
    if (!model) return [];
    const a = model.states[0].k, b = model.states.at(-1).k, step = b - a > 16 ? 2 : 1, out = [];
    for (let k = a; k <= b; k++) if ((k - model.k0) % step === 0) out.push({ v: k, label: String(k) });
    return out;
  });
  let kDom = $derived(model ? [model.states[0].k - 0.7, model.states.at(-1).k + 0.7] : [0, 1]);
  const atLeast = (k) => model && ((model.capDown && k <= model.states[0].k + 0.5) || (model.capUp && model.hasUp && k >= model.states.at(-1).k - 0.5));
  const tempTicks = [100, 200, 300, 400, 500, 600, 700, 800].map((v) => ({ v, label: String(v) }));
  const domOf = (vals, extra = [0]) => {
    const ys = vals.filter(isFinite).concat(extra), a = Math.min(...ys), b = Math.max(...ys), p = 0.08 * (b - a || 1);
    return [a - p, b + p];
  };

  // =====================================================================================
  // Stability: written with the dot on the left, positive = the dot is stable against that reaction
  // =====================================================================================
  let sLogC = $state(0);                  // monomer concentration (M), log10
  let hasStab = $derived(!!(model?.G_dot && model?.G_bulk && model?.G_MLn && model?.G_MX));
  const tLine = (fn) => media.map((m) => ({ m, pts: model.T.map((t, i) => [t, fn(i, m.sv) * uf()]) }));
  let decCurves = $derived(hasStab ? tLine((i, sv) => -decompositionPerCore(model, i, 10 ** sLogC, sv)) : []);
  let dissCurves = $derived(hasStab ? tLine((i, sv) => -bindingPerUnit(model, i, 10 ** sLogC, sv)) : []);
  // dissolution temperature: the dissociation free energy crosses zero
  const tZero = (pts) => {
    for (let j = 1; j < pts.length; j++) if (pts[j - 1][1] > 0 && pts[j][1] <= 0)
      return pts[j - 1][0] + (pts[j][0] - pts[j - 1][0]) * pts[j - 1][1] / (pts[j - 1][1] - pts[j][1]);
    return null;
  };
  const curveLayers = (curves, zeroMarks = false) => [
    ...curves.map(({ m, pts }) => ({ type: "line", pts, color: m.color, dash: m.dash })),
    { type: "hline", y: 0, color: "#94a3b8", dash: "" },
    ...(zeroMarks ? curves.map(({ pts }) => tZero(pts)).map((x, j) => (x ? { type: "vline", x, color: curves[j].m.color, label: `T_{diss} ${Math.round(x)} K` } : null)).filter(Boolean) : [])];
  const curveHover = (curves, label) => (x) => {
    if (!curves.length) return null;
    const i = nearestT(x);
    return { x: model.T[i], title: `T = ${model.T[i]} K`, rows: curves.map(({ m, pts }) => ({ color: m.color, dash: m.dash, text: `${m.label}: ${label} ${fe(pts[i][1] / uf())} ${ul()}` })) };
  };

  // =====================================================================================
  // Ladder: closed system (T, dot concentration, solvent)
  // =====================================================================================
  let lT = $state(300);
  let lLogDot = $state(-6);               // dot concentration (M), log10
  let lI = $derived(model ? iTof(lT) : 0);
  const closedAt = (i, cdot, c, sv) => solveClosed(model, i, cdot, 0, 0, c, null, sv);
  let pops = $derived(model ? combos().map((cb) => ({ ...cb, ...closedAt(lI, 10 ** lLogDot, cb.c, cb.m.sv) })) : []);
  let costs = $derived(model ? media.map((m) => ({ m, cost: unitCosts(model, lI, m.sv) })) : []);
  // the line units must lie below to leave: k_BT ln(1 M / c_free), c_free from the mass balance (S_conf = 0)
  let thresholds = $derived(model ? media.map((m) => ({ m, y: -closedAt(lI, 10 ** lLogDot, false, m.sv).dmu })) : []);
  const thrLayers = (f) => thresholds.map(({ m, y }) => ({ type: "hline", y: y * f, color: "#0f172a", dash: m.key === "gas" && media.length > 1 ? "2 4" : "6 4", labelAt: m.key === "gas" && media.length > 1 ? "start" : "end",
    label: `k_BT ln(1 M / c_{free})${media.length > 1 ? ` (${m.key === "gas" ? "vacuum" : "solution"})` : ""} = ${fe(y)} ${ul()}` }));

  let costLayers = $derived.by(() => {
    if (!model) return [];
    const kT = KB * model.T[lI], f = uf(), n = costs.length, w = n > 1 ? 0.34 : 0.62;
    const L = [];
    costs.forEach(({ m, cost }, j) => {
      const off = n > 1 ? (j ? 0.18 : -0.18) : 0;
      L.push({ type: "bars", bars: model.units.map((u, q) => ({ x: u.k - 0.5 + off, y0: 0, y1: cost[q] * f, w, color: m.color, opacity: 0.85 })) });
      if (model.hasConf && confMode !== "none")
        L.push({ type: "ticks", ticks: model.units.map((u, q) => ({ x: u.k - 0.5 + off, y: (cost[q] - kT * u.dS) * f, w: w + 0.08, color: CONF.true.color })) });
    });
    L.push({ type: "points", pts: model.units.map((u) => [u.k - 0.5, u.dE * f]), color: "#334155", hollow: true, r: 3 });
    L.push(...thrLayers(f));
    L.push({ type: "hline", y: 0, color: "#94a3b8", dash: "" });
    L.push({ type: "vline", x: model.k0, label: "as built" });
    return L;
  });
  let costDom = $derived(model ? domOf([...costs.flatMap((c) => c.cost), ...model.units.map((u) => u.dE),
    ...costs.flatMap((c) => c.cost.map((v, q) => v - KB * model.T[lI] * model.units[q].dS))], [0, ...thresholds.map((t) => t.y)]).map((v) => v * uf()) : [0, 1]);
  const costHover = (x) => {
    const q = model.units.findIndex((u) => Math.abs(u.k - 0.5 - x) < 0.5);
    if (q < 0) return null;
    const u = model.units[q], kT = KB * model.T[lI];
    const which = u.side === "removal" ? `removal ${u.order}` : `addition ${u.order} (its desorption)`;
    return { x: u.k - 0.5, title: `unit between k = ${u.k - 1} and ${u.k} · ${which}`,
             rows: [{ text: `ΔE (vacuum, electronic) ${fe(u.dE)} ${ul()}`, color: "#334155", dash: "2 2" },
                    ...costs.map(({ m, cost }) => ({ color: m.color, text: `${m.label}: ΔG ${fe(cost[q])}${model.hasConf ? `, − TΔS<sub>conf</sub> ${fe(cost[q] - kT * u.dS)}` : ""} ${ul()}` }))] };
  };

  let profLayers = $derived.by(() => {
    if (!model) return [];
    const f = uf();
    return pops.map((P) => ({ type: "line", pts: model.states.map((s, j) => [s.k, P.omega[j] * f]), color: CONF[P.c].color, dash: P.m.dash }))
      .concat([{ type: "vline", x: model.k0, label: "as built" }]);
  });
  let profDom = $derived.by(() => {
    if (!pops.length) return [0, 1];
    const ys = pops.flatMap((P) => P.omega), a = Math.min(...ys), b = Math.min(Math.max(...ys), a + 12);
    return [(a - 0.05 * (b - a + 0.1)) * uf(), (b + 0.05 * (b - a + 0.1)) * uf()];
  });
  let popLayers = $derived.by(() => {
    if (!model) return [];
    const n = pops.length, w = 0.8 / n;
    return pops.map((P, i) => ({ type: "bars", bars: model.states.map((s, j) => ({
      x: s.k - 0.4 + w * (i + 0.5), y0: 0, y1: P.p[j], w: w * 0.9, color: CONF[P.c].color, opacity: P.m.key === "gas" && media.length > 1 ? 0.35 : 0.85 })) }))
      .concat([{ type: "vline", x: model.k0, label: "as built" }]);
  });
  const popHover = (x) => {
    const j = model.states.findIndex((s) => Math.abs(s.k - x) < 0.5);
    if (j < 0) return null;
    const s = model.states[j];
    return { x: s.k, title: `k = ${s.k}${atLeast(s.k) ? " (search limit)" : ""}`,
             rows: pops.map((P) => ({ color: CONF[P.c].color, dash: P.m.dash, text: `${comboLabel(P)}: P = ${P.p[j].toFixed(3)}, Ω − Ω<sub>k0</sub> = ${fe(P.omega[j])} ${ul()}` })) };
  };

  // dilution: <k> against the dot concentration at T
  const DIL = Array.from({ length: 25 }, (_, i) => -9 + 0.25 * i);
  let dilLayers = $derived.by(() => {
    if (!model) return [];
    return combos().map(({ m, c }) => ({ type: "line", color: CONF[c].color, dash: m.dash,
      pts: DIL.map((l) => [l, closedAt(lI, 10 ** l, c, m.sv).kmean]) }))
      .concat([{ type: "hline", y: model.k0, color: "#334155", label: `as built, k = ${model.k0}` }, { type: "vline", x: lLogDot }]);
  });
  let dilDom = $derived.by(() => {
    if (!model) return [0, 1];
    const ys = dilLayers.filter((l) => l.pts).flatMap((l) => l.pts.map((p) => p[1])).concat([model.k0 - 1, model.k0 + 0.5]);
    return [Math.max(model.states[0].k, Math.floor(Math.min(...ys))) - 0.5, Math.min(model.states.at(-1).k, Math.ceil(Math.max(...ys))) + 0.5];
  });
  const dilHover = (x) => {
    const l = Math.round(x * 4) / 4, t = model.T[lI];
    return { x: l, title: `dots ${cTxt(10 ** l)}, ${t} K`,
             rows: combos().map((cb) => { const r = closedAt(lI, 10 ** l, cb.c, cb.m.sv);
               return { color: CONF[cb.c].color, dash: cb.m.dash, text: `${comboLabel(cb)}: ⟨k⟩ = ${fk(r.kmean)}, free ${unitHtml} ${cTxt(r.cFree)}` }; }) };
  };

  // <k> maps over T and the dot concentration, one per medium: T interpolated between the grid points
  const MAP_T = Array.from({ length: 76 }, (_, i) => 50 + 10 * i);
  const MAP_C = Array.from({ length: 61 }, (_, i) => -9 + 0.1 * i);
  let mapConf = $derived(model?.hasConf && confMode !== "none");
  function kRGB(dk, dn, up) {
    const a = dk < 0 ? Math.min(1, -dk / Math.max(dn, 1)) : Math.min(1, dk / Math.max(up, 1));
    const g = [241, 240, 236], to = dk < 0 ? [214, 69, 65] : [42, 120, 214], e = Math.pow(a, 0.75);
    return g.map((v, j) => Math.round(v + (to[j] - v) * e));
  }
  const kColor = (dk, dn, up) => `rgb(${kRGB(dk, dn, up).join(",")})`;
  let maps = $derived.by(() => {
    if (!model || typeof document === "undefined") return [];
    const k0 = model.k0, dn = k0 - model.states[0].k, up = model.states.at(-1).k - k0;
    return media.map((m) => {
      const cv = document.createElement("canvas");
      cv.width = MAP_T.length; cv.height = MAP_C.length;
      const ctx = cv.getContext("2d"), img = ctx.createImageData(cv.width, cv.height);
      const grid = closedKmeanGrid(model, MAP_T, MAP_C, mapConf, m.sv);
      MAP_T.forEach((_, i) => MAP_C.forEach((l, j) => {
        const [r, g, b] = kRGB(grid[i * MAP_C.length + j] - k0, dn, up), o = 4 * ((cv.height - 1 - j) * cv.width + i);
        img.data[o] = r; img.data[o + 1] = g; img.data[o + 2] = b; img.data[o + 3] = 255;
      }));
      ctx.putImageData(img, 0, 0);
      return { m, href: cv.toDataURL() };
    });
  });
  const mapHover = (m) => (x) => {
    const t = Math.round(x / 5) * 5, mt = modelAt(model, t);
    return { x: t, title: `${m.label}, T = ${t} K`,
             rows: [-9, -7, -6, -5, -3].map((l) => {
               const r = solveClosed(mt, 0, 10 ** l, 0, 0, mapConf, null, m.sv);
               return { text: `dots ${cTxt(10 ** l)}: ⟨k⟩ = ${fk(r.kmean)}` };
             }) };
  };

  // site map: removal free energy of each symmetry-distinct unit of the intact dot, grouped by where it sits
  let sT = $state(300);
  let sLogDot = $state(-6);
  let sI = $derived(model ? iTof(sT) : 0);
  let siteOrder = $derived.by(() => {
    if (!model) return [];
    const cls = model.siteClasses.map((c, i) => ({ i, g: siteGroup(c), dE: c.dE }));
    const gmin = new Map();
    for (const c of cls) gmin.set(c.g, Math.min(gmin.get(c.g) ?? Infinity, c.dE));
    return cls.sort((a, b) => gmin.get(a.g) - gmin.get(b.g) || a.g.localeCompare(b.g) || a.dE - b.dE);
  });
  let siteVals = $derived(model ? media.map((m) => ({ m, g: siteCosts(model, sI, m.sv) })) : []);
  let siteThr = $derived(model ? media.map((m) => ({ m, y: -closedAt(sI, 10 ** sLogDot, false, m.sv).dmu })) : []);
  let siteLayers = $derived.by(() => {
    if (!model?.siteClasses.length) return [];
    const f = uf(), n = siteVals.length, w = n > 1 ? 0.36 : 0.66, L = [];
    siteOrder.forEach((c, r) => { if (r && siteOrder[r - 1].g !== c.g) L.push({ type: "vline", x: r + 0.5, color: "#94a3b8", dash: "" }); });
    siteVals.forEach(({ m, g }, j) => L.push({ type: "bars", bars: siteOrder.map((c, r) => ({
      x: r + 1 + (n > 1 ? (j ? 0.19 : -0.19) : 0), y0: 0, y1: g[c.i] * f, w, color: m.color, opacity: model.siteClasses[c.i].molecular === false ? 0.55 : 0.9 })) }));
    siteThr.forEach(({ m, y }) => L.push({ type: "hline", y: y * f, color: "#0f172a", dash: m.key === "gas" && media.length > 1 ? "2 4" : "6 4", labelAt: m.key === "gas" && media.length > 1 ? "start" : "end",
      label: `k_BT ln(1 M / c_{free})${media.length > 1 ? ` (${m.key === "gas" ? "vacuum" : "solution"})` : ""}` }));
    L.push({ type: "hline", y: 0, color: "#94a3b8", dash: "" });
    return L;
  });
  let siteGroups = $derived.by(() => {
    const out = [];
    siteOrder.forEach((c, r) => { if (!out.length || out.at(-1).g !== c.g) out.push({ g: c.g, a: r + 1, b: r + 1 }); else out.at(-1).b = r + 1; });
    return out;
  });
  const siteHover = (x) => {
    const r = Math.round(x) - 1, c = siteOrder[r];
    if (!c) return null;
    const cl = model.siteClasses[c.i];
    return { x: r + 1, title: `${c.g} · ×${cl.mult} equivalent${cl.molecular === false ? " · non-molecular unit" : ""}`,
             rows: [{ text: `ΔE ${fe(cl.dE)} ${ul()}` }, ...siteVals.map(({ m, g }) => ({ color: m.color, dash: m.dash, text: `${m.label}: ΔG(${model.T[sI]} K, 1 M) ${fe(g[c.i])} ${ul()}` }))] };
  };

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
        <button class={seg(medium === k)} disabled={!sm && k !== "gas"} title={!sm && k !== "gas" ? "No DFT solvation for this dot yet" : ""}
                onclick={() => (medium = k)}>{l}</button>
      {/each}
    </div>
    {#if sm && medium !== "gas"}
      {@render slider("ε", eps, 1.5, 12, 0.01, (v) => (eps = v), `${eps.toFixed(2)}${solvName() ? ` · ${solvName()}` : ""}`)}
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
{#snippet reaction(text)}
  <div class="text-base md:text-lg text-slate-900 font-medium bg-slate-50 border border-slate-200 rounded-xl px-4 py-2 self-start">{@html text}</div>
{/snippet}
{#snippet panel(title, caption)}
  <h4 class="text-xs font-bold text-slate-500 uppercase tracking-widest mb-1">{@html title}</h4>
  {#if caption}<p class="text-[11px] text-slate-500 mb-1 max-w-4xl leading-relaxed">{@html caption}</p>{/if}
{/snippet}

{#if status === "loading"}
  <div class="flex items-center justify-center h-full text-brand-600 text-sm font-bold">
    <div class="w-6 h-6 border-4 border-brand-500 border-t-transparent rounded-full animate-spin mr-3"></div>Loading…
  </div>
{:else if status === "error"}
  <div class="flex items-center justify-center h-full text-slate-500 text-sm">Could not load {detUrl}.</div>
{:else if model}
  <div class="flex flex-col gap-6 pb-6">
    <!-- header: the dot and the shared settings -->
    <div class="flex flex-col gap-2">
      <div class="text-sm text-slate-700">
        <span class="font-bold text-slate-900">{@html dotHtml}</span> · Z-type unit {@html unitHtml} ·
        ladder from k = {model.states[0].k} to {model.states.at(-1).k}
        ({model.nDown} removals{model.hasUp ? `, ${model.nUp} additions` : ", no addition site"})
      </div>
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
      </div>
      <p class="text-[11px] text-slate-500 leading-relaxed max-w-4xl">
        Free energies of MACE-MH-1 (harmonic vibrations, rigid rotor, ideal solutes at 1 M) for every species.
        {#if sm}
          Solution: ΔG<sub>solv</sub>(ε) = (1 − 1/ε) E<sub>∞</sub> + (α+γ)(ε) S, COSMO on the DFT PBE potential with the cavity on the
          density isosurface (cubes of the intact, the fully stripped{sm.upExtrapolated ? "" : " and the fully added"} dot and the monomers;
          the other ladder states linear in k{sm.upExtrapolated ? "; the additions extend the removal line" : ""}). Dashed: vacuum; solid: solution.
        {:else}
          No DFT solvation for this dot yet: vacuum only.
        {/if}
      </p>
    </div>

    <!-- stability -->
    {#if hasStab}
      {@render section("Stability", `against bulk ${coreHtml} and against the free molecules`)}
      <div class="flex flex-wrap items-center gap-x-6 gap-y-2">
        {@render mediumCtl()}
        {@render slider(`${coreHtml} and ${unitHtml}`, sLogC, -6, 0, 0.25, (v) => (sLogC = v), cTxt(10 ** sLogC))}
      </div>
      <div class="flex flex-col gap-2">
        {@render panel("Decomposition into bulk", "")}
        {@render reaction(`${dotHtml} → ${model.m} ${coreHtml}(bulk) + ${model.k0} ${unitHtml}`)}
        <p class="text-[11px] text-slate-500 max-w-4xl">ΔG<sub>r</sub> per {@html coreHtml}: positive, the dot is stable against bulk {@html coreHtml} plus dissolved {@html unitHtml}; negative, it would rather grow into the bulk (it is metastable, as small dots are).</p>
        <Plot layers={curveLayers(decCurves)} xDomain={[50, 800]} yDomain={domOf(decCurves.flatMap((c) => c.pts.map((p) => p[1])))} xTicks={tempTicks}
              xLabel="T (K)" yLabel={`ΔG_{r} per ${sub(model.core)} (${ul()})`} height={300} hover={curveHover(decCurves, "ΔG<sub>r</sub>")}
              ariaLabel="decomposition into bulk" legend={mediaLegend()} />
      </div>
      <div class="flex flex-col gap-2">
        {@render panel("Dissociation into free molecules", "")}
        {@render reaction(`${dotHtml} → ${model.m} ${coreHtml} + ${model.k0} ${unitHtml}`)}
        <p class="text-[11px] text-slate-500 max-w-4xl">ΔG<sub>r</sub> per unit ({model.m + model.k0} units): positive, the dot holds together; T<sub>diss</sub> where it crosses zero. The molecular {@html coreHtml} is a reference state, not the species of a synthesis.</p>
        <Plot layers={curveLayers(dissCurves, true)} xDomain={[50, 800]} yDomain={domOf(dissCurves.flatMap((c) => c.pts.map((p) => p[1])))} xTicks={tempTicks}
              xLabel="T (K)" yLabel={`ΔG_{r} per unit (${ul()})`} height={300} hover={curveHover(dissCurves, "ΔG<sub>r</sub>")}
              ariaLabel="dissociation into free molecules" legend={mediaLegend()} />
      </div>
    {/if}

    <!-- ladder -->
    {@render section(`${unitHtml} ladder`, "purified dots in a solvent: every free unit comes from the dots")}
    {@render reaction(`${dotHtml} ⇌ [${coreHtml}]<sub>${model.m}</sub>[${unitHtml}]<sub>${model.k0}−n</sub> + n ${unitHtml}, &nbsp; (k<sub>0</sub> − ⟨k⟩) c<sub>dot</sub> = c<sub>free</sub>`)}
    <div class="flex flex-wrap items-center gap-x-6 gap-y-2">
      {@render mediumCtl()}
      {@render slider("T", lT, 50, 800, 25, (v) => (lT = v), `${lT} K`)}
      {@render slider("dots", lLogDot, -9, -3, 0.25, (v) => (lLogDot = v), cTxt(10 ** lLogDot))}
    </div>
    <div class="flex flex-wrap gap-4 text-xs text-slate-700">
      {#each pops as P}
        <span class="flex items-center gap-1.5"><svg width="16" height="6"><line x1="0" x2="16" y1="3" y2="3" stroke={CONF[P.c].color} stroke-width="2.5" stroke-dasharray={P.m.dash} /></svg>{@html comboLabel(P)}: ⟨k⟩ = <b class="font-mono">{fk(P.kmean)}</b>, free {@html unitHtml} <b class="font-mono">{@html cTxt(P.cFree)}</b></span>
      {/each}
    </div>
    <p class="text-[11px] text-slate-500 max-w-4xl">States above k<sub>0</sub> still appear: some dots take up units that others release. Adding {@html unitHtml} (a post-treatment) is in the Synthesis and post-synthesis tab.</p>
    <div>
      {@render panel("Cost of each unit (desorption, positive = bound)",
        `Bar between k − 1 and k: the free energy to take one ${unitHtml} off a dot with k units (1 M). Units whose bar lies below the dashed line leave: the line is k<sub>B</sub>T ln(1 M / c<sub>free</sub>) with the free concentration of the mass balance.`)}
      <Plot layers={costLayers} xDomain={kDom} yDomain={costDom} xTicks={kTicks} xLabel={kLabel} yLabel={`ΔG to remove the unit (${ul()})`}
            height={360} hover={costHover} ariaLabel="cost of each unit"
            legend={[...mediaLegend().map((g) => ({ ...g, kind: "box" })),
                     ...(model.hasConf && confMode !== "none" ? [{ kind: "tick", color: CONF.true.color, label: "− TΔS<sub>conf</sub> (ideal)" }] : []),
                     { kind: "circle", color: "#334155", label: "ΔE (vacuum, electronic)" }]} />
    </div>
    <div>
      {@render panel("Free-energy profile at this T and dot concentration", "Ω<sub>k</sub> − Ω<sub>k0</sub>: the lowest state is the most populated.")}
      <Plot layers={profLayers} xDomain={kDom} yDomain={profDom} xTicks={kTicks} xLabel={kLabel} yLabel={`Ω_{k} − Ω_{k0} (${ul()})`}
            height={300} hover={popHover} ariaLabel="free energy profile" legend={legendCombos()} />
    </div>
    <div>
      {@render panel("Populations of the shell states", media.length > 1 ? "Light bars: vacuum; dark: solution." : "")}
      <Plot layers={popLayers} xDomain={kDom} yDomain={[0, 1.05]} xTicks={kTicks} xLabel={kLabel} yLabel={"population P_{k}"}
            height={240} hover={popHover} ariaLabel="populations" legend={confs().map((c) => ({ kind: "box", color: CONF[c].color, label: CONF[c].label }))} />
    </div>
    <div>
      {@render panel(`Dilution at ${model.T[lI]} K`, `⟨k⟩ against the dot concentration: the more dilute, the more units leave.`)}
      <Plot layers={dilLayers} xDomain={[-9, -3]} yDomain={dilDom} xTicks={cDotTicks}
            xLabel="dot concentration (M)" yLabel="⟨k⟩" height={320} hover={dilHover} ariaLabel="dilution" legend={legendCombos()} />
    </div>
    {#each maps as mp}
      <div>
        {@render panel(`Mean shell ⟨k⟩ against T and dot concentration · ${mp.m.label}${mapConf ? ", ideal S<sub>conf</sub>" : ""}`,
          "Dashed cross: the T and concentration set above. Red: units lost.")}
        <Plot layers={[{ type: "image", href: mp.href, x0: 45, x1: 805, y0: -9.05, y1: -2.95 }, { type: "vline", x: model.T[lI], color: "#0f172a", dash: "3 3" }, { type: "hline", y: lLogDot, color: "#0f172a", dash: "3 3" }]}
              xDomain={[45, 805]} yDomain={[-9.05, -2.95]} xTicks={tempTicks} yTicks={cDotTicks}
              xLabel="T (K)" yLabel="dot concentration (M)" height={380} hover={mapHover(mp.m)} ariaLabel="mean shell map" />
        <div class="flex flex-wrap items-center gap-3 text-[11px] text-slate-500 mt-1 px-1">
          {#each [model.states[0].k, Math.round((model.states[0].k + model.k0) / 2), model.k0] as k}
            <span class="flex items-center gap-1"><span class="w-3 h-3 rounded-sm inline-block border border-slate-200" style="background:{kColor(k - model.k0, model.k0 - model.states[0].k, model.states.at(-1).k - model.k0)}"></span>⟨k⟩ = {k}{atLeast(k) ? " or fewer" : ""}</span>
          {/each}
        </div>
      </div>
    {/each}
    {#if model.siteClasses.length}
      <div class="flex flex-col gap-2">
        {@render panel("Binding-site map: first removal, per symmetry-distinct unit",
          `ΔG(1 M) of removing each distinct ${unitHtml} unit of the intact dot: its relaxed ΔE plus the thermal part of the first ladder step${media.some((m) => m.key === "sol") ? ", plus the solvation change of that step in solution" : ""}, grouped by where it sits. Faded bars: non-molecular units (completed with a neighbouring ${model.L}).`)}
        <div class="flex flex-wrap items-center gap-x-6 gap-y-2">
          {@render slider("T", sT, 50, 800, 25, (v) => (sT = v), `${sT} K`)}
          {@render slider("dots", sLogDot, -9, -3, 0.25, (v) => (sLogDot = v), cTxt(10 ** sLogDot))}
        </div>
        <div class="flex text-[11px] font-bold text-slate-500 pl-[58px] pr-[14px]">
          {#each siteGroups as g}<div class="text-center border-b-2 border-slate-200 mx-0.5 pb-0.5" style="flex: {g.b - g.a + 1}">{g.g}</div>{/each}
        </div>
        <Plot layers={siteLayers} xDomain={[0.5, model.siteClasses.length + 0.5]} yDomain={domOf(siteVals.flatMap((s) => s.g), [0, ...siteThr.map((t) => t.y)]).map((v) => v * uf())}
              xTicks={siteOrder.map((_, r) => ({ v: r + 1, label: String(r + 1) }))}
              xLabel="symmetry-distinct unit (cheapest first within each group)" yLabel={`ΔG to remove (${ul()})`} height={320} hover={siteHover} ariaLabel="binding-site map"
              legend={mediaLegend().map((g) => ({ ...g, kind: "box" }))} />
      </div>
    {/if}
  </div>
{/if}
