<script>
  // Native thermodynamics of a dot [MX]_m[ML_n]_k0 from props/detachment.json (Z-type ladder) and
  // props/stability.json, in two views:
  //   view = "thermo":    stability in vacuum and the ladder (+ the report's solution-phase section)
  //   view = "synthesis": formation from MX and ML_n, post-synthesis shell (+ the multi-dot growth page)
  // ML_n is complexed by oleylamine. Vacuum (gas-phase) free energies of MACE-MH-1; formulas in model.js.
  import Plot from "./Plot.svelte";
  import TrimmedPage from "./TrimmedPage.svelte";
  import { OAM_COMPLEXES } from "./oamComplexes.js";
  import { KB, LN10, UNITS, buildModel, populations, dmuOf, amineLnZ, amineSpeciation, solveClosed,
           hasFormation, log10Sstar } from "./model.js";

  // reportUrl: the report page with the solution-phase section (thermo); growthUrl: synthesis.html (synthesis)
  let { detUrl = null, stabUrl = null, view = "thermo", reportUrl = null, growthUrl = null } = $props();

  let model = $state(null);
  let status = $state("idle");
  $effect(() => {
    const d = detUrl, s = stabUrl;
    if (!d) { model = null; status = "none"; return; }
    status = "loading";
    model = null;
    Promise.all([fetch(d).then((r) => { if (!r.ok) throw new Error(r.status); return r.json(); }),
                 s ? fetch(s).then((r) => (r.ok ? r.json() : null)).catch(() => null) : Promise.resolve(null)])
      .then(([det, stab]) => { if (d !== detUrl) return; model = buildModel(det, stab); status = "ready"; })
      .catch(() => { if (d === detUrl) status = "error"; });
  });

  // ---------------- shared settings ----------------
  let unit = $state("eV");
  let confMode = $state("both");          // "both" | "ideal" | "none"
  let delta = $state(0);                  // solution correction per bound OAm (eV)
  const C_NEAT = 3.04;                    // neat oleylamine (M)

  const iTof = (t) => Math.max(0, model.T.findIndex((v) => Math.abs(v - t) < 1e-6));
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
  const confs = () => (confMode === "both" ? [true, false] : [confMode === "ideal"]);
  const CONF = { true: { color: "#eb6834", label: "ideal S<sub>conf</sub> (upper bound)" },
                 false: { color: "#2a78d6", label: "S<sub>conf</sub> = 0 (lower bound)" } };

  let complexes = $derived(model ? OAM_COMPLEXES.complexes[model.unit] || null : null);
  const lnZ = (iT, cOAm) => amineLnZ(complexes, iT, model.T[iT], cOAm, delta);

  // ---------------- labels ----------------
  let dotHtml = $derived(model ? `[${html(model.core || "MX")}]<sub>${model.m ?? "m"}</sub>[${html(model.unit)}]<sub>${model.k0 ?? "k"}</sub>` : "");
  let unitHtml = $derived(model ? html(model.unit) : "");
  let kLabel = $derived(model ? `k, ${sub(model.unit)} units on the dot` : "");
  let kTicks = $derived.by(() => {
    if (!model) return [];
    const a = model.states[0].k, b = model.states.at(-1).k, step = b - a > 16 ? 2 : 1, out = [];
    for (let k = a; k <= b; k++) if ((k - model.k0) % step === 0) out.push({ v: k, label: String(k) });
    return out;
  });
  let kDom = $derived(model ? [model.states[0].k - 0.7, model.states.at(-1).k + 0.7] : [0, 1]);
  const atLeast = (k) => model && ((model.capDown && k <= model.states[0].k + 0.5) || (model.capUp && model.hasUp && k >= model.states.at(-1).k - 0.5));

  // =====================================================================================
  // 1. Ladder
  // =====================================================================================
  let lT = $state(300);
  let lLogC = $state(-6);
  let lI = $derived(model ? iTof(lT) : 0);
  let lDmu = $derived(model ? KB * model.T[lI] * lLogC * LN10 : 0);
  let lPop = $derived(model ? confs().map((c) => ({ c, ...populations(model, lI, lDmu, c) })) : []);

  let costLayers = $derived.by(() => {
    if (!model) return [];
    const t = model.T[lI], kT = KB * t, f = uf(), w = 0.62;
    const L = [];
    L.push({ type: "bars", bars: model.units.map((u) => ({ x: u.k - 0.5, y0: 0, y1: u.cost[lI] * f, w,
             color: CONF.false.color, opacity: confMode === "ideal" ? 0.25 : 0.85 })) });
    if (confMode !== "none" && model.hasConf)
      L.push({ type: "ticks", ticks: model.units.map((u) => ({ x: u.k - 0.5, y: (u.cost[lI] - kT * u.dS) * f, w: w + 0.12, color: CONF.true.color })) });
    L.push({ type: "points", pts: model.units.map((u) => [u.k - 0.5, u.dE * f]), color: "#64748b", hollow: true });
    L.push({ type: "hline", y: -lDmu * f, color: "#0f172a", label: `k_BT ln(1 M / c) = ${fe(-lDmu)} ${ul()}` });
    L.push({ type: "hline", y: 0, color: "#94a3b8", dash: "" });
    L.push({ type: "vline", x: model.k0, label: "as built" });
    return L;
  });
  let costDom = $derived.by(() => {
    if (!model) return [0, 1];
    const f = uf(), ys = model.units.flatMap((u) => [u.cost[lI], u.dE, u.cost[lI] - KB * model.T[lI] * u.dS]).concat([0, -lDmu]);
    const a = Math.min(...ys), b = Math.max(...ys), pad = 0.08 * (b - a || 1);
    return [(a - pad) * f, (b + pad) * f];
  });
  const costHover = (x) => {
    const u = model.units.find((q) => Math.abs(q.k - 0.5 - x) < 0.5);
    if (!u) return null;
    const kT = KB * model.T[lI];
    const which = u.side === "removal" ? `removal ${u.order}` : `addition ${u.order} (shown as its desorption)`;
    return { x: u.k - 0.5, title: `unit between k = ${u.k - 1} and ${u.k} · ${which}`,
             rows: [{ text: `ΔE ${fe(u.dE)} ${ul()}`, color: "#64748b", dash: "2 2" },
                    { text: `ΔG(${model.T[lI]} K) ${fe(u.cost[lI])} ${ul()}`, color: CONF.false.color },
                    ...(model.hasConf ? [{ text: `ΔG − TΔS<sub>conf</sub> ${fe(u.cost[lI] - kT * u.dS)} ${ul()}`, color: CONF.true.color }] : [])] };
  };

  let profLayers = $derived.by(() => {
    if (!model) return [];
    const f = uf();
    return lPop.map((P) => ({ type: "line", pts: model.states.map((s, j) => [s.k, (P.omega[j]) * f]), color: CONF[P.c].color }))
      .concat(lPop.map((P) => ({ type: "points", pts: model.states.map((s, j) => [s.k, P.omega[j] * f]), color: CONF[P.c].color, r: 2.5 })))
      .concat([{ type: "vline", x: model.k0, label: "as built" }]);
  });
  let profDom = $derived.by(() => {
    if (!lPop.length) return [0, 1];
    const ys = lPop.flatMap((P) => P.omega), a = Math.min(...ys), b = Math.min(Math.max(...ys), a + 10);
    return [(a - 0.05 * (b - a + 0.1)) * uf(), (b + 0.05 * (b - a + 0.1)) * uf()];
  });
  let popLayers = $derived.by(() => {
    if (!model) return [];
    const w = lPop.length > 1 ? 0.36 : 0.6;
    return lPop.map((P, i) => ({ type: "bars", bars: model.states.map((s, j) => ({
      x: s.k + (lPop.length > 1 ? (i ? 0.19 : -0.19) : 0), y0: 0, y1: P.p[j], w, color: CONF[P.c].color, opacity: 0.8 })) }))
      .concat([{ type: "vline", x: model.k0, label: "as built" }]);
  });
  const popHover = (x) => {
    const j = model.states.findIndex((s) => Math.abs(s.k - x) < 0.5);
    if (j < 0) return null;
    const s = model.states[j];
    return { x: s.k, title: `k = ${s.k}${atLeast(s.k) ? " (search limit)" : ""}`,
             rows: lPop.map((P) => ({ color: CONF[P.c].color, text: `P = ${P.p[j].toFixed(3)} · Ω − Ω<sub>k0</sub> = ${fe(P.omega[j])} ${ul()}` })) };
  };

  // <k> map over T and log10 c_free
  const MAP_LOGC = Array.from({ length: 61 }, (_, i) => -14 + 0.25 * i);
  let mapConf = $derived(confMode !== "none");
  let mapCells = $derived.by(() => {
    if (!model) return [];
    const cells = [], T = model.T, k0 = model.k0;
    const spanDn = k0 - model.states[0].k, spanUp = model.states.at(-1).k - k0;
    for (let i = 0; i < T.length; i++) {
      const t0 = i ? (T[i - 1] + T[i]) / 2 : T[0] - 12.5, t1 = i < T.length - 1 ? (T[i] + T[i + 1]) / 2 : T[i] + 12.5;
      for (const l of MAP_LOGC) {
        const km = populations(model, i, KB * T[i] * l * LN10, mapConf).kmean;
        cells.push({ x0: t0, x1: t1, y0: l - 0.125, y1: l + 0.125, k: km, color: kColor(km - k0, spanDn, spanUp) });
      }
    }
    return cells;
  });
  function kColor(dk, dn, upSpan) {
    const a = dk < 0 ? Math.min(1, -dk / Math.max(dn, 1)) : Math.min(1, dk / Math.max(upSpan, 1));
    const g = [241, 240, 236], to = dk < 0 ? [214, 69, 65] : [42, 120, 214], e = Math.pow(a, 0.75);
    return `rgb(${g.map((v, j) => Math.round(v + (to[j] - v) * e)).join(",")})`;
  }
  const mapHover = (x) => {
    const i = model.T.reduce((b, t, j) => (Math.abs(t - x) < Math.abs(model.T[b] - x) ? j : b), 0);
    const rows = [-12, -9, -6, -3, 0].map((l) => {
      const km = populations(model, i, KB * model.T[i] * l * LN10, mapConf).kmean;
      return { text: `c = ${cTxt(10 ** l)}: ⟨k⟩ = ${fk(km)}${atLeast(Math.round(km)) ? " (≥ limit)" : ""}` };
    });
    return { x: model.T[i], title: `T = ${model.T[i]} K`, rows };
  };

  // site map: first removal per symmetry-distinct unit
  let siteLayers = $derived.by(() => {
    if (!model?.siteClasses.length) return [];
    const f = uf(), sc = model.siteClasses.slice().sort((a, b) => a.dE - b.dE);
    return [{ type: "bars", bars: sc.map((c, i) => ({ x: i + 1, y0: 0, y1: c.dE * f, w: 0.7, color: c.molecular === false ? "#94a3b8" : "#475569" })) },
            { type: "hline", y: 0, color: "#94a3b8", dash: "" }];
  });
  let siteSorted = $derived(model ? model.siteClasses.slice().sort((a, b) => a.dE - b.dE) : []);
  const siteHover = (x) => {
    const i = Math.round(x) - 1, c = siteSorted[i];
    if (!c) return null;
    return { x: i + 1, title: `site class ${i + 1}`, rows: [{ text: `ΔE ${fe(c.dE)} ${ul()} · ×${c.mult} equivalent${c.site ? " · " + c.site : ""}${c.molecular === false ? " · non-molecular unit" : ""}` }] };
  };

  // =====================================================================================
  // 2. Formation
  // =====================================================================================
  let fT = $state(500);
  let fLogCL = $state(-1);         // total ML_n precursor
  let fLogCX = $state(-1.3);       // X precursor
  let fXi = $state(0.5);           // conversion of the X precursor
  let fOAm = $state(C_NEAT);
  const C_DOT_F = 1e-6;
  let fI = $derived(model ? iTof(fT) : 0);
  let fLeft = $derived(10 ** fLogCL - fXi * 10 ** fLogCX);           // ML_n left in solution (M)
  const fFree = (i, cOAm = fOAm) => Math.max(fLeft, 1e-30) / Math.exp(lnZ(i, cOAm));
  let fRows = $derived.by(() => {
    if (!model || !hasFormation(model)) return [];
    return confs().map((c) => ({ c, S: log10Sstar(model, fI, fFree(fI), C_DOT_F, c),
                                  k: populations(model, fI, dmuOf(fFree(fI), model.T[fI]), c).kmean }));
  });
  let sStarLayers = $derived.by(() => {
    if (!model || !hasFormation(model) || fLeft <= 0) return [];
    const L = [];
    for (const [cOAm, dash] of [[fOAm, ""], ...(fOAm > 0 && complexes ? [[0, "6 4"]] : [])])
      for (const c of confs())
        L.push({ type: "line", color: CONF[c].color, dash, pts: model.T.map((t, i) => [t, t >= F_DOM[0] ? log10Sstar(model, i, fFree(i, cOAm), C_DOT_F, c) : NaN]) });
    L.push({ type: "vline", x: model.T[fI], label: `${model.T[fI]} K` });
    return L;
  });
  let sStarDom = $derived.by(() => {
    const ys = sStarLayers.filter((l) => l.pts).flatMap((l) => l.pts.map((p) => p[1])).filter(isFinite);
    if (!ys.length) return [0, 1];
    const a = Math.min(...ys), b = Math.max(...ys);
    return [Math.floor(a - 0.2), Math.ceil(b + 0.2)];
  });
  const sStarHover = (x) => {
    const i = model.T.reduce((b, t, j) => (Math.abs(t - x) < Math.abs(model.T[b] - x) ? j : b), 0);
    return { x: model.T[i], title: `T = ${model.T[i]} K · free ${unitHtml} ${cTxt(fFree(i))}`,
             rows: confs().map((c) => ({ color: CONF[c].color, text: `log<sub>10</sub> S* = ${log10Sstar(model, i, fFree(i), C_DOT_F, c).toFixed(2)}` })) };
  };
  let fShellLayers = $derived.by(() => {
    if (!model || fLeft <= 0) return [];
    const L = [];
    for (const [cOAm, dash] of [[fOAm, ""], ...(fOAm > 0 && complexes ? [[0, "6 4"]] : [])])
      for (const c of confs())
        L.push({ type: "line", color: CONF[c].color, dash, pts: model.T.map((t, i) => [t, t >= F_DOM[0] ? populations(model, i, dmuOf(fFree(i, cOAm), t), c).kmean : NaN]) });
    L.push({ type: "hline", y: model.k0, color: "#334155", label: `as built, k = ${model.k0}` });
    L.push({ type: "vline", x: model.T[fI], label: `${model.T[fI]} K` });
    return L;
  });
  const fShellHover = (x) => {
    const i = model.T.reduce((b, t, j) => (Math.abs(t - x) < Math.abs(model.T[b] - x) ? j : b), 0);
    return { x: model.T[i], title: `T = ${model.T[i]} K`,
             rows: confs().map((c) => ({ color: CONF[c].color, text: `⟨k⟩ = ${fk(populations(model, i, dmuOf(fFree(i), model.T[i]), c).kmean)}` })) };
  };
  // stability in vacuum (props/stability.json): per MX against bulk MX + ML_n, per unit against the monomers
  const oneLine = (vals, color) => (vals ? [{ type: "line", color, pts: model.T.map((t, i) => [t, vals[i] * uf()]) },
                                           { type: "hline", y: 0, color: "#94a3b8", dash: "" }] : []);
  let decLayers = $derived(model?.decomposition ? oneLine(model.decomposition.map((v) => v / model.m), "#475569") : []);
  let bindLayers = $derived(model?.bindPerUnit ? oneLine(model.bindPerUnit, "#475569") : []);
  const domOf = (layers) => {
    const ys = layers.filter((l) => l.pts).flatMap((l) => l.pts.map((p) => p[1])).filter(isFinite).concat([0]);
    const a = Math.min(...ys), b = Math.max(...ys), p = 0.08 * (b - a || 1);
    return [a - p, b + p];
  };
  const valueHover = (layers, label) => (x) => {
    const l = layers.find((q) => q.type === "line");
    if (!l) return null;
    const p = l.pts.reduce((b, q) => (Math.abs(q[0] - x) < Math.abs(b[0] - x) ? q : b), l.pts[0]);
    return { x: p[0], title: `T = ${p[0]} K`, rows: [{ color: l.color, text: `${label} ${fe(p[1] / uf())} ${ul()}` }] };
  };

  // =====================================================================================
  // 3. Post-synthesis (closed system)
  // =====================================================================================
  let pT = $state(300);
  let pLogDot = $state(-6);
  let pOAm = $state(0);
  let pLogAdd = $state(-12);       // -12 = nothing added
  let pWash = $state(0.9);
  let pI = $derived(model ? iTof(pT) : 0);
  let pAdd = $derived(pLogAdd <= -11.99 ? 0 : 10 ** pLogAdd);
  const closed = (i, cdot, cOAm, c, add = pAdd, kRef = null) => solveClosed(model, i, cdot, add, lnZ(i, cOAm), c, kRef);
  let pRows = $derived(model ? confs().map((c) => ({ c, ...closed(pI, 10 ** pLogDot, pOAm, c) })) : []);
  let firstUnit = $derived.by(() => {
    if (!model) return null;
    const u = model.units.find((q) => q.k === model.k0);
    if (!u) return null;
    const t = model.T[pI], kT = KB * t;
    const thr = kT * Math.log(1 / 10 ** pLogDot) + kT * lnZ(pI, pOAm);
    return { dG: u.cost[pI], dGc: u.cost[pI] - kT * u.dS, thr };
  });
  const LOGDOT = Array.from({ length: 25 }, (_, i) => -9 + 0.25 * i);
  let dilLayers = $derived.by(() => {
    if (!model) return [];
    const L = [];
    for (const [cOAm, dash] of [[pOAm, ""], ...(pOAm > 0 ? [[0, "6 4"]] : complexes ? [[C_NEAT, "2 3"]] : [])])
      for (const c of confs())
        L.push({ type: "line", color: CONF[c].color, dash, pts: LOGDOT.map((l) => [l, closed(pI, 10 ** l, cOAm, c).kmean]) });
    L.push({ type: "hline", y: model.k0, color: "#334155", label: `as built, k = ${model.k0}` });
    L.push({ type: "vline", x: pLogDot });
    return L;
  });
  let dilAlt = $derived(pOAm > 0 ? "no OAm" : complexes ? "neat OAm (3.04 M)" : "");
  const OAMS = Array.from({ length: 31 }, (_, i) => (C_NEAT * i) / 30);
  let oamLayers = $derived.by(() => {
    if (!model || !complexes) return [];
    const L = confs().map((c) => ({ type: "line", color: CONF[c].color, pts: OAMS.map((a) => [a, closed(pI, 10 ** pLogDot, a, c).kmean]) }));
    L.push({ type: "hline", y: model.k0, color: "#334155", label: `as built, k = ${model.k0}` });
    L.push({ type: "vline", x: pOAm });
    return L;
  });
  // washing: each wash removes a fraction of the dissolved ML_n (free and complexed); the dots re-equilibrate
  let washLayers = $derived.by(() => {
    if (!model) return [];
    const cdot = 10 ** pLogDot, L = [];
    for (const c of confs()) {
      let r = closed(pI, cdot, pOAm, c), pts = [[0, r.kmean]];
      for (let n = 1; n <= 8; n++) {
        const left = (1 - pWash) * Math.exp(lnZ(pI, pOAm)) * r.cFree;
        r = closed(pI, cdot, pOAm, c, left, r.kmean);
        pts.push([n, r.kmean]);
      }
      L.push({ type: "line", color: CONF[c].color, pts: pts.flatMap(([n, k], j) => (j ? [[n - 0.5, pts[j - 1][1]], [n - 0.5, k], [n, k]] : [[0, k]])) });
      L.push({ type: "points", color: CONF[c].color, pts, r: 3 });
    }
    L.push({ type: "hline", y: model.k0, color: "#334155", label: `as built, k = ${model.k0}` });
    return L;
  });
  let annealLayers = $derived.by(() => {
    if (!model) return [];
    const L = confs().map((c) => ({ type: "line", color: CONF[c].color, pts: model.T.map((t, i) => [t, closed(i, 10 ** pLogDot, pOAm, c).kmean]) }));
    L.push({ type: "hline", y: model.k0, color: "#334155", label: `as built, k = ${model.k0}` });
    L.push({ type: "vline", x: model.T[pI], label: `${model.T[pI]} K` });
    return L;
  });
  const lineHover = (layers, fmtX) => (x) => {
    const lines = layers.filter((l) => l.type === "line");
    if (!lines.length) return null;
    const rows = lines.map((l) => {
      const p = l.pts.reduce((b, q) => (Math.abs(q[0] - x) < Math.abs(b[0] - x) ? q : b), l.pts[0]);
      return { color: l.color, dash: l.dash, text: `⟨k⟩ = ${fk(p[1])}`, x: p[0] };
    });
    return { x: rows[0].x, title: fmtX(rows[0].x), rows };
  };
  let speciesF = $derived(model && complexes ? amineSpeciation(complexes, fI, model.T[fI], fOAm, delta) : null);
  let speciesP = $derived(model && complexes ? amineSpeciation(complexes, pI, model.T[pI], pOAm, delta) : null);

  const tempTicks = [100, 200, 300, 400, 500, 600, 700, 800].map((v) => ({ v, label: String(v) }));
  const F_DOM = [275, 800];                       // formation: synthesis temperatures
  const fTicks = [300, 400, 500, 600, 700, 800].map((v) => ({ v, label: String(v) }));
  // <k> axis fitted to the curves, at least k0 +- 2, within the ladder
  function kDomOf(layers) {
    const ys = layers.filter((l) => l.type === "line").flatMap((l) => l.pts.map((p) => p[1])).filter(isFinite).concat([model.k0 - 2, model.k0 + 2]);
    const lo = Math.max(model.states[0].k, Math.floor(Math.min(...ys))) - 0.5, hi = Math.min(model.states.at(-1).k, Math.ceil(Math.max(...ys))) + 0.5;
    return [lo, hi];
  }
  const sliderCls = "w-40 accent-brand-600";
  const seg = (on) => `px-3 py-1.5 rounded-lg font-bold transition-all ${on ? "bg-white text-slate-900 shadow-sm" : "text-slate-500 hover:text-slate-800"}`;
</script>

{#snippet slider(label, value, min, max, step, set, text)}
  <label class="flex items-center gap-2 text-xs text-slate-600">
    <span class="font-bold text-slate-500 whitespace-nowrap">{@html label}</span>
    <input type="range" {min} {max} {step} value={value} oninput={(e) => set(+e.currentTarget.value)} class={sliderCls} />
    <span class="font-mono text-slate-900 min-w-[5.5rem]">{@html text}</span>
  </label>
{/snippet}

{#snippet badge(text, tone = "slate")}
  <span class="px-2 py-0.5 rounded border text-[10px] font-bold {tone === 'amber' ? 'bg-amber-50 text-amber-800 border-amber-200' : 'bg-slate-50 text-slate-600 border-slate-200'}">{@html text}</span>
{/snippet}

{#snippet section(title, sub)}
  <div class="flex items-baseline gap-3 flex-wrap border-b border-slate-200 pb-2 mt-2">
    <h3 class="text-base font-bold text-slate-900">{@html title}</h3>
    <span class="text-xs text-slate-500">{@html sub}</span>
  </div>
{/snippet}

{#snippet panel(title, caption)}
  <h4 class="text-xs font-bold text-slate-500 uppercase tracking-widest mb-1">{@html title}</h4>
  {#if caption}<p class="text-[11px] text-slate-500 mb-1 max-w-4xl leading-relaxed">{@html caption}</p>{/if}
{/snippet}

{#if status === "loading"}
  <div class="flex items-center justify-center h-full text-brand-600 text-sm font-bold">
    <div class="w-6 h-6 border-4 border-brand-500 border-t-transparent rounded-full animate-spin mr-3"></div>Loading the ligand ladder…
  </div>
{:else if status === "error"}
  <div class="flex items-center justify-center h-full text-slate-500 text-sm">Could not load {detUrl}.</div>
{:else}
  <div class="flex flex-col gap-6 pb-6">
    {#if model}
      <!-- Header and shared settings -->
      <div class="flex flex-col gap-2">
        <div class="text-sm text-slate-700">
          <span class="font-bold text-slate-900">{@html dotHtml}</span>
          · Z-type unit {@html unitHtml} ·
          ladder from k = {model.states[0].k} to {model.states.at(-1).k}
          ({model.nDown} removals{model.hasUp ? `, ${model.nUp} additions` : ", no addition site"})
        </div>
        <div class="flex flex-wrap gap-2">
          {@render badge("vacuum (gas phase), 1 M standard state", "amber")}
          {@render badge("solution-phase ladder pending")}
          {#if !model.hasConf}{@render badge("no configurational entropy in this ladder", "amber")}{/if}
        </div>
        <div class="flex flex-wrap items-center gap-4 text-xs mt-1">
          <div class="flex gap-1 bg-slate-100 p-1 rounded-xl border border-slate-200">
            {#each Object.entries(UNITS) as [k, u]}<button class={seg(unit === k)} onclick={() => (unit = k)}>{u[1]}</button>{/each}
          </div>
          {#if model.hasConf}
            <span class="font-bold text-slate-400">S<sub>conf</sub></span>
            <div class="flex gap-1 bg-slate-100 p-1 rounded-xl border border-slate-200">
              {#each [["both", "both bounds"], ["ideal", "ideal count"], ["none", "zero"]] as [k, l]}<button class={seg(confMode === k)} onclick={() => (confMode = k)}>{l}</button>{/each}
            </div>
          {/if}
          {#if view === "synthesis" && complexes}
            {@render slider("OAm solution correction", delta, -0.4, 0.2, 0.01, (v) => (delta = v), `${fe(delta)} ${ul()} per OAm`)}
          {/if}
        </div>
        <p class="text-[11px] text-slate-500 leading-relaxed max-w-4xl">
          k counts the {@html unitHtml} units on the dot; the builder made k = {model.k0}.
          {#if model.capDown || model.capUp}The search stops after 12 steps each way, so the end states mean "at least 12 removed / added".{/if}
          {#if view === "synthesis"}
            {#if complexes}Oleylamine binds free {@html unitHtml} as {@html complexes.map((c) => html(c.label)).join(", ")}
              (MACE-MH-1 omat_pbe, vacuum). The correction adds a constant per bound OAm, for what vacuum misses in solution
              (it overestimates the entropy lost by the amine; negative values strengthen binding).
            {:else}No oleylamine complexes computed for {@html unitHtml}: the amine is left out.{/if}
          {/if}
        </p>
      </div>
    {/if}

    {#if view === "thermo" && model}
      <!-- Stability in vacuum -->
      {#if model.decomposition || model.bindPerUnit}
        {@render section("Stability in vacuum", `MACE-MH-1, ideal gases at 1 M, bulk ${html(model.core || "MX")} with phonons`)}
        {#if model.decomposition}
          <div>
            {@render panel(`Decomposition free energy, per ${html(model.core || "MX")}`,
              `ΔG<sub>dec</sub> = [G(dot) − ${model.m} G(${html(model.core)}, bulk) − ${model.k0} G°(${unitHtml})] / ${model.m}: positive, the dot is less stable than bulk ${html(model.core)} plus free ${unitHtml} (its excess surface free energy).`)}
            <Plot layers={decLayers} xDomain={[50, 800]} yDomain={domOf(decLayers)} xTicks={tempTicks} xLabel="T (K)" yLabel={`ΔG_dec per ${sub(model.core)} (${ul()})`}
                  width={1000} height={300} hover={valueHover(decLayers, "ΔG<sub>dec</sub>")} ariaLabel="decomposition free energy" />
          </div>
        {/if}
        {#if model.bindPerUnit}
          <div>
            {@render panel("Binding free energy, per unit",
              `ΔG<sub>bind</sub> = [G(dot) − ${model.m} G°(${html(model.core)}) − ${model.k0} G°(${unitHtml})] / ${model.m + model.k0}, against the free molecules at 1 M: negative, the dot holds together.`)}
            <Plot layers={bindLayers} xDomain={[50, 800]} yDomain={domOf(bindLayers)} xTicks={tempTicks} xLabel="T (K)" yLabel={`ΔG_bind per unit (${ul()})`}
                  width={1000} height={300} hover={valueHover(bindLayers, "ΔG<sub>bind</sub>")} ariaLabel="binding free energy" />
          </div>
        {/if}
      {/if}

      <!-- Ladder -->
      {@render section(`${unitHtml} ladder`, `cost of every unit, and the shell against a reservoir of free ${unitHtml}`)}
      <div class="flex flex-wrap gap-x-6 gap-y-2">
        {@render slider("T", lT, 50, 800, 25, (v) => (lT = v), `${lT} K`)}
        {@render slider(`free ${unitHtml}`, lLogC, -14, 1, 0.25, (v) => (lLogC = v), cTxt(10 ** lLogC))}
      </div>
      <div class="flex flex-wrap gap-4 text-xs text-slate-700">
        <span>Δμ = k<sub>B</sub>T ln(c / 1 M) = <b class="font-mono">{fe(lDmu)} {ul()}</b></span>
        {#each lPop as P}<span class="flex items-center gap-1.5"><span class="w-3 h-1 inline-block" style="background:{CONF[P.c].color}"></span>⟨k⟩ = <b class="font-mono">{fk(P.kmean)}</b>{atLeast(Math.round(P.kmean)) ? " (at the search limit)" : ""}</span>{/each}
      </div>
      <div>
        {@render panel("Cost of each unit (desorption, positive = bound)",
          `Bar between k − 1 and k: the free energy to take one ${unitHtml} off a dot with k units. Left of "as built": the units the dot has; right: units added (their binding read as desorption). Units whose bar lies below the dashed line leave at this concentration.`)}
        <Plot layers={costLayers} xDomain={kDom} yDomain={costDom} xTicks={kTicks} xLabel={kLabel} yLabel={`ΔG to remove the unit (${ul()})`}
              width={1000} height={360} hover={costHover} ariaLabel="cost of each unit"
              legend={[{ kind: "box", color: CONF.false.color, label: `ΔG(${model.T[lI]} K)` },
                       ...(model.hasConf && confMode !== "none" ? [{ kind: "tick", color: CONF.true.color, label: "ΔG − TΔS<sub>conf</sub> (ideal)" }] : []),
                       { kind: "circle", color: "#64748b", label: "ΔE" },
                       { kind: "line", color: "#0f172a", dash: "5 4", label: "k<sub>B</sub>T ln(1 M / c): units below this line leave" }]} />
      </div>
      <div>
        {@render panel("Free-energy profile at this T and concentration",
          `Ω<sub>k</sub> − Ω<sub>k0</sub> = A<sub>k</sub> − (k − k<sub>0</sub>) Δμ − T (S<sub>k</sub> − S<sub>k0</sub>): the lowest state is the most populated.`)}
        <Plot layers={profLayers} xDomain={kDom} yDomain={profDom} xTicks={kTicks} xLabel={kLabel} yLabel={`Ω_k − Ω_k0 (${ul()})`}
              width={1000} height={300} hover={popHover} ariaLabel="free energy profile"
              legend={lPop.map((P) => ({ color: CONF[P.c].color, label: CONF[P.c].label }))} />
      </div>
      <div>
        {@render panel("Populations of the shell states", "")}
        <Plot layers={popLayers} xDomain={kDom} yDomain={[0, 1.05]} xTicks={kTicks} xLabel={kLabel} yLabel="population P_k"
              width={1000} height={240} hover={popHover} ariaLabel="populations"
              legend={lPop.map((P) => ({ color: CONF[P.c].color, label: CONF[P.c].label }))} />
      </div>
      <div>
        {@render panel(`Mean shell ⟨k⟩ against T and free ${unitHtml} (${mapConf ? "ideal S<sub>conf</sub>" : "S<sub>conf</sub> = 0"})`,
          "Dashed cross: the T and concentration set above. Hover for ⟨k⟩ at several concentrations.")}
        <Plot layers={[{ type: "cells", cells: mapCells }, { type: "vline", x: model.T[lI], color: "#0f172a", dash: "3 3" }, { type: "hline", y: lLogC, color: "#0f172a", dash: "3 3" }]}
              xDomain={[37.5, 812.5]} yDomain={[-14.125, 1.125]} xTicks={tempTicks}
              yTicks={[-14, -12, -10, -8, -6, -4, -2, 0].map((v) => ({ v, label: v === 0 ? "1 M" : pow10(v) }))}
              xLabel="T (K)" yLabel={`free ${sub(model.unit)} (M)`} width={1000} height={420} hover={mapHover} ariaLabel="mean shell map" />
        <div class="flex flex-wrap items-center gap-3 text-[11px] text-slate-500 mt-1 px-1">
          {#each [model.states[0].k, Math.round((model.states[0].k + model.k0) / 2), model.k0, ...(model.hasUp ? [Math.round((model.k0 + model.states.at(-1).k) / 2), model.states.at(-1).k] : [])] as k}
            <span class="flex items-center gap-1"><span class="w-3 h-3 rounded-sm inline-block border border-slate-200" style="background:{kColor(k - model.k0, model.k0 - model.states[0].k, model.states.at(-1).k - model.k0)}"></span>k = {k}{atLeast(k) ? " or beyond" : ""}</span>
          {/each}
          <span>red: units lost · blue: units gained</span>
        </div>
      </div>
      {#if model.siteClasses.length}
        <div>
          {@render panel("Site map: first removal, per symmetry-distinct unit", "")}
          <Plot layers={siteLayers} xDomain={[0.3, model.siteClasses.length + 0.7]} yDomain={(() => { const v = siteSorted.map((c) => c.dE * uf()).concat([0]); const a = Math.min(...v), b = Math.max(...v); return [a - 0.08 * (b - a), b + 0.08 * (b - a)]; })()}
                xTicks={siteSorted.length > 40 ? null : siteSorted.map((_, i) => ({ v: i + 1, label: String(i + 1) }))}
                xLabel="symmetry-distinct unit, cheapest first" yLabel={`ΔE to remove (${ul()})`} width={1000} height={300} hover={siteHover} ariaLabel="site map"
                legend={[{ kind: "box", color: "#475569", label: `molecular unit (M with ${model.n} own ${model.L})` }, { kind: "box", color: "#94a3b8", label: `non-molecular (completed with a neighbouring ${model.L})` }]} />
        </div>
      {/if}
    {/if}

    {#if view === "thermo" && reportUrl}
      {@render section("In solution", "the report's solution-phase model: Generalized Born solvation, removal ladder only, no configurational entropy; it moves here natively once the model includes both")}
      <TrimmedPage src={reportUrl} keep={(t) => t.startsWith("Stability & ligands in solution")} header={false} autosize title="Stability in solution" />
    {/if}

    {#if view === "synthesis" && model}
      <!-- Formation -->
      {@render section("Formation", `the dot as built, from MX monomers and ${unitHtml}, in oleylamine`)}
      <p class="text-[11px] text-slate-500 leading-relaxed max-w-4xl">
        {@html unitHtml} + X precursor → {@html html(model.core || "MX")} (conversion ξ, irreversible: sets the monomer supply);
        m {@html html(model.core || "MX")} + k {@html unitHtml} ⇌ {@html dotHtml}; {@html unitHtml} + p OAm ⇌ {@html unitHtml}(OAm)<sub>p</sub>.
        The monomer enters as its supersaturation against bulk {@html html(model.core || "MX")}, S = c / c<sub>sat</sub>, so the free energy of the
        molecular monomer cancels. S* is the supersaturation at which this exact dot is in equilibrium with the monomers (dots at 1 µM).
        No claim about which dot a synthesis actually makes.
      </p>
      {#if hasFormation(model)}
        <div class="flex flex-wrap gap-x-6 gap-y-2">
          {@render slider("T", fT, 300, 650, 25, (v) => (fT = v), `${fT} K`)}
          {@render slider(`total ${unitHtml}`, fLogCL, -4, 0.3, 0.05, (v) => (fLogCL = v), cTxt(10 ** fLogCL))}
          {@render slider(`X precursor`, fLogCX, -4, 0.3, 0.05, (v) => (fLogCX = v), cTxt(10 ** fLogCX))}
          {@render slider("conversion ξ", fXi, 0, 1, 0.01, (v) => (fXi = v), fXi.toFixed(2))}
          {#if complexes}{@render slider("OAm", fOAm, 0, C_NEAT, 0.02, (v) => (fOAm = v), fOAm >= C_NEAT - 1e-9 ? "neat (3.04 M)" : `${fOAm.toFixed(2)} M`)}{/if}
        </div>
        {#if fLeft <= 0}
          <div class="text-xs text-amber-800 bg-amber-50 border border-amber-200 rounded-lg px-3 py-2">The conversion uses up all the {@html unitHtml}: ξ·C<sub>X</sub> ≥ C<sub>{@html unitHtml}</sub>. Lower ξ or the X precursor.</div>
        {:else}
          <div class="flex flex-wrap gap-4 text-xs text-slate-700">
            <span>{@html unitHtml} left: <b class="font-mono">{@html cTxt(fLeft)}</b>, free: <b class="font-mono">{@html cTxt(fFree(fI))}</b></span>
            {#if speciesF}<span class="text-slate-500">({@html speciesF.species.filter((s) => s.frac > 0.005).map((s) => `${s.label === "free" ? "free" : html(s.label)} ${(100 * s.frac).toFixed(0)}%`).join(", ")})</span>{/if}
            {#each fRows as r}<span class="flex items-center gap-1.5"><span class="w-3 h-1 inline-block" style="background:{CONF[r.c].color}"></span>log<sub>10</sub> S* = <b class="font-mono">{r.S.toFixed(2)}</b>, ⟨k⟩ = <b class="font-mono">{fk(r.k)}</b></span>{/each}
          </div>
          <div>
            {@render panel("Critical supersaturation of the dot as built", "Above S*, monomers assemble into this dot; below, it dissolves. Dashed: without oleylamine.")}
            <Plot layers={sStarLayers} xDomain={F_DOM} yDomain={sStarDom} xTicks={fTicks} xLabel="T (K)" yLabel="log₁₀ S*" width={1000} height={320}
                  hover={sStarHover} ariaLabel="critical supersaturation"
                  legend={[...confs().map((c) => ({ color: CONF[c].color, label: CONF[c].label })), ...(fOAm > 0 && complexes ? [{ color: "#64748b", dash: "6 4", label: "without OAm" }] : [])]} />
          </div>
          <div>
            {@render panel("Shell at these conditions (open system)", `⟨k⟩ with the free ${unitHtml} left after conversion and complexation: does the dot keep the shell it was built with?`)}
            <Plot layers={fShellLayers} xDomain={F_DOM} yDomain={kDomOf(fShellLayers)} xTicks={fTicks} xLabel="T (K)" yLabel="⟨k⟩" width={1000} height={320}
                  hover={fShellHover} ariaLabel="shell during formation"
                  legend={[...confs().map((c) => ({ color: CONF[c].color, label: CONF[c].label })), ...(fOAm > 0 && complexes ? [{ color: "#64748b", dash: "6 4", label: "without OAm" }] : [])]} />
          </div>
        {/if}
      {:else}
        <div class="text-xs text-slate-500">No props/stability.json for this dot: the formation free energies are missing.</div>
      {/if}

      <!-- Post-synthesis -->
      {@render section("Post-synthesis", `dots as built, washed, diluted, stored or annealed; every free ${unitHtml} comes from the dots (or is added)`)}
      <div class="flex flex-wrap gap-x-6 gap-y-2">
        {@render slider("T", pT, 50, 800, 25, (v) => (pT = v), `${pT} K`)}
        {@render slider("dots", pLogDot, -9, -3, 0.25, (v) => (pLogDot = v), cTxt(10 ** pLogDot))}
        {#if complexes}{@render slider("OAm", pOAm, 0, C_NEAT, 0.02, (v) => (pOAm = v), pOAm >= C_NEAT - 1e-9 ? "neat (3.04 M)" : pOAm === 0 ? "none" : `${pOAm.toFixed(2)} M`)}{/if}
        {@render slider(`added ${unitHtml}`, pLogAdd, -12, -1, 0.25, (v) => (pLogAdd = v), pAdd ? cTxt(pAdd) : "none")}
        {@render slider("wash removes", pWash, 0.5, 0.999, 0.001, (v) => (pWash = v), `${(100 * pWash).toFixed(1)}%`)}
      </div>
      <div class="flex flex-wrap gap-4 text-xs text-slate-700">
        {#each pRows as r}<span class="flex items-center gap-1.5"><span class="w-3 h-1 inline-block" style="background:{CONF[r.c].color}"></span>⟨k⟩ = <b class="font-mono">{fk(r.kmean)}</b>, free {@html unitHtml} <b class="font-mono">{@html cTxt(r.cFree)}</b></span>{/each}
        {#if speciesP && pOAm > 0}<span class="text-slate-500">(dissolved {@html unitHtml}: {@html speciesP.species.filter((s) => s.frac > 0.005).map((s) => `${s.label === "free" ? "free" : html(s.label)} ${(100 * s.frac).toFixed(0)}%`).join(", ")})</span>{/if}
      </div>
      {#if firstUnit}
        {@const lose0 = firstUnit.dG < firstUnit.thr}
        {@const lose1 = firstUnit.dGc < firstUnit.thr}
        <div class="text-xs text-slate-700 bg-slate-50 border border-slate-200 rounded-lg px-3 py-2 max-w-4xl leading-relaxed">
          <b>First unit</b> (nothing added): the dot loses units if the cost of the first one is below
          k<sub>B</sub>T ln(1 M / c<sub>dot</sub>){#if pOAm > 0 && complexes} + k<sub>B</sub>T ln(1 + Σ<sub>p</sub> K<sub>p</sub> c<sub>OAm</sub><sup>p</sup>){/if} = <b class="font-mono">{fe(firstUnit.thr)} {ul()}</b>.
          Cost ΔG = <b class="font-mono">{fe(firstUnit.dG)}</b>{model.hasConf ? "" : ` ${ul()}`}{#if model.hasConf}, ΔG − TΔS<sub>conf</sub> = <b class="font-mono">{fe(firstUnit.dGc)} {ul()}</b>{/if}:
          <span class="font-bold {lose0 || lose1 ? 'text-rose-700' : 'text-emerald-700'}">
            {lose0 && lose1 ? "the first unit leaves" : !lose0 && !lose1 ? "the first unit stays" : "borderline: it depends on the configurational entropy"}</span>.
          The first unit alone is a quick test; the curves below solve the whole ladder, where cheaper later units can leave together.
        </div>
      {/if}
      <div>
        {@render panel("Dilution", `⟨k⟩ against the dot concentration: the units that leave stay in solution, free or bound to oleylamine.`)}
        <Plot layers={dilLayers} xDomain={[-9, -3]} yDomain={kDomOf(dilLayers)}
              xTicks={[-9, -8, -7, -6, -5, -4, -3].map((v) => ({ v, label: v === -6 ? "1 µM" : v === -3 ? "1 mM" : v === -9 ? "1 nM" : pow10(v) }))}
              xLabel="dot concentration (M)" yLabel="⟨k⟩" width={1000} height={320} ariaLabel="dilution"
              hover={lineHover(dilLayers, (x) => `dots ${cTxt(10 ** x)}`)}
              legend={[...confs().map((c) => ({ color: CONF[c].color, label: CONF[c].label })), ...(dilAlt ? [{ color: "#64748b", dash: pOAm > 0 ? "6 4" : "2 3", label: dilAlt }] : [])]} />
      </div>
      {#if complexes}
        <div>
          {@render panel("Oleylamine", `⟨k⟩ against the amine concentration at this T and dot concentration: ligand stripping by complexation.`)}
          <Plot layers={oamLayers} xDomain={[0, C_NEAT]} yDomain={kDomOf(oamLayers)} xLabel="oleylamine (M), 3.04 M = neat" yLabel="⟨k⟩" width={1000} height={320}
                ariaLabel="oleylamine" hover={lineHover(oamLayers, (x) => `OAm ${x.toFixed(2)} M`)}
                legend={confs().map((c) => ({ color: CONF[c].color, label: CONF[c].label }))} />
        </div>
      {/if}
      <div>
        {@render panel("Washing cycles", `Each wash removes ${(100 * pWash).toFixed(1)}% of the dissolved ${unitHtml} (free and complexed); the dots, at the same concentration and in the same medium, re-equilibrate.`)}
        <Plot layers={washLayers} xDomain={[-0.4, 8.4]} yDomain={kDomOf(washLayers)} xTicks={[0, 1, 2, 3, 4, 5, 6, 7, 8].map((v) => ({ v, label: v ? String(v) : "start" }))}
              xLabel="wash" yLabel="⟨k⟩" width={1000} height={320} ariaLabel="washing"
              hover={lineHover(washLayers, (x) => (Math.round(x) ? `after wash ${Math.round(x)}` : "start"))}
              legend={confs().map((c) => ({ color: CONF[c].color, label: CONF[c].label }))} />
      </div>
      <div>
        {@render panel("Annealing (closed system)", "⟨k⟩ against T with nothing added: desorption is self-limiting, the released units raise the free concentration.")}
        <Plot layers={annealLayers} xDomain={[50, 800]} yDomain={kDomOf(annealLayers)} xTicks={tempTicks} xLabel="T (K)" yLabel="⟨k⟩" width={1000} height={320}
              ariaLabel="annealing" hover={lineHover(annealLayers, (x) => `T = ${Math.round(x)} K`)}
              legend={confs().map((c) => ({ color: CONF[c].color, label: CONF[c].label }))} />
      </div>
      <p class="text-[11px] text-slate-500 leading-relaxed max-w-4xl">
        Mass balance: (k<sub>0</sub> − ⟨k⟩) c<sub>dot</sub> + C<sub>added</sub> = c<sub>free</sub> (1 + Σ<sub>p</sub> K<sub>p</sub> c<sub>OAm</sub><sup>p</sup>).
        Equilibrium only: a washed dot may stay out of equilibrium for a long time; these are driving forces.
      </p>
    {/if}

    {#if view === "synthesis" && growthUrl}
      {@render section("Growth and size selection", "several dots of this material and their stripped states in equilibrium with the MX and ML<sub>n</sub> monomers (solution-phase model of the report)")}
      <TrimmedPage src={growthUrl} autosize title="Synthesis thermodynamics" />
    {/if}

    {#if !model && !(view === "synthesis" && growthUrl) && !(view === "thermo" && reportUrl)}
      <div class="text-sm text-slate-500">No ladder data for this dot.</div>
    {/if}
  </div>
{/if}
