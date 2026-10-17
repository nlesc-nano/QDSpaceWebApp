<script>
  // Structure and vibrations tab: the relaxed structure at each level of theory (MACE-MH-1, DFT), and the
  // vibrations and thermochemistry of the MACE-MH-1 Hessian with IR / Raman intensities (DFT or g-xTB).
  // Data: props/{relax,structure,hessian,vibspec,spectra,report}.json and the relaxed geometries.
  import Plot from "./Plot.svelte";
  import { parseXyz, bonds, cnCounts, displacements, sameAtoms, histogram } from "./geom.js";

  let { files, geoms } = $props();         // findFiles(rec), geometries(rec)

  let data = $state({});
  let xyz = $state({});                    // {mace, dft: [..]} parsed geometries
  let status = $state("loading");
  $effect(() => {
    const f = files, g = geoms;
    status = "loading";
    const get = (u) => (u ? fetch(u).then((r) => (r.ok ? r.json() : null)).catch(() => null) : Promise.resolve(null));
    const txt = (u) => (u ? fetch(u).then((r) => (r.ok ? r.text() : null)).then((t) => (t ? parseXyz(t) : null)).catch(() => null) : Promise.resolve(null));
    Promise.all([get(f.relax), get(f.structure), get(f.hessian), get(f.vibspec), get(f.spectra), get(f.report),
                 txt(g.mace), ...g.dft.map((d) => txt(d.url))])
      .then(([relax, structure, hessian, vibspec, spectra, report, mace, ...dft]) => {
        if (f !== files) return;
        data = { relax, structure, hessian, vibspec, spectra, report };
        xyz = { mace, dft };
        status = "ready";
      });
  });

  // ---------------- levels of theory ----------------
  const MACE = { color: "#059669", label: "MACE-MH-1" };
  const DFTC = { color: "#0284c7" };
  let level = $state("both");              // "mace" | "dft" | "both"
  let dftI = $state(0);
  let dftG = $derived(geoms.dft[dftI] ? { ...geoms.dft[dftI], g: xyz.dft?.[dftI] } : null);
  let levels = $derived([
    ...(level !== "dft" ? [{ key: "mace", ...MACE, g: xyz.mace }] : []),
    ...(level !== "mace" && dftG?.g ? [{ key: "dft", color: DFTC.color, label: dftG.label, g: dftG.g }] : [])].filter((l) => l.g));
  const seg = (on) => `px-3 py-1.5 rounded-lg text-xs font-bold transition-all ${on ? "bg-white text-slate-900 shadow-sm" : "text-slate-500 hover:text-slate-800"}`;
  const fmt = (v, n = 2) => (v === null || v === undefined || !isFinite(v) ? "–" : Number(v).toFixed(n));
  const tex = (s) => String(s).replace(/\$\\Delta /g, "Δ").replace(/\\mathrm\{([^}]*)\}/g, "$1").replace(/_\{?([A-Za-z0-9,]+)\}?/g, "<sub>$1</sub>")
    .replace(/\^\{?(-?[0-9]+)\}?/g, "<sup>$1</sup>").replace(/\\alpha/g, "α").replace(/\\langle/g, "⟨").replace(/\\rangle/g, "⟩")
    .replace(/\\frac\{1\}\{2\}/g, "½").replace(/\\Sigma/g, "Σ").replace(/\\hbar\\omega/g, "ħω").replace(/\\Delta/g, "Δ").replace(/[$\\]/g, "");

  // ---------------- structure ----------------
  let roles = $derived(data.structure?.role || null);
  let cutoffs = $derived(data.structure?.cutoffs_A || null);
  let analysed = $derived(cutoffs ? levels.map((l) => {
    const r = sameAtoms(l.g, xyz.mace) ? roles : null;
    const b = bonds(l.g, cutoffs, r);
    return { ...l, b, cn: cnCounts(l.g, b.neighbours), roles: r };
  }) : []);
  let d0 = $derived(data.structure?.summary?.bulk_bond_A ?? null);

  // key numbers per level
  let maceCards = $derived.by(() => {
    const r = data.relax?.summary, s = data.structure?.summary;
    if (!r) return [];
    return [["Relaxation energy", `${fmt(r.relaxation_energy_meV_atom, 0)} meV/atom`, "(E<sub>relaxed</sub> − E<sub>start</sub>) / N"],
            ["RMSD start → relaxed", `${fmt(r.rmsd_A)} Å`, "after optimal superposition"],
            ["Bond topology", `${s?.bonds_broken ?? "–"} broken / ${s?.bonds_formed ?? "–"} formed`, "relaxed against start bond graph"],
            ["BFGS steps", `${r.steps}`, `to ${fmt(r.fmax_eV_A, 3)} eV/Å`]];
  });
  let dftCards = $derived.by(() => {
    if (!dftG?.g || !xyz.mace) return [];
    const same = sameAtoms(dftG.g, xyz.mace);
    const disp = same ? displacements(xyz.mace, dftG.g) : null;
    const rmsd = disp ? Math.sqrt(disp.reduce((a, d) => a + d * d, 0) / disp.length) : null;
    const mean = (g, p) => { const b = bonds(g, cutoffs, null).bonds[p] || []; return b.length ? b.reduce((a, x) => a + x.d, 0) / b.length : null; };
    const pairs = cutoffs ? Object.keys(cutoffs) : [];
    return [["RMSD from MACE-MH-1", disp ? `${fmt(rmsd, 3)} Å` : "–", disp ? "same atoms, after optimal superposition" : "different atom order"],
            ["Largest displacement", disp ? `${fmt(Math.max(...disp))} Å` : "–", "from the MACE-MH-1 geometry"],
            ...pairs.map((p) => { const m = mean(xyz.mace, p), d = mean(dftG.g, p);
              return [`Mean ${p} bond`, `${fmt(d, 3)} Å`, `MACE-MH-1 ${fmt(m, 3)} Å (${d && m ? (d - m >= 0 ? "+" : "−") + Math.abs(d - m).toFixed(3) : "–"})`]; })];
  });

  // relaxation: MACE energy trace / DFT displacements from the MACE geometry
  let relaxTab = $state("mace");
  let traceLayers = $derived.by(() => {
    const tr = data.relax?.energy_trace_eV;
    if (!tr) return [];
    const e1 = tr.at(-1), n = xyz.mace?.sym.length || 1;
    return [{ type: "line", color: MACE.color, pts: tr.map((e, i) => [i, ((e - e1) / n) * 1000]) }, { type: "hline", y: 0, color: "#94a3b8", dash: "" }];
  });
  let dispLayers = $derived.by(() => {
    if (!dftG?.g || !sameAtoms(dftG.g, xyz.mace)) return [];
    const d = displacements(xyz.mace, dftG.g), r = roles || d.map(() => "atom");
    const ROLE = { core: "#334155", surface: "#d97706", ligand: "#16a34a", atom: DFTC.color };
    return Object.keys(ROLE).filter((k) => r.includes(k)).map((k) => ({
      type: "bars", bars: histogram(d.filter((_, i) => r[i] === k), 0, 0.2, 0.005).map((h) => ({ x: h.x, y0: 0, y1: h.c, w: 0.0045, color: ROLE[k], opacity: 0.7 })), key: k, color: ROLE[k] }));
  });

  // bond lengths: one histogram per pair (Cd-Se split into core and shell when the roles are known)
  let bondPanels = $derived.by(() => {
    if (!analysed.length) return [];
    const out = [];
    for (const p of Object.keys(cutoffs)) {
      const split = analysed.some((a) => a.roles) && p.split("-").every((e) => (data.structure?.bulk_cn || {})[e]);
      const groups = split ? [["core", true], ["shell", false]] : [["all", null]];
      for (const [gname, core] of groups) {
        const series = analysed.map((a) => ({ a, vals: a.b.bonds[p].filter((x) => core === null || x.core === core).map((x) => x.d) }));
        const all = series.flatMap((s) => s.vals);
        if (!all.length) continue;
        const lo = Math.floor((Math.min(...all) - 0.05) * 50) / 50, hi = Math.ceil((Math.max(...all) + 0.05) * 50) / 50;
        const n = series.length, w = 0.02 / (n + 0.5);
        out.push({ title: `${p}${gname !== "all" ? ` · ${gname}` : ""}`, lo, hi,
          layers: [...series.map((s, j) => ({ type: "bars", bars: histogram(s.vals, lo, hi, 0.02).map((h) => ({ x: h.x - 0.01 + w * (j + 0.75), y0: 0, y1: h.c, w, color: s.a.color, opacity: 0.85 })) })),
                   ...(d0 && p.split("-").every((e) => (data.structure?.bulk_cn || {})[e]) ? [{ type: "vline", x: d0, color: "#0f172a", label: `bulk ${d0.toFixed(3)} Å` }] : [])],
          stats: series.map((s) => ({ a: s.a, n: s.vals.length, mean: s.vals.reduce((x, y) => x + y, 0) / (s.vals.length || 1) })) });
      }
    }
    return out;
  });

  // coordination numbers per element
  let cnPanels = $derived.by(() => {
    if (!analysed.length) return [];
    const els = [...new Set(analysed.flatMap((a) => Object.keys(a.cn)))];
    return els.map((e) => {
      const cns = [...new Set(analysed.flatMap((a) => Object.keys(a.cn[e] || {}).map(Number)))].sort((x, y) => x - y);
      const n = analysed.length, w = 0.8 / n;
      return { e, cns, layers: analysed.map((a, j) => ({ type: "bars", bars: cns.map((c) => ({ x: c - 0.4 + w * (j + 0.5), y0: 0, y1: a.cn[e]?.[c] || 0, w: w * 0.9, color: a.color, opacity: 0.85 })) })),
               ymax: Math.max(1, ...analysed.flatMap((a) => Object.values(a.cn[e] || {}))), bulk: data.structure?.bulk_cn?.[e] };
    });
  });

  // ---------------- vibrations ----------------
  let vibSource = $derived(data.vibspec?.summary?.source || data.vibspec?.provenance?.engine || "");
  let vibDft = $derived(vibSource.startsWith("CP2K"));
  let vibCards = $derived.by(() => {
    const h = data.hessian?.summary, v = data.vibspec?.summary;
    if (!h) return [];
    return [["Imaginary modes", `${h.n_imaginary}`, "none at a true minimum"],
            ["ZPE", `${fmt(h.zpe_eV, 3)} eV`, "½ Σ<sub>k</sub> ħω<sub>k</sub>"],
            ["S<sub>vib</sub>, 300 K", `${fmt(h.S_vib_300K_meV_K)} meV/K`, "quasi-RRHO"],
            ["F<sub>vib</sub>, 300 K", `${fmt(h.F_vib_300K_eV, 3)} eV`, "U<sub>vib</sub> − TS<sub>vib</sub>"],
            ...(v ? [["Point group", v.point_group, "of the relaxed dot"],
                     ["Polarisability α<sub>iso</sub>", `${fmt(v.alpha_iso_A3, 0)} Å³`, vibDft ? "DFT PBE, finite field" : ""],
                     ["Strongest IR peak", `${fmt(v.top_ir?.[0]?.[0], 0)} cm⁻¹ (${v.top_ir?.[0]?.[1]})`, `${fmt(v.top_ir?.[0]?.[2], 0)} km/mol`],
                     ["Strongest Raman", `${fmt(v.top_raman?.[0]?.[0], 0)} cm⁻¹ (${v.top_raman?.[0]?.[1]})`, `${fmt(v.top_raman?.[0]?.[2], 0)} Å⁴/amu`]] : [])];
  });
  let vdosBy = $state("element");
  const ECOL = ["#7c3aed", "#db2777", "#16a34a", "#ea580c", "#0891b2"];
  const RCOL = { core: "#2a78d6", surface: "#eb6834", ligand: "#1baf7a" };
  let vdos = $derived.by(() => {
    const v = data.spectra?.vdos;
    if (!v) return { layers: [], legend: [] };
    const proj = v.projected?.[vdosBy] || {};
    const keys = Object.keys(proj);
    const col = (k, i) => (vdosBy === "role" ? RCOL[k] || ECOL[i] : ECOL[i % ECOL.length]);
    return { layers: [{ type: "line", color: "#94a3b8", dash: "4 3", pts: v.cm1.map((x, i) => [x, v.total[i]]) },
                      ...keys.map((k, i) => ({ type: "line", color: col(k, i), pts: v.cm1.map((x, j) => [x, proj[k][j]]) }))],
             legend: [{ color: "#94a3b8", dash: "4 3", label: "total" }, ...keys.map((k, i) => ({ color: col(k, i), label: k }))],
             ymax: Math.max(...v.total) * 1.05, xmax: v.cm1.at(-1) };
  });
  let fwhm = $state(8);
  const lorentz = (modes, key, fw, xmax) => {
    const g = fw / 2, out = [];
    for (let x = 0; x <= xmax; x += 1) {
      let s = 0;
      for (const m of modes) s += (m[key] || 0) * (g / Math.PI) / ((x - m.nu) ** 2 + g * g);
      out.push([x, s]);
    }
    return out;
  };
  let modes = $derived(data.vibspec?.modes || []);
  let xmaxVib = $derived(Math.ceil(((data.hessian?.summary?.highest_cm1 || 300) + 30) / 50) * 50);
  const spec = (key, color) => {
    if (!modes.length) return { layers: [], ymax: 1 };
    const sticks = modes.filter((m) => (m[key] || 0) > 0);
    const smax = Math.max(...sticks.map((m) => m[key]), 1e-9);
    const env = lorentz(modes, key, fwhm, xmaxVib), emax = Math.max(...env.map((p) => p[1]), 1e-9);
    const bulk = data.vibspec?.summary?.bulk_gamma_optical_cm1?.[0];
    return { layers: [{ type: "bars", bars: sticks.map((m) => ({ x: m.nu, y0: 0, y1: m[key] / smax, w: 0.8, color: "#94a3b8", opacity: 0.8 })) },
                      { type: "line", color, pts: env.map(([x, y]) => [x, y / emax]) },
                      ...(bulk ? [{ type: "vline", x: bulk, color: "#0f172a", label: `bulk TO ${bulk.toFixed(0)} cm⁻¹` }] : [])], ymax: 1.08 };
  };
  let ir = $derived(spec("ir_km_mol", "#dc2626"));
  let raman = $derived(spec("raman_A4_amu", "#7c3aed"));
  const modeHover = (key, unit) => (x) => {
    if (!modes.length) return null;
    const m = modes.reduce((b, q) => (Math.abs(q.nu - x) < Math.abs(b.nu - x) ? q : b), modes[0]);
    return { x: m.nu, title: `${m.nu.toFixed(1)} cm⁻¹ · ${m.irrep}${m.class ? " · " + m.class : ""}`,
             rows: [{ text: `${key === "ir_km_mol" ? "IR" : "Raman"} ${fmt(m[key], 2)} ${unit}` },
                    { text: `core ${(100 * (m.core || 0)).toFixed(0)}% · surface ${(100 * (m.surface || 0)).toFixed(0)}% · ligand ${(100 * (m.ligand || 0)).toFixed(0)}%` }] };
  };
  let thermo = $derived(data.spectra?.thermo || null);
</script>

{#snippet cards(list, tone)}
  <div class="grid grid-cols-2 md:grid-cols-4 gap-3">
    {#each list as [label, value, note]}
      <div class="rounded-xl border p-3 {tone}">
        <div class="text-[11px] text-slate-500">{@html label}</div>
        <div class="text-lg font-bold text-slate-900">{@html value}</div>
        {#if note}<div class="text-[11px] text-slate-500 leading-snug">{@html note}</div>{/if}
      </div>
    {/each}
  </div>
{/snippet}
{#snippet section(title, subtitle)}
  <div class="flex items-baseline gap-3 flex-wrap border-b border-slate-200 pb-2 mt-2">
    <h3 class="text-base font-bold text-slate-900">{@html title}</h3>
    <span class="text-xs text-slate-500">{@html subtitle}</span>
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
{:else}
  <div class="flex flex-col gap-6 pb-6">
    <!-- structure -->
    {@render section("Structure and relaxation", "the relaxed dot at each level of theory")}
    <div class="flex flex-wrap items-center gap-3 text-xs">
      <span class="font-bold text-slate-400">Level</span>
      <div class="flex gap-1 bg-slate-100 p-1 rounded-xl border border-slate-200">
        <button class={seg(level === "mace")} onclick={() => (level = "mace")}>MACE-MH-1</button>
        {#if geoms.dft.length}<button class={seg(level === "dft")} onclick={() => (level = "dft")}>DFT</button>
          <button class={seg(level === "both")} onclick={() => (level = "both")}>both</button>{/if}
      </div>
      {#if geoms.dft.length > 1 && level !== "mace"}
        <div class="flex gap-1 bg-slate-100 p-1 rounded-xl border border-slate-200">
          {#each geoms.dft as d, i}<button class={seg(dftI === i)} onclick={() => (dftI = i)}>{d.short}{d.chain ? "" : " (earlier)"}</button>{/each}
        </div>
      {/if}
    </div>
    {#if level !== "dft"}{@render cards(maceCards, "bg-emerald-50/40 border-emerald-100")}{/if}
    {#if level !== "mace" && dftCards.length}
      <div class="text-xs font-bold text-sky-800">{dftG?.label}, relaxed from the MACE-MH-1 geometry</div>
      {@render cards(dftCards, "bg-sky-50/40 border-sky-100")}
    {/if}

    <div>
      <div class="flex items-center gap-3 mb-1">
        <h4 class="text-xs font-bold text-slate-500 uppercase tracking-widest">Relaxation</h4>
        <div class="flex gap-1 bg-slate-100 p-1 rounded-xl border border-slate-200">
          <button class={seg(relaxTab === "mace")} onclick={() => (relaxTab = "mace")}>MACE-MH-1</button>
          {#if dispLayers.length}<button class={seg(relaxTab === "dft")} onclick={() => (relaxTab = "dft")}>{dftG?.label}</button>{/if}
        </div>
      </div>
      {#if relaxTab === "mace" || !dispLayers.length}
        <p class="text-[11px] text-slate-500 mb-1">Energy per atom above the relaxed structure along the BFGS relaxation from the builder's geometry.</p>
        <Plot layers={traceLayers} xDomain={[0, Math.max(1, (data.relax?.energy_trace_eV?.length || 2) - 1)]}
              yDomain={[-2, Math.max(10, ...(traceLayers[0]?.pts.map((p) => p[1]) || [10])) * 1.05]} xLabel="BFGS step" yLabel={"E − E_{relaxed} (meV/atom)"} height={260} ariaLabel="relaxation" />
      {:else}
        <p class="text-[11px] text-slate-500 mb-1">How far each atom moves from the MACE-MH-1 to the {dftG?.label} geometry (after optimal superposition), by where it sits.</p>
        <Plot layers={dispLayers} xDomain={[0, 0.2]} yDomain={[0, Math.max(1, ...dispLayers.flatMap((l) => l.bars.map((b) => b.y1))) * 1.08]}
              xLabel="displacement (Å)" yLabel="atoms" height={260} ariaLabel="displacements"
              legend={dispLayers.map((l) => ({ kind: "box", color: l.color, label: l.key }))} />
      {/if}
    </div>

    {#each bondPanels as bp}
      <div>
        {@render panel(`Bond lengths · ${bp.title}`, "")}
        <div class="flex flex-wrap gap-4 text-xs text-slate-600 mb-1">
          {#each bp.stats as s}<span class="flex items-center gap-1.5"><span class="w-3 h-3 rounded-sm inline-block" style="background:{s.a.color}"></span>{s.a.label}: {s.n} bonds, mean <b class="font-mono">{s.mean.toFixed(3)} Å</b></span>{/each}
        </div>
        <Plot layers={bp.layers} xDomain={[bp.lo, bp.hi]} yDomain={[0, Math.max(1, ...bp.layers.filter((l) => l.bars).flatMap((l) => l.bars.map((b) => b.y1))) * 1.1]}
              xLabel="bond length (Å)" yLabel="bonds" height={240} ariaLabel="bond lengths" />
      </div>
    {/each}
    {#if cnPanels.length}
      <div>
        {@render panel("Coordination numbers", "Neighbours within the bond cutoffs, per element; the dashed line is the bulk coordination.")}
        <div class="grid grid-cols-1 md:grid-cols-3 gap-4">
          {#each cnPanels as cp}
            <div>
              <div class="text-xs font-bold text-slate-700 mb-1">{cp.e}</div>
              <Plot layers={[...cp.layers, ...(cp.bulk ? [{ type: "vline", x: cp.bulk, color: "#0f172a", label: "bulk" }] : [])]}
                    xDomain={[Math.min(...cp.cns) - 0.7, Math.max(...cp.cns, cp.bulk || 0) + 0.7]} yDomain={[0, cp.ymax * 1.1]}
                    xTicks={[...new Set([...cp.cns, ...(cp.bulk ? [cp.bulk] : [])])].sort((a, b) => a - b).map((v) => ({ v, label: String(v) }))}
                    xLabel="coordination number" yLabel="atoms" height={220} ariaLabel="coordination" />
            </div>
          {/each}
        </div>
      </div>
    {/if}
    <div class="flex flex-wrap gap-4 text-[11px] text-slate-500">
      {#each levels as l}<span class="flex items-center gap-1.5"><span class="w-3 h-3 rounded-sm inline-block" style="background:{l.color}"></span>{l.label}</span>{/each}
      <span>Bonds: opposite-element pairs closer than {#if cutoffs}{Object.entries(cutoffs).map(([p, c]) => `${p} ${c.toFixed(2)} Å`).join(", ")}{/if} (the structure step's cutoffs).</span>
    </div>

    <!-- vibrations -->
    {@render section("Vibrations and thermochemistry", "harmonic modes of the MACE-MH-1 Hessian")}
    <p class="text-sm text-slate-700 max-w-4xl leading-relaxed">
      The normal modes come from the MACE-MH-1 Hessian (finite differences) at the MACE-MH-1 minimum, with translations and rotations
      projected out; the thermochemistry treats them as harmonic oscillators, with Grimme's quasi-RRHO entropy for the soft modes.
      {#if vibDft}IR and Raman are a hybrid: the peak positions are the MACE-MH-1 frequencies, the intensities the DFT PBE Born charges and
        Raman tensors (13 finite-field single points) projected on those modes.
      {:else if modes.length}IR and Raman intensities: g-xTB dipole and polarisability derivatives along the MACE-MH-1 modes.{/if}
      The theory is in the documentation.
    </p>
    {@render cards(vibCards, "bg-slate-50 border-slate-100")}
    {#if vdos.layers.length}
      <div>
        <div class="flex items-center gap-3 mb-1">
          <h4 class="text-xs font-bold text-slate-500 uppercase tracking-widest">Vibrational density of states</h4>
          <div class="flex gap-1 bg-slate-100 p-1 rounded-xl border border-slate-200">
            <button class={seg(vdosBy === "element")} onclick={() => (vdosBy = "element")}>by element</button>
            <button class={seg(vdosBy === "role")} onclick={() => (vdosBy = "role")}>by site</button>
          </div>
        </div>
        <Plot layers={vdos.layers} xDomain={[0, vdos.xmax]} yDomain={[0, vdos.ymax]} xLabel={"frequency (cm^{-1})"} yLabel="VDOS" height={280}
              ariaLabel="vibrational density of states" legend={vdos.legend} />
      </div>
    {/if}
    {#if modes.length}
      <div class="flex flex-wrap items-center gap-3 text-xs">
        <span class="font-bold text-slate-500">Broadening (FWHM)</span>
        <input type="range" min="2" max="30" step="1" value={fwhm} oninput={(e) => (fwhm = +e.currentTarget.value)} class="w-40 accent-brand-600" />
        <span class="font-mono text-slate-900">{fwhm} cm⁻¹</span>
      </div>
      <div>
        {@render panel("IR spectrum", `Sticks: each mode's intensity; line: Lorentzian broadening. Hover for the mode's symmetry and where it lives.`)}
        <Plot layers={ir.layers} xDomain={[0, xmaxVib]} yDomain={[0, ir.ymax]} xLabel={"frequency (cm^{-1})"} yLabel="IR intensity (normalised)" height={280}
              hover={modeHover("ir_km_mol", "km/mol")} ariaLabel="IR spectrum" />
      </div>
      <div>
        {@render panel("Raman spectrum (non-resonant)", "")}
        <Plot layers={raman.layers} xDomain={[0, xmaxVib]} yDomain={[0, raman.ymax]} xLabel={"frequency (cm^{-1})"} yLabel="Raman activity (normalised)" height={280}
              hover={modeHover("raman_A4_amu", "Å⁴/amu")} ariaLabel="Raman spectrum" />
      </div>
    {/if}
    {#if thermo}
      <div>
        {@render panel("Vibrational thermochemistry", "")}
        <Plot layers={[{ type: "line", color: "#2563eb", pts: thermo.T.map((t, i) => [t, thermo.U_vib_eV[i]]) },
                       { type: "line", color: "#dc2626", pts: thermo.T.map((t, i) => [t, thermo.F_vib_eV[i]]) }]}
              xDomain={[thermo.T[0], thermo.T.at(-1)]} yDomain={[Math.min(...thermo.F_vib_eV) * 1.05, Math.max(...thermo.U_vib_eV) * 1.05]}
              xLabel="T (K)" yLabel="energy (eV)" height={260} ariaLabel="vibrational energies"
              legend={[{ color: "#2563eb", label: "U<sub>vib</sub>" }, { color: "#dc2626", label: "F<sub>vib</sub> = U<sub>vib</sub> − TS<sub>vib</sub>" }]} />
        <Plot layers={[{ type: "line", color: "#7c3aed", pts: thermo.T.map((t, i) => [t, thermo.S_vib_meV_K[i]]) },
                       { type: "line", color: "#0891b2", pts: thermo.T.map((t, i) => [t, thermo.Cv_meV_K[i]]) }]}
              xDomain={[thermo.T[0], thermo.T.at(-1)]} yDomain={[0, Math.max(...thermo.S_vib_meV_K) * 1.05]}
              xLabel="T (K)" yLabel="meV/K" height={260} ariaLabel="vibrational entropy and heat capacity"
              legend={[{ color: "#7c3aed", label: "S<sub>vib</sub> (quasi-RRHO)" }, { color: "#0891b2", label: "C<sub>v</sub>" }]} />
      </div>
    {/if}
  </div>
{/if}
