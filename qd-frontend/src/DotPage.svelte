<script>
  // One dot's page (#/qd/<id>/<tab>): the 3D viewer and what the dot is on top, then one tab per kind of
  // property. Library records come from library_index.json; computed data from recordFiles.js.
  import Viewer from "./Viewer.svelte";
  import ThermoView from "./ligand/ThermoView.svelte";
  import SynthesisView from "./ligand/SynthesisView.svelte";
  import TrimmedPage from "./ligand/TrimmedPage.svelte";
  import StructureView from "./ligand/StructureView.svelte";
  import { findFiles, dftPages, geometries } from "./ligand/recordFiles.js";
  import { buildModel, solventModel, unitCosts, decompositionPerCore, UNITS } from "./ligand/model.js";

  let { id, tab = "thermo", onTab = () => {}, onBack = () => {} } = $props();

  let rec = $state(null);
  let status = $state("loading");
  let files = $state(null);
  let key = $state(null);                // key numbers from the ladder and stability data

  $effect(() => {
    const want = id;
    status = "loading";
    rec = null; files = null; key = null;
    fetch("/library_index.json").then((r) => r.json()).then(async (index) => {
      const r = (index.structures || []).find((s) => s.id === want);
      if (want !== id) return;
      if (!r) { status = "missing"; return; }
      rec = r;
      status = "ready";
      const f = await findFiles(r);
      if (want !== id) return;
      files = f;
      loadKeyNumbers(f, want);
      pickStage("ml");
    }).catch(() => { status = "error"; });
  });

  async function loadKeyNumbers(f, want) {
    if (!f.det) return;
    const get = (u) => (u ? fetch(u).then((r) => (r.ok ? r.json() : null)).catch(() => null) : Promise.resolve(null));
    const [det, stab, cubes] = await Promise.all([get(f.det), get(f.stab), get(f.cubes)]);
    if (want !== id || !det) return;
    const m = buildModel(det, stab), sm = solventModel(cubes), i = m.T.indexOf(300);
    const sol = sm ? { sm, eps: 2.38 } : null, first = m.units.findIndex((u) => u.dk === 0);
    key = {
      unit: m.unit, k0: m.k0, m: m.m, core: m.core,
      removal: first >= 0 ? { gas: unitCosts(m, i, null)[first], sol: sol ? unitCosts(m, i, sol)[first] : null } : null,
      // dot -> m MX(bulk) + k0 ML_n: reaction free energy per MX, positive = the dot is stable
      dec: m.G_dot ? { gas: -decompositionPerCore(m, i, 1, null), sol: sol ? -decompositionPerCore(m, i, 1, sol) : null } : null,
    };
  }

  // ---------------- viewer ----------------
  // modes: ML relaxed (MACE-MH-1), DFT relaxed (one per functional), orbitals of the DFT chain
  let mode = $state("ml");
  let dftIndex = $state(0);
  let orbital = $state("");
  let xyz = $state("");
  let dataUrl = $state("");
  let loadingStage = $state(false);
  let mlStage = $derived(rec?.stages.find((s) => s.stage === "geo_opt" && s.code === "MACE") || rec?.stages.find((s) => s.stage === "start"));
  let dftStages = $derived((rec?.stages || []).filter((s) => s.stage === "geo_opt" && s.code && s.code !== "MACE")
    .map((s) => ({ s, label: `DFT ${s.functional || "(functional unknown)"}${s.file.endsWith("/dft/relaxed.xyz") ? " · this pipeline" : ""}` }))
    .sort((a, b) => b.s.file.endsWith("/dft/relaxed.xyz") - a.s.file.endsWith("/dft/relaxed.xyz")));
  let orbitalUrls = $derived(rec?.dft?.orbitals ? Object.fromEntries(Object.entries(rec.dft.orbitals).map(([k, p]) => [k, `/${p}`])) : null);
  const ORBITALS = ["HOMO-1", "HOMO", "LUMO", "LUMO+1"];

  async function load(st) {
    if (!st) return;
    loadingStage = true;
    try {
      const text = (await (await fetch(`/${st.file}`)).text()).replace(/\r\n/g, "\n");
      xyz = text; dataUrl = `/${st.file}`;
    } catch { xyz = ""; dataUrl = ""; } finally { loadingStage = false; }
  }
  async function pickStage(m, i = dftIndex) {
    mode = m; orbital = "";
    if (m === "ml") await load(mlStage);
    else if (m === "dft") { dftIndex = i; await load(dftStages[i]?.s); }
  }
  async function pickOrbital(label) {
    const chain = dftStages.findIndex((d) => d.s.file.endsWith("/dft/relaxed.xyz"));
    if (chain >= 0) await load(dftStages[chain].s);
    mode = "orbitals"; orbital = label;
  }

  // ---------------- what the dot is ----------------
  const parts = (f) => String(f).split(/(\d+)/).filter(Boolean).map((t) => ({ t, sub: /^\d+$/.test(t) }));
  const html = (s) => String(s).replace(/(\d+)/g, "<sub>$1</sub>");
  const fe = (v) => (v === null || v === undefined || !isFinite(v) ? "–" : `${v >= 0 ? "+" : "−"}${Math.abs(v).toFixed(2)}`);
  let facets = $derived(rec?.origin?.recipe?.facets || []);
  let dotTitle = $derived(key?.core ? `[${html(key.core)}]<sub>${key.m}</sub>[${html(key.unit)}]<sub>${key.k0}</sub>` : null);
  let excitons = $derived(rec?.dft?.summary || null);
  // the same QDEX run with the trap filter on (published next to webapp.json)
  let excitonsTrap = $state(null);
  $effect(() => {
    const w = rec?.dft?.webapp;
    excitonsTrap = null;
    if (!w) return;
    fetch(`/${w.replace(/webapp\.json$/, "webapp_traps_filtered.json")}`)
      .then((r) => (r.ok && (r.headers.get("content-type") || "").includes("json") ? r.json() : null))
      .then((d) => { if (rec?.dft?.webapp === w) excitonsTrap = d?.summary || null; }).catch(() => {});
  });
  let pages = $derived(rec ? dftPages(rec) : { electronic: [], excited: [] });

  // ---------------- tabs ----------------
  let tabs = $derived([
    ["structure", "Structure and vibrations", !!(files?.relax || files?.groundState)],
    ["thermo", "Thermodynamics", !!files?.det],
    ["synthesis", "Synthesis and post-synthesis", !!files?.det],
    ["electronic", "Electronic structure", pages.electronic.length > 0],
    ["excited", "Excited states", pages.excited.length > 0],
  ].filter(([, , ok]) => ok));
  let active = $derived(tabs.some(([k]) => k === tab) ? tab : tabs[0]?.[0]);

  // electronic structure / excited states: functional and SF / SOC
  let fn = $state(0);
  let spin = $state("sf");
  let traps = $state("included");        // excited states: "included" | "filtered"
  let group = $derived((active === "excited" ? pages.excited : pages.electronic)[fn] || (active === "excited" ? pages.excited : pages.electronic)[0]);
  let hasTrapPages = $derived(active === "excited" && !!(group?.sf_trap || group?.soc_trap));
  let pageUrl = $derived(group ? (hasTrapPages && traps === "filtered" ? group[`${spin}_trap`] : null) || group[spin] || group.sf || group.soc : null);

  const pill = (on) => `px-3 py-1.5 rounded-full text-xs font-bold border transition-all ${on ? "bg-brand-600 text-white border-brand-600 shadow-sm" : "bg-white text-slate-600 border-slate-200 hover:border-brand-300 hover:text-brand-700"}`;
  const seg = (on) => `px-3 py-1.5 rounded-lg text-xs font-bold transition-all ${on ? "bg-white text-slate-900 shadow-sm" : "text-slate-500 hover:text-slate-800"}`;
</script>

{#if status === "loading"}
  <div class="flex items-center justify-center py-32 text-brand-600 text-sm font-bold">
    <div class="w-6 h-6 border-4 border-brand-500 border-t-transparent rounded-full animate-spin mr-3"></div>Loading…
  </div>
{:else if status !== "ready"}
  <div class="max-w-3xl mx-auto py-24 text-center text-slate-500">
    <p class="font-bold mb-2">{status === "missing" ? `No dot ${id} in the library.` : "The library index could not be loaded."}</p>
    <button class="text-brand-700 font-bold" onclick={onBack}>Back to the library</button>
  </div>
{:else}
  <main class="max-w-[1500px] mx-auto px-4 md:px-6 py-6 flex flex-col gap-6">
    <!-- header -->
    <div class="flex flex-wrap items-end justify-between gap-3">
      <div class="flex flex-col gap-1">
        <button class="text-xs font-bold text-slate-500 hover:text-brand-700 self-start" onclick={onBack}>← Library</button>
        <h1 class="font-heading text-2xl md:text-3xl font-bold text-slate-900">
          {#if dotTitle}{@html dotTitle}{:else}{#each parts(rec.formula) as p}{#if p.sub}<sub>{p.t}</sub>{:else}{p.t}{/if}{/each}{/if}
        </h1>
        <div class="text-sm text-slate-500">
          {#each parts(rec.formula) as p}{#if p.sub}<sub>{p.t}</sub>{:else}{p.t}{/if}{/each}
          · {rec.material} {rec.phase} · {rec.size?.d_nm?.toFixed(2)} nm · {rec.surface}{rec.reconstruction_pattern ? ` (${rec.reconstruction_pattern})` : ""}
        </div>
      </div>
      <div class="flex flex-wrap gap-1.5">
        {#if mlStage?.code === "MACE"}<span class="px-2 py-0.5 rounded border text-[11px] font-bold bg-emerald-50 text-emerald-800 border-emerald-200">MACE-MH-1</span>{/if}
        {#each [...new Set(dftStages.map((d) => `DFT ${d.s.functional || "(functional unknown)"}`))] as m}<span class="px-2 py-0.5 rounded border text-[11px] font-bold bg-sky-50 text-sky-800 border-sky-200">{m}</span>{/each}
        {#if rec.dft?.method?.qdex_commit || pages.excited.length}<span class="px-2 py-0.5 rounded border text-[11px] font-bold bg-violet-50 text-violet-800 border-violet-200">QDEX</span>{/if}
      </div>
    </div>

    <!-- viewer and what the dot is -->
    <div class="grid grid-cols-1 lg:grid-cols-[minmax(0,1.75fr)_minmax(0,1fr)] gap-6">
      <div class="bg-white rounded-[1.5rem] border border-slate-100 shadow-sm p-4 flex flex-col gap-3">
        <div class="flex flex-wrap items-center gap-2">
          <div class="flex gap-1 bg-slate-100 p-1 rounded-xl border border-slate-200">
            <button class={seg(mode === "ml")} onclick={() => pickStage("ml")}>ML relaxed</button>
            {#if dftStages.length}<button class={seg(mode === "dft")} onclick={() => pickStage("dft")}>DFT relaxed</button>{/if}
            <button class={seg(false)} disabled title="Next step: the binding-site map">Binding map</button>
            {#if orbitalUrls}<button class={seg(mode === "orbitals")} onclick={() => pickOrbital(orbital || "HOMO")}>Orbitals</button>{/if}
          </div>
          {#if mode === "dft" && dftStages.length > 1}
            <select class="text-xs border border-slate-200 rounded-lg px-2 py-1.5" value={dftIndex} onchange={(e) => pickStage("dft", +e.currentTarget.value)}>
              {#each dftStages as d, i}<option value={i}>{d.label}</option>{/each}
            </select>
          {:else if mode === "dft"}
            <span class="text-xs text-slate-500">{dftStages[0]?.label}</span>
          {/if}
          {#if mode === "orbitals" && orbitalUrls}
            <div class="flex gap-1">
              {#each ORBITALS.filter((o) => orbitalUrls[o]) as o}<button class={pill(orbital === o)} onclick={() => pickOrbital(o)}>{o}</button>{/each}
            </div>
          {/if}
        </div>
        <div class="relative h-[460px] rounded-xl overflow-hidden">
          <Viewer {xyz} isMD={false} {dataUrl} sizeMetrics={rec.size ?? null} orbitals={orbitalUrls} orbital={mode === "orbitals" ? orbital : ""}
                  onOrbital={(o) => pickOrbital(o)} activeViewer={mode === "orbitals" ? "3dmol" : "molstar"} />
          {#if loadingStage}<div class="absolute inset-0 bg-slate-900/30 flex items-center justify-center"><div class="w-8 h-8 border-4 border-brand-500 border-t-transparent rounded-full animate-spin"></div></div>{/if}
        </div>
      </div>

      <div class="bg-white rounded-[1.5rem] border border-slate-100 shadow-sm p-5 flex flex-col gap-4 text-sm">
        <h2 class="font-heading text-lg font-bold text-slate-900">About this dot</h2>
        <dl class="grid grid-cols-[auto_1fr] gap-x-4 gap-y-2">
          <dt class="text-slate-500">Core</dt>
          <dd class="font-medium text-slate-900">{#if rec.core}{@html Object.entries(rec.core).map(([e, n]) => `${e}<sub>${n}</sub>`).join("")}{/if} · {rec.centre}-centred · d = {rec.size?.d_nm?.toFixed(2)} nm</dd>
          <dt class="text-slate-500">Facets</dt>
          <dd class="text-slate-900">{#if facets.length}{facets.map((f) => `{${f.hkl}}`).join(" · ")}{:else}–{/if}</dd>
          <dt class="text-slate-500">Ligand shell</dt>
          <dd class="text-slate-900">{#if key}{@html `${key.k0} ${html(key.unit)}`} ({@html Object.entries(rec.ligands || {}).map(([e, n]) => `${n} ${e}`).join(", ")}){:else}{@html Object.entries(rec.ligands || {}).map(([e, n]) => `${e}<sub>${n}</sub>`).join(" ")}{/if}</dd>
          <dt class="text-slate-500">Surface</dt>
          <dd class="text-slate-900">{rec.surface}{rec.reconstruction_pattern ? ` · ${rec.reconstruction_pattern}` : ""}</dd>
        </dl>
        <div class="border-t border-slate-100 pt-3 flex flex-col gap-2">
          <div class="flex items-center gap-2"><span class="text-xs font-bold text-slate-700">Lowest exciton</span>
            <span class="px-1.5 py-0.5 rounded border text-[10px] font-bold bg-sky-50 text-sky-800 border-sky-200">DFT {rec.dft?.method?.functional || "PBE"} + QDEX</span></div>
          <div class="grid grid-cols-2 gap-3">
            {#each [["spin free", "spin_free"], ["with SOC", "soc"]] as [label, k]}
              <div class="bg-slate-50 rounded-xl p-3">
                <div class="text-[11px] text-slate-500">{label}</div>
                <div class="text-lg font-bold text-slate-900">{excitons?.[k]?.lowest_ev ? `${excitons[k].lowest_ev.toFixed(2)} eV` : "–"}
                  {#if excitons?.[k]?.lowest_f !== undefined}<span class="text-xs font-medium text-slate-500">f = {excitons[k].lowest_f.toFixed(2)}</span>{/if}</div>
                {#if excitonsTrap?.[k]?.lowest_ev}<div class="text-xs text-slate-500">surface traps filtered out: <b class="text-slate-700">{excitonsTrap[k].lowest_ev.toFixed(2)} eV</b>, f = {excitonsTrap[k].lowest_f.toFixed(2)}</div>{/if}
              </div>
            {/each}
          </div>
        </div>
        {#if key}
          <div class="border-t border-slate-100 pt-3 flex flex-col gap-2">
            <div class="flex items-center gap-2 flex-wrap"><span class="text-xs font-bold text-slate-700">Stability, 300 K, 1 M</span>
              <span class="px-1.5 py-0.5 rounded border text-[10px] font-bold bg-emerald-50 text-emerald-800 border-emerald-200">MACE-MH-1 free energies</span>
              {#if key.dec?.sol !== null}<span class="px-1.5 py-0.5 rounded border text-[10px] font-bold bg-teal-50 text-teal-800 border-teal-200">toluene: DFT PBE COSMO</span>{/if}</div>
            <div class="grid grid-cols-2 gap-3">
              <div class="bg-slate-50 rounded-xl p-3">
                <div class="text-[11px] text-slate-500 leading-snug">{@html dotTitle} → {key.m} {@html html(key.core)}(bulk) + {key.k0} {@html html(key.unit)}, per {@html html(key.core)}</div>
                <div class="text-lg font-bold text-slate-900">{fe(key.dec?.gas)} <span class="text-xs font-medium text-slate-500">eV vacuum</span></div>
                {#if key.dec?.sol !== null}<div class="text-sm font-bold text-teal-800">{fe(key.dec.sol)} <span class="text-xs font-medium text-slate-500">eV toluene</span></div>{/if}
              </div>
              <div class="bg-slate-50 rounded-xl p-3">
                <div class="text-[11px] text-slate-500 leading-snug">first {@html html(key.unit)} removal: k = {key.k0} → {key.k0 - 1} + {@html html(key.unit)}</div>
                <div class="text-lg font-bold text-slate-900">{fe(key.removal?.gas)} <span class="text-xs font-medium text-slate-500">eV vacuum</span></div>
                {#if key.removal?.sol !== null}<div class="text-sm font-bold text-teal-800">{fe(key.removal.sol)} <span class="text-xs font-medium text-slate-500">eV toluene</span></div>{/if}
              </div>
            </div>
            <p class="text-[11px] text-slate-500">Reaction free energies with the dot on the left: positive, the dot is stable against that reaction.</p>
          </div>
        {/if}
      </div>
    </div>

    <!-- tabs -->
    <div class="sticky top-[60px] z-20 bg-slate-50/95 backdrop-blur py-2 -mx-1 px-1 flex flex-wrap gap-2">
      {#each tabs as [k, label]}<button class={pill(active === k)} onclick={() => onTab(k)}>{label}</button>{/each}
    </div>

    <div class="bg-white rounded-[1.5rem] border border-slate-100 shadow-sm p-5 md:p-6 min-h-[600px]">
      {#if active === "structure"}
        {#if files.relax}
          <StructureView {files} geoms={geometries(rec)} />
        {:else}
          <TrimmedPage src={files.groundState} keep={(t) => !t.startsWith("Stability")} autosize title="Structure and vibrations" />
        {/if}
      {:else if active === "thermo"}
        <ThermoView detUrl={files.det} stabUrl={files.stab} cubesUrl={files.cubes} />
      {:else if active === "synthesis"}
        <SynthesisView detUrl={files.det} stabUrl={files.stab} cubesUrl={files.cubes} />
      {:else if active === "electronic" || active === "excited"}
        {@const list = active === "excited" ? pages.excited : pages.electronic}
        <div class="flex flex-wrap items-center gap-3 mb-3 text-xs">
          <div class="flex gap-1 bg-slate-100 p-1 rounded-xl border border-slate-200">
            {#each [["sf", "spin free"], ["soc", "spin–orbit"]] as [s, l]}
              <button class={seg(spin === s)} disabled={!group?.[s]} onclick={() => (spin = s)}>{l}</button>
            {/each}
          </div>
          {#if hasTrapPages}
            <span class="font-bold text-slate-400">Surface traps</span>
            <div class="flex gap-1 bg-slate-100 p-1 rounded-xl border border-slate-200">
              {#each [["included", "included"], ["filtered", "filtered out"]] as [t, l]}
                <button class={seg(traps === t)} onclick={() => (traps = t)}>{l}</button>
              {/each}
            </div>
          {/if}
          {#if list.length > 1}
            <span class="font-bold text-slate-400">Functional</span>
            <div class="flex gap-1 bg-slate-100 p-1 rounded-xl border border-slate-200">
              {#each list as g, i}<button class={seg((list[fn] ? fn : 0) === i)} onclick={() => (fn = i)}>{g.code} {g.functional}{g.chain ? "" : " (earlier run)"}</button>{/each}
            </div>
          {:else if group}
            <span class="text-slate-500">{group.code} {group.functional}</span>
          {/if}
        </div>
        {#if pageUrl}
          {#key pageUrl}
            <iframe title={active === "excited" ? "Excited states" : "Electronic structure"} src={pageUrl}
                    class="w-full h-[85vh] rounded-[1rem] border border-slate-200 bg-white"
                    sandbox="allow-scripts allow-same-origin allow-popups allow-modals" referrerpolicy="no-referrer"></iframe>
          {/key}
        {/if}
      {/if}
    </div>
  </main>
{/if}
