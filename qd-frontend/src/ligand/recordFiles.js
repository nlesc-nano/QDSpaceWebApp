// Where a library record's computed data lives, for the dot page and the library preview.
// Builder records keep the MACE-MH-1 results in <record>/props/ and the DFT chain in <record>/dft/
// (import_dft_results.py); older DFT runs come from the legacy tree, merged by ingest_library.py.
// The dev server answers a missing file with index.html, so JSON is recognised by its content type.

export const recordDir = (r) => {
  const start = (r?.stages || []).find((s) => s.stage === "start" && !s.dft_start)?.file;
  return start ? `/${start.slice(0, start.lastIndexOf("/"))}` : null;
};

const isJson = (u) => (u ? fetch(u, { method: "HEAD", cache: "no-store" })
  .then((res) => res.ok && (res.headers.get("content-type") || "").includes("json")).catch(() => false) : Promise.resolve(false));
const getJson = (u) => (u ? fetch(u, { cache: "no-store" }).then((res) =>
  res.ok && (res.headers.get("content-type") || "").includes("json") ? res.json() : null).catch(() => null) : Promise.resolve(null));

// props/ files: the ladder, stability, DFT solvation cubes; how IR/Raman and solvation were computed
export async function findFiles(r) {
  const props = r?.properties || {}, dir = recordDir(r);
  const url = (n) => (props[n] ? `/${props[n]}` : dir ? `${dir}/props/${n}` : null);
  const names = ["relax.json", "structure.json", "hessian.json", "vibspec.json", "spectra.json", "report.json"];
  const [det, stab, cubes, field, solv, ...more] = await Promise.all([isJson(url("detachment.json")), isJson(url("stability.json")),
    isJson(url("solvation_cubes.json")), getJson(url("field.json")), getJson(url("solvation.json")), ...names.map((n) => isJson(url(n)))]);
  const model = solv?.summary?.model || "";
  const found = Object.fromEntries(names.map((n, i) => [n.replace(".json", ""), more[i] ? url(n) : null]));
  return {
    ...found,                     // relax, structure, hessian, vibspec, spectra, report
    det: det ? url("detachment.json") : null, stab: stab ? url("stability.json") : null,
    cubes: cubes ? url("solvation_cubes.json") : null,
    groundState: props["ground_state.html"] ? `/${props["ground_state.html"]}` : null,
    methods: { vib: field?.summary?.source?.startsWith("CP2K") ? "dft" : field || solv ? "xtb" : null,
               solv: model.startsWith("COSMO") ? "cosmo" : model ? "gb" : null },
  };
}

// DFT pages grouped by functional: {electronic: [{functional, code, sf, soc}], excited: [...]}.
// Pages of the published chain (<record>/dft/) take the chain's method; legacy pages the functional of
// the legacy relaxation in the same folder tree.
export function dftPages(r) {
  const props = r?.properties || {};
  const chain = r?.dft?.method || {};
  const stages = (r?.stages || []).filter((s) => s.stage === "geo_opt" && s.code && s.code !== "MACE");
  const methodOf = (path) => {
    if (path.includes("/dft/")) return { functional: chain.functional || "DFT", code: chain.code || "CP2K", chain: true };
    const top = path.split("/").slice(0, 4).join("/");
    const st = stages.find((s) => !s.file.endsWith("/dft/relaxed.xyz") && s.file.startsWith(top))
            || stages.find((s) => !s.file.endsWith("/dft/relaxed.xyz"));
    return { functional: st?.functional || "earlier DFT", code: st?.code || "CP2K", chain: false };
  };
  const group = (kinds) => {
    const out = new Map();
    for (const [kind, names] of Object.entries(kinds)) {
      for (const n of names) {
        const p = props[n];
        if (!p) continue;
        const m = methodOf(p), key = `${m.functional}|${m.chain}`;
        if (!out.has(key)) out.set(key, { ...m, sf: null, soc: null, sf_trap: null, soc_trap: null });
        const g = out.get(key);
        if (!g[kind]) g[kind] = `/${p}`;
      }
    }
    return [...out.values()].sort((a, b) => b.chain - a.chain);     // the published chain first
  };
  return {
    electronic: group({ sf: ["fuzzy_dashboard_sf.html", "plot.html", "plot.html.gz"], soc: ["fuzzy_dashboard_soc.html"] }),
    // sf_trap / soc_trap: the same run with the trap filter on (transitions with surface trap states excluded)
    excited: group({ sf: ["exciton_analysis_sf.html"], soc: ["exciton_analysis_soc.html"],
                     sf_trap: ["exciton_analysis_sf_traps_filtered.html"], soc_trap: ["exciton_analysis_soc_traps_filtered.html"] }),
  };
}

// Geometries of a record: the start, the MACE-MH-1 relaxation and every DFT relaxation (the published
// chain first), as {label, short, url, chain}
export function geometries(r) {
  const st = r?.stages || [];
  const start = st.find((s) => s.stage === "start" && !s.dft_start);
  const mace = st.find((s) => s.stage === "geo_opt" && s.code === "MACE");
  const dft = st.filter((s) => s.stage === "geo_opt" && s.code && s.code !== "MACE")
    .map((s) => ({ label: `DFT ${s.functional || "(functional unknown)"}`, short: s.functional || "DFT", url: `/${s.file}`,
                   chain: s.file.endsWith("/dft/relaxed.xyz") }))
    .sort((a, b) => b.chain - a.chain);
  return { start: start ? `/${start.file}` : null, mace: mace ? `/${mace.file}` : null, dft };
}
