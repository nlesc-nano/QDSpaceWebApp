<script>
  import Viewer from './Viewer.svelte';

  // --- Common Oxidation States for QD Elements & Ligands ---
  const OXIDATION_STATES = {
    'Cd': 2, 'Zn': 2, 'Pb': 2, 'Hg': 2, 'Mg': 2, 'Ca': 2, 'Sr': 2, 'Ba': 2,
    'In': 3, 'Ga': 3, 'Al': 3, 'Bi': 3, 
    'Cs': 1, 'Rb': 1, 'K': 1, 'Na': 1, 'Li': 1, 'Ag': 1, 'Cu': 1, 'Au': 1,
    'S': -2, 'Se': -2, 'Te': -2, 'O': -2,
    'P': -3, 'As': -3, 'N': -3, 'Sb': -3,
    'F': -1, 'Cl': -1, 'Br': -1, 'I': -1
  };

  // ==========================================
  // STATE: Inputs & Configuration
  // ==========================================
  let coreFile = $state(null);
  let radius = $state(20);
  let aspect = $state([1.0, 1.0, 1.0]);
  let activePreset = $state('Sphere');
  let coreFacets = $state([{ id: crypto.randomUUID(), hkl: '100', gamma: 1.0 }]);
  let shells = $state([]);
  
  let passivateExpanded = $state(false);
  let capDist = $state('random');
  let anionicLigands = $state([]);
  let showCationic = $state(false);
  let cationicLigands = $state([]);
  
  let runMode = $state('quiet');
  let posQ = $state('remove');

  // ==========================================
  // STATE: UI & Outputs
  // ==========================================
  let isBuilding = $state(false);
  let logs = $state("[status] Ready.\n");
  let xyzData = $state("");
  let finalResult = $state(null);
  let logContainer;

  let lastUnpassivatedXyz = $state("");
  let tooltipText = $state("");
  let tooltipPos = $state({ x: 0, y: 0 });

  // --- DYNAMIC STOICHIOMETRY & RATIOS PARSER ---
  let calculatedStoich = $derived.by(() => {
    if (!xyzData) return {};
    const lines = xyzData.trim().split('\n');
    const nAtoms = parseInt(lines[0].trim());
    if (isNaN(nAtoms)) return {};
    
    const sto = {};
    for (let i = 2; i < 2 + nAtoms && i < lines.length; i++) {
      const parts = lines[i].trim().split(/\s+/);
      if (parts.length >= 4) {
        const el = parts[0];
        sto[el] = (sto[el] || 0) + 1;
      }
    }
    return sto;
  });

  // Calculate Elemental Ratios
  let calculatedRatios = $derived.by(() => {
    const sto = calculatedStoich;
    const elements = Object.keys(sto);
    const ratios = {};
    
    for (let i = 0; i < elements.length; i++) {
      for (let j = 0; j < elements.length; j++) {
        if (i !== j) {
          const el1 = elements[i];
          const el2 = elements[j];
          ratios[`${el1}/${el2}`] = (sto[el1] / sto[el2]).toFixed(3);
        }
      }
    }
    return ratios;
  });

  // Advanced Charge Calculator
  let calculatedCharge = $derived.by(() => {
    if (!xyzData) return 0;
    let charge = 0;

    // 1. If we have the API response breakdown, it's the safest way to calculate
    if (finalResult?.grouped_counts) {
      const core = finalResult.grouped_counts.core?.by_element || {};
      const shell = finalResult.grouped_counts.shell?.by_element || {};

      // A) Calculate core + shell inorganic charges
      for (const [el, count] of Object.entries(core)) {
         charge += (OXIDATION_STATES[el] || 0) * count;
      }
      for (const [el, count] of Object.entries(shell)) {
         charge += (OXIDATION_STATES[el] || 0) * count;
      }

      // B) Calculate complex organic ligand charges (-1 for anionic, +1 for cationic)
      if (finalResult.ligand_detail) {
         const anions = Object.values(finalResult.ligand_detail.anionic || {}).reduce((a,b)=>a+b, 0);
         const cations = Object.values(finalResult.ligand_detail.cationic || {}).reduce((a,b)=>a+b, 0);
         charge -= anions;
         charge += cations;
      }

      // C) Add leftover/initial Dummy atoms (like Cl, Br, Rb) present in the XYZ
      // This ensures dummy-passivated cores are correctly calculated as 0
      const inorganicElements = { ...core };
      for (const [el, count] of Object.entries(shell)) {
          inorganicElements[el] = (inorganicElements[el] || 0) + count;
      }

      const knownDummies = ['Cl', 'Br', 'I', 'F', 'Na', 'K', 'Rb', 'Cs'];
      for (const dummy of knownDummies) {
          if (calculatedStoich[dummy]) {
              const coreCount = inorganicElements[dummy] || 0;
              const extraDummies = calculatedStoich[dummy] - coreCount; // Only count dummies NOT part of the core
              if (extraDummies > 0) {
                  charge += (OXIDATION_STATES[dummy] || 0) * extraDummies;
              }
          }
      }
      return charge;
    }

    // 2. Fallback if API fails: Ignore purely organic elements
    const qdCoreElements = ['Cd', 'Zn', 'Pb', 'Hg', 'In', 'Ga', 'Al', 'S', 'Se', 'Te', 'P', 'As', 'Sb', 'Ag', 'Cu', 'Au', 'Cs', 'Rb', 'K', 'Na'];
    for (const [el, count] of Object.entries(calculatedStoich)) {
       if (qdCoreElements.includes(el) || ['Cl', 'Br', 'I', 'F'].includes(el)) {
           charge += (OXIDATION_STATES[el] || 0) * count;
       }
    }
    return charge;
  });

  // ==========================================
  // HELPERS & UI LOGIC
  // ==========================================
  function applyPreset(name) {
    activePreset = name;
    const presets = {
      Sphere:   { aspect: [1, 1, 1],       facets: [{ hkl: '100', gamma: 1.0 }] },
      Platelet: { aspect: [1, 1, 0.2],     facets: [{ hkl: '100', gamma: 0.8 }] },
      Rod:      { aspect: [1, 0.2, 0.2],   facets: [{ hkl: '100', gamma: 0.8 }] }
    };
    if (presets[name]) {
      aspect = [...presets[name].aspect];
      coreFacets = presets[name].facets.map(f => ({ id: crypto.randomUUID(), ...f }));
    }
  }

  function handleCoreDrop(e) {
    e.preventDefault();
    if (e.dataTransfer.files && e.dataTransfer.files[0]) {
      coreFile = e.dataTransfer.files[0];
    }
  }

  function addShell() {
    shells.push({
      id: crypto.randomUUID(),
      file: null,
      aspect: [1.0, 1.0, 1.0],
      facets: coreFacets.map(f => ({ id: crypto.randomUUID(), hkl: f.hkl, gamma: f.gamma }))
    });
  }

  $effect(() => {
    if (logs && logContainer) {
      logContainer.scrollTop = logContainer.scrollHeight;
    }
  });

  function downloadXYZ() {
    if (!xyzData) return;
    const name = finalResult?.download_name || 'final.xyz';
    const a = document.createElement('a');
    a.href = URL.createObjectURL(new Blob([xyzData], { type: 'text/plain' }));
    a.download = name;
    a.click();
  }

  function onReset() {
    coreFile = null;
    xyzData = "";
    finalResult = null;
    lastUnpassivatedXyz = "";
    logs = "[status] Reset.\n";
    shells = [];
    anionicLigands = [];
    cationicLigands = [];
    applyPreset('Sphere');
  }

  // ==========================================
  // API BUILD STREAM
  // ==========================================
  async function onBuild(skipCoreBuild = false) {
    if (!coreFile) { alert('Please upload a core .cif file.'); return; }
    if (skipCoreBuild && !lastUnpassivatedXyz) {
      alert('No built structure found. Please run a full build first.');
      return;
    }

    isBuilding = true;
    finalResult = null;
    xyzData = "";
    logs = `[status] Sending ${skipCoreBuild ? 'repassivation' : 'build'} job to /builder/api/build_stream …\n`;

    const formData = new FormData();
    formData.append('files', coreFile);
    for (const shell of shells) {
      if (shell.file) formData.append('files', shell.file);
    }

    const options = {
      radius_A: radius,
      core_cif_filename: coreFile.name,
      aspect: aspect,
      facets: coreFacets.map(f => ({ hkl: f.hkl, gamma: f.gamma })),
      shells: shells.filter(s => s.file).map(s => ({
        material_cif: s.file.name,
        aspect: s.aspect,
        facets: s.facets.map(f => ({ hkl: f.hkl, gamma: f.gamma }))
      })),
      cap_distribution: capDist,
      cap_anionic_jobs: anionicLigands.map(l => ({ smiles: l.smiles, ratio: l.ratio, dummy: l.dummy })),
      cap_cationic_jobs: cationicLigands.map(l => ({ smiles: l.smiles, ratio: l.ratio, dummy: l.dummy })),
      skip_core_build: skipCoreBuild,
      xyz_unpassivated: skipCoreBuild ? lastUnpassivatedXyz : null
    };

    formData.append('options', JSON.stringify(options));
    formData.append('positive_q_mode', posQ);
    formData.append('mode', runMode);

    const IS_LOCAL = window.location.hostname === "localhost" || window.location.hostname === "127.0.0.1";
    const BUILD_STREAM_URL = IS_LOCAL 
      ? "http://127.0.0.1:8000/builder/api/build_stream" 
      : "/builder/api/build_stream";

    try {
      const res = await fetch(BUILD_STREAM_URL, { method: 'POST', body: formData });
      if (!res.ok) {
        const errorText = await res.text();
        throw new Error(`HTTP ${res.status}: ${errorText}`);
      }

      const reader = res.body.getReader();
      const decoder = new TextDecoder();
      let buf = '';

      while (true) {
        const { done, value } = await reader.read();
        if (done) break;
        buf += decoder.decode(value, { stream: true });
        
        let nl;
        while ((nl = buf.indexOf('\n')) >= 0) {
          const line = buf.slice(0, nl).trim();
          buf = buf.slice(nl + 1);
          if (!line) continue;
          
          try {
            const evt = JSON.parse(line);
            if (evt.event === 'log') logs += evt.line + "\n";
            else if (evt.event === 'cmd') logs += `[cmd] ${evt.line}\n`;
            else if (evt.event === 'result') finalResult = evt;
          } catch (e) {
            logs += line + "\n";
          }
        }
      }

      if (finalResult && finalResult.status === 'success') {
        xyzData = finalResult.xyz_passivated || finalResult.xyz || finalResult.xyz_unpassivated || "";
        if (!skipCoreBuild) {
            lastUnpassivatedXyz = finalResult.xyz_unpassivated || finalResult.xyz_passivated || xyzData;
        }
        if (finalResult.last_command) logs += `[cmd][final] ${finalResult.last_command}\n`;
        logs += "\n[status] Rendered.\n";
      } else {
        logs += "[error] Build failed.\n";
      }
    } catch (err) {
      logs += `[error] fetch failed: ${err.message}\n`;
    } finally {
      isBuilding = false;
    }
  }
</script>

<svelte:window onmousemove={(e) => tooltipPos = { x: e.pageX + 15, y: e.pageY + 15 }} />
{#if tooltipText}
  <div class="fixed bg-slate-900/90 backdrop-blur-sm text-white px-3 py-1.5 rounded-lg text-sm z-50 pointer-events-none shadow-xl font-medium" style="left: {tooltipPos.x}px; top: {tooltipPos.y}px;">
    {tooltipText}
  </div>
{/if}

<div class="flex flex-col lg:flex-row h-[calc(100vh-64px)] overflow-hidden font-sans bg-slate-50">
  
  <aside class="w-full lg:w-[380px] p-6 bg-slate-50 overflow-y-auto flex-shrink-0 border-r border-slate-200 flex flex-col gap-6">
    
    <div class="px-2 pt-2 mb-2">
      <h1 class="font-heading text-2xl font-bold flex items-center gap-3 mb-2 text-slate-900 tracking-tight">
        <img src="/assets/logos/QD_Builder_mini.png" alt="QD Builder" class="h-9 w-auto mix-blend-multiply">
        Quantum Dot Builder
      </h1>
      <p class="text-sm text-slate-600 font-medium leading-relaxed">
        Upload a .cif, set parameters, and build your nanocrystal.
      </p>
    </div>

    <div class="bg-white border border-slate-100 shadow-sm rounded-[1.5rem] p-6">
      <h2 class="font-heading text-lg font-bold text-slate-900 mb-4">1) Core Structure</h2>
      <label class="border-2 border-dashed border-slate-200 hover:bg-brand-50 hover:border-brand-400 hover:text-brand-700 transition-all p-6 rounded-2xl text-center cursor-pointer block group"
             ondragover={(e) => e.preventDefault()}
             ondrop={handleCoreDrop}>
        <input type="file" class="hidden" accept=".cif" onchange={(e) => coreFile = e.target.files[0]}>
        <div class="mt-2 text-sm text-slate-500 font-medium group-hover:text-brand-600 transition-colors">
          {#if coreFile}
            <span class="font-bold text-brand-600">{coreFile.name}</span>
          {:else}
            <span class="font-bold text-brand-600">Click to upload</span> or drag & drop a .cif file
          {/if}
        </div>
      </label>
    </div>

    <div class="bg-white border border-slate-100 shadow-sm rounded-[1.5rem] p-6">
      <h2 class="font-heading text-lg font-bold text-slate-900 mb-4">2) Size & Aspect</h2>
      <label class="block text-xs font-bold text-slate-400 uppercase tracking-widest mb-2" for="radius_input">Final Radius (Å)</label>
      <input id="radius_input" type="number" bind:value={radius} step="1" min="1" class="w-full mb-5 border-none ring-1 ring-slate-200 rounded-xl px-4 py-2.5 text-sm bg-slate-50 focus:ring-2 focus:ring-brand-400 outline-none font-medium">
      
      <span class="block text-xs font-bold text-slate-400 uppercase tracking-widest mb-2">Core Aspect Ratio</span>
      <div class="flex space-x-1 rounded-xl p-1 bg-slate-100 my-2">
        {#each ['Sphere', 'Platelet', 'Rod', 'Custom'] as preset}
          <button class="w-full rounded-lg py-1.5 text-sm font-bold transition-all {activePreset === preset ? 'bg-white text-brand-600 shadow-sm' : 'text-slate-500 hover:text-slate-700 hover:bg-slate-200/50'}"
                  onclick={() => applyPreset(preset)}>{preset}</button>
        {/each}
      </div>
      <div class="grid grid-cols-3 gap-3 mt-3">
        <input type="number" bind:value={aspect[0]} step="0.1" class="w-full text-center border-none ring-1 ring-slate-200 rounded-xl py-2 bg-slate-50 text-sm font-medium focus:ring-2 focus:ring-brand-400 outline-none" oninput={() => activePreset = 'Custom'}>
        <input type="number" bind:value={aspect[1]} step="0.1" class="w-full text-center border-none ring-1 ring-slate-200 rounded-xl py-2 bg-slate-50 text-sm font-medium focus:ring-2 focus:ring-brand-400 outline-none" oninput={() => activePreset = 'Custom'}>
        <input type="number" bind:value={aspect[2]} step="0.1" class="w-full text-center border-none ring-1 ring-slate-200 rounded-xl py-2 bg-slate-50 text-sm font-medium focus:ring-2 focus:ring-brand-400 outline-none" oninput={() => activePreset = 'Custom'}>
      </div>
    </div>

    <div class="bg-white border border-slate-100 shadow-sm rounded-[1.5rem] p-6">
      <h2 class="font-heading text-lg font-bold text-slate-900 mb-4">3) Core Facets</h2>
      <div class="grid grid-cols-[1fr_1fr_32px] gap-3 items-center text-[10px] text-slate-400 uppercase font-extrabold tracking-widest mb-3 px-1">
        <span>Miller Index</span> <span>Energy (γ)</span> <span></span>
      </div>
      {#each coreFacets as facet, i (facet.id)}
        <div class="grid grid-cols-[1fr_1fr_32px] gap-3 mb-3">
          <input type="text" bind:value={facet.hkl} class="border-none ring-1 ring-slate-200 rounded-xl px-2 py-2 text-sm text-center bg-slate-50 focus:ring-2 focus:ring-brand-400 outline-none font-medium min-w-0">
          <input type="number" bind:value={facet.gamma} step="0.1" class="border-none ring-1 ring-slate-200 rounded-xl px-2 py-2 text-sm text-center bg-slate-50 focus:ring-2 focus:ring-brand-400 outline-none font-medium min-w-0">
          <button class="bg-red-50 text-red-600 hover:bg-red-100 hover:text-red-700 rounded-xl flex items-center justify-center font-bold transition-colors" onclick={() => coreFacets.splice(i, 1)}>✕</button>
        </div>
      {/each}
      <button class="w-full bg-slate-50 hover:bg-slate-100 text-brand-600 py-2.5 rounded-xl mt-2 text-sm font-bold border border-slate-200 transition-colors" onclick={() => coreFacets.push({ id: crypto.randomUUID(), hkl: '100', gamma: 1.0 })}>+ Add Facet</button>
    </div>

    <div class="bg-white border border-slate-100 shadow-sm rounded-[1.5rem] p-6">
      <h2 class="font-heading text-lg font-bold text-slate-900 mb-4">4) Shells</h2>
      {#each shells as shell, sIndex (shell.id)}
        <div class="border border-brand-200 bg-brand-50/40 p-4 rounded-[1.25rem] mb-4">
          <div class="flex justify-between items-center mb-4">
            <span class="font-bold text-brand-800 text-sm bg-brand-100 px-3 py-1 rounded-lg">Shell {sIndex + 1}</span>
            <button class="text-red-500 hover:text-red-700 text-xs font-bold transition-colors" onclick={() => shells.splice(sIndex, 1)}>Remove</button>
          </div>
          
          <span class="block text-[10px] font-extrabold text-brand-700 uppercase tracking-widest mb-2">Material CIF</span>
          <input type="file" accept=".cif" class="w-full text-xs mb-4 file:mr-4 file:py-1.5 file:px-4 file:rounded-full file:border-0 file:text-xs file:font-bold file:bg-brand-100 file:text-brand-700 hover:file:bg-brand-200 transition-all text-slate-600" onchange={(e) => shell.file = e.target.files[0]}>
          
          <span class="block text-[10px] font-extrabold text-brand-700 uppercase tracking-widest mb-2">Aspect Ratio</span>
          <div class="grid grid-cols-3 gap-2 mb-4">
            <input type="number" bind:value={shell.aspect[0]} step="0.1" class="border-none ring-1 ring-brand-200 rounded-lg px-2 py-1.5 text-xs text-center bg-white focus:ring-2 focus:ring-brand-400 outline-none font-medium min-w-0">
            <input type="number" bind:value={shell.aspect[1]} step="0.1" class="border-none ring-1 ring-brand-200 rounded-lg px-2 py-1.5 text-xs text-center bg-white focus:ring-2 focus:ring-brand-400 outline-none font-medium min-w-0">
            <input type="number" bind:value={shell.aspect[2]} step="0.1" class="border-none ring-1 ring-brand-200 rounded-lg px-2 py-1.5 text-xs text-center bg-white focus:ring-2 focus:ring-brand-400 outline-none font-medium min-w-0">
          </div>
          
          <span class="block text-[10px] font-extrabold text-brand-700 uppercase tracking-widest mb-2">Facets</span>
          {#each shell.facets as sfacet, fi (sfacet.id)}
            <div class="grid grid-cols-[1fr_1fr_32px] gap-2 mb-2">
              <input type="text" bind:value={sfacet.hkl} class="border-none ring-1 ring-brand-200 rounded-lg px-2 py-1.5 text-xs text-center bg-white focus:ring-2 focus:ring-brand-400 outline-none font-medium min-w-0">
              <input type="number" bind:value={sfacet.gamma} step="0.1" class="border-none ring-1 ring-brand-200 rounded-lg px-2 py-1.5 text-xs text-center bg-white focus:ring-2 focus:ring-brand-400 outline-none font-medium min-w-0">
              <button class="bg-red-100 text-red-600 hover:bg-red-200 rounded-lg flex items-center justify-center text-xs font-bold transition-colors" onclick={() => shell.facets.splice(fi, 1)}>✕</button>
            </div>
          {/each}
          <button class="w-full bg-white hover:bg-brand-50 border border-brand-200 text-brand-700 text-xs py-2 mt-2 rounded-lg font-bold transition-colors" onclick={() => shell.facets.push({ id: crypto.randomUUID(), hkl: '100', gamma: 1.0 })}>+ Add Shell Facet</button>
        </div>
      {/each}
      <button class="w-full bg-slate-50 hover:bg-slate-100 border border-slate-200 text-brand-600 py-2.5 rounded-xl text-sm font-bold transition-colors" onclick={addShell}>+ Add Shell</button>
    </div>

    <div class="bg-white border border-slate-100 shadow-sm rounded-[1.5rem] p-6">
      <h2 class="font-heading text-lg font-bold text-slate-900 mb-4">5) Passivation</h2>
      {#if !passivateExpanded}
        <button class="w-full bg-slate-50 hover:bg-slate-100 border border-slate-200 text-accent-600 py-2.5 rounded-xl text-sm font-bold transition-colors" onclick={() => passivateExpanded = true}>+ Add Ligands</button>
      {:else}
        <div class="border border-accent-200 p-5 rounded-[1.25rem] bg-accent-50/30">
          <span class="block text-[10px] font-extrabold text-accent-700 uppercase tracking-widest mb-3">Anionic Ligands</span>
          <div class="grid grid-cols-[1fr_3.5rem_3.5rem_28px] gap-2 text-[10px] text-accent-600 uppercase font-bold tracking-widest mb-2 px-1">
            <span>SMILES</span><span>Dummy</span><span>Ratio</span><span></span>
          </div>
          {#each anionicLigands as lig, i (lig.id)}
            <div class="grid grid-cols-[1fr_3.5rem_3.5rem_28px] gap-2 mb-3">
              <input type="text" bind:value={lig.smiles} class="border-none ring-1 ring-accent-200 rounded-lg px-2 py-2 text-sm bg-white focus:ring-2 focus:ring-accent-400 outline-none font-medium min-w-0">
              <input type="text" bind:value={lig.dummy} class="border-none ring-1 ring-accent-200 rounded-lg px-1 py-2 text-sm text-center bg-white focus:ring-2 focus:ring-accent-400 outline-none font-medium min-w-0">
              <input type="number" bind:value={lig.ratio} step="0.1" class="border-none ring-1 ring-accent-200 rounded-lg px-1 py-2 text-sm text-center bg-white focus:ring-2 focus:ring-accent-400 outline-none font-medium min-w-0">
              <button class="bg-red-100 text-red-600 hover:bg-red-200 rounded-lg flex items-center justify-center text-xs font-bold transition-colors" onclick={() => anionicLigands.splice(i, 1)}>✕</button>
            </div>
          {/each}
          <button class="w-full bg-white hover:bg-accent-50 border border-accent-200 text-xs py-2 rounded-lg text-accent-700 mb-6 font-bold transition-colors" onclick={() => anionicLigands.push({ id: crypto.randomUUID(), smiles: '', ratio: 1.0, dummy: 'Cl' })}>+ Add Anion</button>

          <span class="block text-[10px] font-extrabold text-accent-700 uppercase tracking-widest mb-3">Distribution</span>
          <div class="flex gap-4 text-sm mb-6 bg-white p-2 rounded-xl border border-accent-100">
            <label class="flex items-center gap-1.5 cursor-pointer"><input type="radio" bind:group={capDist} value="uniform" class="accent-accent-600"> <span class="font-medium text-slate-700">Uniform</span></label>
            <label class="flex items-center gap-1.5 cursor-pointer"><input type="radio" bind:group={capDist} value="segmented" class="accent-accent-600"> <span class="font-medium text-slate-700">Segmented</span></label>
            <label class="flex items-center gap-1.5 cursor-pointer"><input type="radio" bind:group={capDist} value="random" class="accent-accent-600"> <span class="font-medium text-slate-700">Random</span></label>
          </div>

          {#if !showCationic}
            <button class="w-full bg-white hover:bg-accent-50 border border-accent-200 text-xs py-2 rounded-lg text-accent-700 font-bold transition-colors" onclick={() => showCationic = true}>+ Add Cationic Ligand</button>
          {:else}
            <span class="block text-[10px] font-extrabold text-accent-700 uppercase tracking-widest mb-3 pt-4 border-t border-accent-200">Cationic Ligands</span>
            {#each cationicLigands as lig, i (lig.id)}
              <div class="grid grid-cols-[1fr_3.5rem_3.5rem_28px] gap-2 mb-3">
                <input type="text" bind:value={lig.smiles} class="border-none ring-1 ring-accent-200 rounded-lg px-2 py-2 text-sm bg-white focus:ring-2 focus:ring-accent-400 outline-none font-medium min-w-0">
                <input type="text" bind:value={lig.dummy} class="border-none ring-1 ring-accent-200 rounded-lg px-1 py-2 text-sm text-center bg-white focus:ring-2 focus:ring-accent-400 outline-none font-medium min-w-0">
                <input type="number" bind:value={lig.ratio} step="0.1" class="border-none ring-1 ring-accent-200 rounded-lg px-1 py-2 text-sm text-center bg-white focus:ring-2 focus:ring-accent-400 outline-none font-medium min-w-0">
                <button class="bg-red-100 text-red-600 hover:bg-red-200 rounded-lg flex items-center justify-center text-xs font-bold transition-colors" onclick={() => cationicLigands.splice(i, 1)}>✕</button>
              </div>
            {/each}
            <button class="w-full bg-white hover:bg-accent-50 border border-accent-200 text-xs py-2 rounded-lg text-accent-700 font-bold transition-colors" onclick={() => cationicLigands.push({ id: crypto.randomUUID(), smiles: '', ratio: 1.0, dummy: 'Rb' })}>+ Add Cation</button>
          {/if}
        </div>
      {/if}

      {#if lastUnpassivatedXyz}
        <hr class="my-5 border-slate-100" />
        <button class="w-full bg-accent-600 hover:bg-accent-700 text-white font-bold py-3 rounded-xl text-sm transition-all shadow-glow active:scale-95 flex items-center justify-center gap-2" 
                onclick={() => onBuild(true)} disabled={isBuilding}>
          {isBuilding ? 'Processing...' : '⚡ Repassivate Current Structure'}
        </button>
      {/if}
    </div>

    <div class="bg-white border border-slate-100 shadow-sm rounded-[1.5rem] p-6 mb-8">
      <h2 class="font-heading text-lg font-bold text-slate-900 mb-5">6) Build Settings</h2>
      <div class="mb-5">
        <span class="block text-xs font-bold text-slate-400 uppercase tracking-widest mb-3">Run Mode</span>
        <div class="flex gap-4">
          <label class="flex items-center gap-1.5 cursor-pointer text-sm font-medium text-slate-700"><input type="radio" bind:group={runMode} value="quiet" class="accent-brand-600"> Turbo (quiet)</label>
          <label class="flex items-center gap-1.5 cursor-pointer text-sm font-medium text-slate-700"><input type="radio" bind:group={runMode} value="live-tty" class="accent-brand-600"> Live logs</label>
        </div>
      </div>
      <div class="mb-6">
        <span class="block text-xs font-bold text-slate-400 uppercase tracking-widest mb-3">Positive-Q Strategy</span>
        <div class="flex gap-4">
          <label class="flex items-center gap-1.5 cursor-pointer text-sm font-medium text-slate-700"><input type="radio" bind:group={posQ} value="remove" class="accent-brand-600"> Remove (default)</label>
          <label class="flex items-center gap-1.5 cursor-pointer text-sm font-medium text-slate-700"><input type="radio" bind:group={posQ} value="add" class="accent-brand-600"> Add</label>
        </div>
      </div>
      <div class="flex gap-3">
        <button class="flex-1 bg-brand-600 hover:bg-brand-700 text-white font-bold py-3 rounded-xl shadow-glow transition-all active:scale-95" 
                onclick={() => onBuild(false)} disabled={isBuilding}>
          {isBuilding ? 'Building...' : 'Build Nanocrystal'}
        </button>
        <button class="bg-slate-100 hover:bg-slate-200 text-slate-700 font-bold py-3 px-6 rounded-xl transition-all" onclick={onReset}>Reset</button>
      </div>
      {#if finalResult}
        <button class="w-full bg-emerald-600 hover:bg-emerald-700 text-white font-bold py-3 rounded-xl mt-4 shadow-sm transition-all active:scale-95" onclick={downloadXYZ}>
          Download {finalResult.download_name || 'final.xyz'}
        </button>
      {/if}
    </div>

  </aside>
     
  <main class="flex-grow h-full flex flex-col p-6 gap-6 bg-slate-50 overflow-y-auto">
    
    <div class="flex flex-col gap-6 flex-shrink-0 h-[calc(100vh-112px)] min-h-[700px]">
      
      <div class="relative flex-1 min-h-[400px] bg-white rounded-[1.5rem] p-4 border border-slate-100 shadow-sm flex flex-col">
        <div class="flex justify-between items-center mb-3 px-2">
          <h2 class="font-heading font-bold text-xl text-slate-900">Live Render</h2>
        </div>
        <div class="flex-1 bg-slate-50 rounded-[1rem] border border-slate-200 overflow-hidden relative shadow-inner">
          <Viewer xyz={xyzData} />
          {#if isBuilding}
            <div class="absolute inset-0 bg-slate-900/60 backdrop-blur-sm flex flex-col items-center justify-center text-white z-10">
              <div class="w-12 h-12 border-4 border-brand-500 border-t-transparent rounded-full animate-spin mb-4"></div>
              <p class="font-bold text-lg tracking-wide">Processing Architecture...</p>
            </div>
          {/if}
        </div>
      </div>
      
      <div class="h-1/3 grid grid-cols-1 md:grid-cols-4 gap-6 min-h-[280px]">
        
        <div class="bg-white border border-slate-100 shadow-sm rounded-[1.5rem] p-6 md:col-span-1 overflow-y-auto h-full">
          <h2 class="font-heading text-lg font-bold text-slate-900 mb-4 border-b border-slate-100 pb-2">Stoichiometry</h2>

          {#if Object.keys(calculatedStoich).length > 0}
            
            <h3 class="text-xs font-bold text-slate-400 uppercase tracking-widest mb-3">Total Stoichiometry</h3>
            <div class="flex flex-wrap gap-2 mb-6">
              {#each Object.entries(calculatedStoich) as [elem, count]}
                <span class="bg-brand-50 text-brand-800 border border-brand-100 px-3 py-1.5 rounded-xl font-bold text-sm">
                  {elem}: {count}
                </span>
              {/each}
            </div>

            <h3 class="text-xs font-bold text-slate-400 uppercase tracking-widest mb-3">Elemental Ratios</h3>
            <div class="flex flex-wrap gap-2 mb-6">
              {#each Object.entries(calculatedRatios) as [ratio, val]}
                <span class="bg-accent-50 text-accent-800 border border-accent-100 px-2.5 py-1 rounded-lg font-mono text-xs font-bold">
                  <span class="text-accent-600">{ratio}</span>: {val}
                </span>
              {/each}
            </div>

            <h3 class="text-xs font-bold text-slate-400 uppercase tracking-widest mb-3">Charge Balance</h3>
            <div class="mb-6">
              <div class="flex items-center justify-between px-3 py-2.5 rounded-lg border {calculatedCharge === 0 ? 'bg-emerald-50 border-emerald-200 text-emerald-800' : 'bg-red-50 border-red-200 text-red-800'} font-bold text-xs shadow-sm w-full">
                <span>Total Charge:</span>
                <span class="text-sm">{calculatedCharge > 0 ? '+' : ''}{calculatedCharge} e</span>
              </div>
              {#if calculatedCharge !== 0}
                <p class="text-[10px] text-red-600 mt-1.5 font-bold flex items-center gap-1">
                  <svg class="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z"></path></svg>
                  Structure is not charge-neutral.
                </p>
              {/if}
            </div>

          {/if}

          {#if finalResult?.grouped_counts}
            <h3 class="text-xs font-bold text-slate-400 uppercase tracking-widest mb-3 border-t border-slate-100 pt-4">System Breakdown</h3>
            <div class="text-sm">
              <p class="font-extrabold text-slate-900 text-base mb-3">Total Atoms: {finalResult.grouped_counts.total_atoms}</p>
              
              <p class="text-[10px] font-bold text-slate-400 uppercase tracking-widest mb-2">Core ({finalResult.grouped_counts.core.total})</p>
              <div class="flex flex-wrap gap-2 mb-4">
                {#each Object.entries(finalResult.grouped_counts.core.by_element || {}) as [el, count]}
                  <span class="bg-indigo-50 text-indigo-800 border border-indigo-200 px-2 py-0.5 rounded font-medium text-xs">{el}: {count}</span>
                {/each}
              </div>

              {#if finalResult.ligand_detail && finalResult.ligand_detail.total > 0}
                <p class="text-[10px] font-bold text-slate-400 uppercase tracking-widest mb-2">Ligands ({finalResult.ligand_detail.total})</p>
                {#each Object.entries(finalResult.ligand_detail.anionic || {}) as [smiles, count]}
                  {#if count > 0}<div class="flex justify-between text-xs font-mono bg-slate-50 border border-slate-100 p-2 rounded-lg mb-1.5"><span class="truncate pr-2 text-slate-600">{smiles}</span><span class="font-bold text-slate-900">{count}</span></div>{/if}
                {/each}
                {#each Object.entries(finalResult.ligand_detail.cationic || {}) as [smiles, count]}
                  {#if count > 0}<div class="flex justify-between text-xs font-mono bg-slate-50 border border-slate-100 p-2 rounded-lg mb-1.5"><span class="truncate pr-2 text-slate-600">{smiles}</span><span class="font-bold text-slate-900">{count}</span></div>{/if}
                {/each}
              {:else if finalResult.grouped_counts.shell?.total > 0}
                 <p class="text-[10px] font-bold text-slate-400 uppercase tracking-widest mb-2">Shell ({finalResult.grouped_counts.shell.total})</p>
                 <div class="flex flex-wrap gap-2 mb-4">
                   {#each Object.entries(finalResult.grouped_counts.shell.by_element || {}) as [el, count]}
                     <span class="bg-blue-50 text-blue-800 border border-blue-200 px-2 py-0.5 rounded font-medium text-xs">{el}: {count}</span>
                   {/each}
                 </div>
              {/if}
            </div>
          {:else if !xyzData}
            <div class="h-32 flex items-center justify-center">
              <p class="text-sm text-slate-400 italic">No structure built yet.</p>
            </div>
          {/if}
        </div>

        <div class="bg-slate-900 border border-slate-800 rounded-[1.5rem] p-0 md:col-span-3 flex flex-col overflow-hidden shadow-inner relative h-full">
          <div class="bg-slate-800/80 backdrop-blur-sm text-slate-400 text-[10px] px-4 py-2 font-mono font-bold uppercase tracking-widest border-b border-slate-700/50 flex justify-between items-center absolute w-full top-0 left-0 z-10">
            <span class="flex items-center gap-2"><div class="w-2 h-2 rounded-full bg-emerald-500"></div> miniCAT Engine Terminal</span>
            {#if isBuilding}<span class="text-brand-400 animate-pulse">Running process...</span>{/if}
          </div>
          <pre bind:this={logContainer} class="text-xs text-emerald-400 font-mono p-5 pt-12 h-full overflow-y-auto whitespace-pre-wrap leading-relaxed selection:bg-emerald-900 selection:text-emerald-100">{logs}</pre>
        </div>

      </div>
    </div>
  </main>
</div>

