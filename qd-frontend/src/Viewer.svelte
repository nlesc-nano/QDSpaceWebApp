<script>
  import { Structure } from "matterviz";
  import { Trajectory } from "matterviz/trajectory";
  import { parse_xyz } from "matterviz/structure/parse";
  import { onMount } from "svelte";

  let props = $props();

  // Parse the structure synchronously (used for non-MD single-frame mode)
  let structure = $derived(props.xyz && !props.isMD ? parse_xyz(props.xyz) : null);

  // Only disable bonds if the structure is truly massive (>5000 atoms) to protect the browser
  let isMassive = $derived(structure ? structure.num_atoms > 5000 : false);

  let activeViewer = $state(props.isMD ? "matterviz" : "3dmol");
  
  let container3dmol = $state(null);
  let containerNgl = $state(null);
  let containerMolstar = $state(null);
  
  let viewer3dmol = null;
  let stageNgl = null;
  let viewerMolstar = null;
  
  let scriptLoaded3dmol = $state(false);
  let scriptLoadedNgl = $state(false);
  let scriptLoadedMolstar = $state(false);

  // Helper to load external scripts asynchronously
  function loadScript(src, checkGlobal) {
    return new Promise((resolve) => {
      if (window[checkGlobal]) {
        resolve(true);
        return;
      }
      const script = document.createElement("script");
      script.src = src;
      script.async = true;
      script.onload = () => resolve(true);
      script.onerror = () => resolve(false);
      document.head.appendChild(script);
    });
  }

  // Render 3Dmol structure
  function render3dmol() {
    if (!scriptLoaded3dmol || !container3dmol || !props.xyz || props.isMD) return;
    container3dmol.innerHTML = "";

    if (window.$3Dmol) {
      viewer3dmol = window.$3Dmol.createViewer(container3dmol, {
        backgroundColor: "#0f172a" // Slate-900 matching terminal
      });

      // Add Model
      viewer3dmol.addModel(props.xyz, "xyz");

      // Styles based on system size
      const atomCount = props.xyz.split("\n")[0] ? parseInt(props.xyz.split("\n")[0]) : 0;
      if (atomCount > 4000) {
        viewer3dmol.setStyle({}, { line: {} });
      } else {
        viewer3dmol.setStyle({}, {
          sphere: { scale: 0.28, colorscheme: "Jmol" },
          stick: { radius: 0.1, colorscheme: "Jmol" }
        });
      }

      viewer3dmol.zoomTo();
      viewer3dmol.render();
    }
  }

  // Render NGL Viewer structure
  function renderNgl() {
    if (!scriptLoadedNgl || !containerNgl || !props.xyz || props.isMD) return;
    containerNgl.innerHTML = "";

    if (window.NGL) {
      stageNgl = new window.NGL.Stage(containerNgl, { backgroundColor: "#0f172a" });
      const blob = new Blob([props.xyz], { type: "text/plain" });
      const url = URL.createObjectURL(blob);

      stageNgl.loadFile(url, { ext: "xyz" }).then(function (o) {
        const atomCount = props.xyz.split("\n")[0] ? parseInt(props.xyz.split("\n")[0]) : 0;
        if (atomCount > 4000) {
          o.addRepresentation("line");
        } else {
          o.addRepresentation("ball+stick", { multipleBond: "off", scale: 2.0 });
        }
        stageNgl.autoView();
        URL.revokeObjectURL(url);
      });
    }
  }

  // Render Molstar structure
  function renderMolstar() {
    if (!scriptLoadedMolstar || !containerMolstar || !props.xyz || props.isMD) return;
    containerMolstar.innerHTML = "";

    if (window.molstar) {
      window.molstar.Viewer.create(containerMolstar, {
        layoutIsExpanded: false,
        layoutShowControls: false,
        layoutShowRemoteState: false,
        layoutShowSequence: false,
        layoutShowLog: false,
        viewportShowExpand: false,
        viewportShowSelectionMode: false,
        viewportShowAnimation: false,
        collapseLeftControls: true,
      }).then(v => {
        viewerMolstar = v;
        const blob = new Blob([props.xyz], { type: "text/plain" });
        const url = URL.createObjectURL(blob);
        viewerMolstar.loadAll({ url: url, format: 'xyz' }).then(() => {
          URL.revokeObjectURL(url);
        });
      });
    }
  }

  // React to changes in loading status, xyz data, and activeViewer choice
  $effect(() => {
    if (activeViewer === "3dmol" && scriptLoaded3dmol && props.xyz && container3dmol) {
      render3dmol();
    }
  });

  $effect(() => {
    if (activeViewer === "ngl" && scriptLoadedNgl && props.xyz && containerNgl) {
      renderNgl();
    }
  });

  $effect(() => {
    if (activeViewer === "molstar" && scriptLoadedMolstar && props.xyz && containerMolstar) {
      renderMolstar();
    }
  });

  onMount(() => {
    // Parallel scripts load
    loadScript("https://cdnjs.cloudflare.com/ajax/libs/3dmol/2.2.0/3Dmol-min.js", "$3Dmol").then(ok => {
      if (ok) scriptLoaded3dmol = true;
    });

    loadScript("https://cdnjs.cloudflare.com/ajax/libs/ngl/2.0.0-dev.37/ngl.js", "NGL").then(ok => {
      if (ok) scriptLoadedNgl = true;
    });

    loadScript("https://cdn.jsdelivr.net/npm/molstar@4.3.0/build/viewer/molstar.js", "molstar").then(ok => {
      if (ok) scriptLoadedMolstar = true;
    });

    return () => {
      if (viewer3dmol) viewer3dmol.clear();
      if (stageNgl) stageNgl.dispose();
    };
  });
</script>

<svelte:head>
  <link rel="stylesheet" type="text/css" href="https://cdn.jsdelivr.net/npm/molstar@4.3.0/build/viewer/molstar.css" />
</svelte:head>

<div class="relative w-full h-full">

  {#if !props.isMD && props.xyz}
  <!-- Interactive Selector Toolbar -->
  <div class="absolute top-4 right-4 z-20 flex gap-1 bg-slate-900/80 backdrop-blur-md p-1 rounded-lg border border-slate-700/50">
    <button class="px-2.5 py-1 rounded text-[10px] font-bold transition-all {activeViewer === '3dmol' ? 'bg-brand-600 text-white shadow-sm' : 'text-slate-400 hover:text-slate-200'}"
            onclick={() => activeViewer = '3dmol'}>
      3Dmol
    </button>
    <button class="px-2.5 py-1 rounded text-[10px] font-bold transition-all {activeViewer === 'ngl' ? 'bg-brand-600 text-white shadow-sm' : 'text-slate-400 hover:text-slate-200'}"
            onclick={() => activeViewer = 'ngl'}>
      NGL
    </button>
    <button class="px-2.5 py-1 rounded text-[10px] font-bold transition-all {activeViewer === 'molstar' ? 'bg-brand-600 text-white shadow-sm' : 'text-slate-400 hover:text-slate-200'}"
            onclick={() => activeViewer = 'molstar'}>
      Mol*
    </button>
    <button class="px-2.5 py-1 rounded text-[10px] font-bold transition-all {activeViewer === 'matterviz' ? 'bg-brand-600 text-white shadow-sm' : 'text-slate-400 hover:text-slate-200'}"
            onclick={() => activeViewer = 'matterviz'}>
      MatterViz
    </button>
  </div>
  {/if}

  {#if !props.xyz}
    <div class="p-4 flex items-center justify-center h-full text-slate-500 font-medium bg-slate-900 rounded-[1.5rem]">
      No structure loaded
    </div>
  {:else if props.isMD && props.dataUrl}
    {#key props.dataUrl}
      <Trajectory
        data_url={props.dataUrl}
        auto_play={true}
        fps={10}
        structure_props={{
          performance_mode: "speed",
          scene_props: {
            atom_radius: 1.8,
            same_size_atoms: false,
            show_bonds: "always",
          },
          color_scheme: "Jmol",
        }}
        display_mode="structure"
        style="height:100%"
      />
    {/key}
  {:else if activeViewer === "3dmol"}
    <div bind:this={container3dmol} class="w-full h-full rounded-[1.5rem] overflow-hidden">
      {#if !scriptLoaded3dmol}
        <div class="flex items-center justify-center h-full text-slate-500 font-medium bg-slate-900">
          Loading 3Dmol library...
        </div>
      {/if}
    </div>
  {:else if activeViewer === "ngl"}
    <div bind:this={containerNgl} class="w-full h-full rounded-[1.5rem] overflow-hidden">
      {#if !scriptLoadedNgl}
        <div class="flex items-center justify-center h-full text-slate-500 font-medium bg-slate-900">
          Loading NGL Viewer...
        </div>
      {/if}
    </div>
  {:else if activeViewer === "molstar"}
    <div bind:this={containerMolstar} class="w-full h-full rounded-[1.5rem] overflow-hidden relative">
      {#if !scriptLoadedMolstar}
        <div class="flex items-center justify-center h-full text-slate-500 font-medium bg-slate-900">
          Loading Mol* library...
        </div>
      {/if}
    </div>
  {:else if structure}
    {#key structure}
      <Structure
        {structure}
        color_scheme="Jmol"
        performance_mode={isMassive}
        scene_props={{
          atom_radius: 1.8,
          same_size_atoms: false,
          show_bonds: "always",
        }}
        style="height:100%"
      />
    {/key}
  {/if}

  {#if props.sizeMetrics}
    <div class="absolute top-4 left-4 bg-white/80 backdrop-blur-md border border-slate-200 p-3 rounded-lg shadow-soft pointer-events-none z-10">
      <h4 class="text-[10px] font-mono font-bold text-slate-500 uppercase tracking-widest mb-1.5">
        Effective Size (Core)
      </h4>
      <div class="grid grid-cols-2 gap-x-4 gap-y-1 text-sm">
        
        <span class="text-slate-600">Radius</span>
        <span class="font-mono text-brand-600 font-medium text-right">
          {Number((props.sizeMetrics.R_eff_hull / 10).toFixed(2))} nm
        </span>
        
        <span class="text-slate-600">Diameter</span>
        <span class="font-mono text-brand-600 font-medium text-right">
          {Number((props.sizeMetrics.diameter_hull / 10).toFixed(2))} nm
        </span>

      </div>
    </div>
  {/if}

</div>
