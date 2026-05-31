<script>
  import { Structure } from "matterviz";
  import { Trajectory } from "matterviz/trajectory";
  import { parse_xyz } from "matterviz/structure/parse";
  import { onMount } from "svelte";

  let { xyz = "", isMD = false, dataUrl = "", sizeMetrics = null } = $props();

  // Parse the structure synchronously (used for non-MD single-frame mode)
  let structure = $derived(xyz && !isMD ? parse_xyz(xyz) : null);

  // Only disable bonds if the structure is truly massive (>5000 atoms) to protect the browser
  let isMassive = $derived(structure ? structure.num_atoms > 5000 : false);

  let activeViewer = $state(isMD ? "matterviz" : "3dmol");
  let container3dmol = $state(null);
  let viewer3dmol = null;
  let scriptLoaded = $state(false);

  // Render 3Dmol structure
  function render3dmol() {
    if (!scriptLoaded || !container3dmol || !xyz || isMD) return;

    // Clean up container
    container3dmol.innerHTML = "";

    if (window.$3Dmol) {
      viewer3dmol = window.$3Dmol.createViewer(container3dmol, {
        backgroundColor: "#0f172a" // Slate-900 matching terminal
      });

      // Add Model
      viewer3dmol.addModel(xyz, "xyz");

      // Styles based on system size
      const atomCount = xyz.split("\n")[0] ? parseInt(xyz.split("\n")[0]) : 0;
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

  // React to changes in loading status, xyz data, and activeViewer choice
  $effect(() => {
    if (activeViewer === "3dmol" && scriptLoaded && xyz && container3dmol) {
      render3dmol();
    }
  });

  onMount(() => {
    if (window.$3Dmol) {
      scriptLoaded = true;
      return;
    }

    const script = document.createElement("script");
    script.src = "https://cdnjs.cloudflare.com/ajax/libs/3dmol/2.2.0/3Dmol-min.js";
    script.async = true;
    script.onload = () => {
      scriptLoaded = true;
    };
    document.head.appendChild(script);

    return () => {
      if (viewer3dmol) {
        viewer3dmol.clear();
      }
    };
  });
</script>

<div class="relative w-full h-full">

  {#if !isMD && xyz}
  <!-- Interactive Selector Toolbar -->
  <div class="absolute top-4 right-4 z-20 flex gap-1 bg-slate-900/80 backdrop-blur-md p-1 rounded-lg border border-slate-700/50">
    <button class="px-2.5 py-1 rounded text-[10px] font-bold transition-all {activeViewer === '3dmol' ? 'bg-brand-600 text-white shadow-sm' : 'text-slate-400 hover:text-slate-200'}"
            onclick={() => activeViewer = '3dmol'}>
      3Dmol (Fast)
    </button>
    <button class="px-2.5 py-1 rounded text-[10px] font-bold transition-all {activeViewer === 'matterviz' ? 'bg-brand-600 text-white shadow-sm' : 'text-slate-400 hover:text-slate-200'}"
            onclick={() => activeViewer = 'matterviz'}>
      MatterViz
    </button>
  </div>
  {/if}

  {#if isMD && dataUrl}
    {#key dataUrl}
      <Trajectory
        data_url={dataUrl}
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
      {#if !scriptLoaded}
        <div class="flex items-center justify-center h-full text-slate-500 font-medium bg-slate-900">
          Loading 3Dmol library...
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
  {:else}
    <div class="p-4 flex items-center justify-center h-full text-slate-500 font-medium">
      No structure loaded
    </div>
  {/if}

  {#if sizeMetrics}
    <div class="absolute top-4 left-4 bg-white/80 backdrop-blur-md border border-slate-200 p-3 rounded-lg shadow-soft pointer-events-none z-10">
      <h4 class="text-[10px] font-mono font-bold text-slate-500 uppercase tracking-widest mb-1.5">
        Effective Size (Core)
      </h4>
      <div class="grid grid-cols-2 gap-x-4 gap-y-1 text-sm">
        
        <span class="text-slate-600">Radius</span>
        <span class="font-mono text-brand-600 font-medium text-right">
          {Number((sizeMetrics.R_eff_hull / 10).toFixed(2))} nm
        </span>
        
        <span class="text-slate-600">Diameter</span>
        <span class="font-mono text-brand-600 font-medium text-right">
          {Number((sizeMetrics.diameter_hull / 10).toFixed(2))} nm
        </span>

      </div>
    </div>
  {/if}

</div> ```


