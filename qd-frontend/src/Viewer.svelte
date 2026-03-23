<script>
  import { Structure } from "matterviz";
  import { Trajectory } from "matterviz/trajectory";
  import { parse_xyz } from "matterviz/structure/parse";

  let { xyz = "", isMD = false, dataUrl = "", sizeMetrics = null } = $props();

  // Parse the structure synchronously (used for non-MD single-frame mode)
  let structure = $derived(xyz && !isMD ? parse_xyz(xyz) : null);

  // Only disable bonds if the structure is truly massive (>5000 atoms) to protect the browser
  let isMassive = $derived(structure ? structure.num_atoms > 5000 : false);
</script>

<div class="relative w-full h-full">

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


