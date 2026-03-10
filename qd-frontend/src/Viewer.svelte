<script>
  import { Structure } from "matterviz";
  import { Trajectory } from "matterviz/trajectory";
  import { parse_xyz } from "matterviz/structure/parse";

  let { xyz = "", isMD = false, dataUrl = "" } = $props();

  // Parse the structure synchronously (used for non-MD single-frame mode)
  let structure = $derived(xyz && !isMD ? parse_xyz(xyz) : null);

  // Only disable bonds if the structure is truly massive (>5000 atoms) to protect the browser
  let isMassive = $derived(structure ? structure.num_atoms > 5000 : false);
</script>

{#if isMD && dataUrl}
  <!-- MD Trajectory mode: use MatterViz Trajectory with speed optimizations -->
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
  <div
    class="p-4 flex items-center justify-center h-full text-slate-500 font-medium"
  >
    No structure loaded
  </div>
{/if}
