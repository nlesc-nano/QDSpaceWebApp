<script>
  import Header from './Header.svelte';
  import Home from './Home.svelte';
  import About from './About.svelte';
  import Contact from './Contact.svelte';
  import Builder from './Builder.svelte';
  import Library from './Library.svelte';

  let currentRoute = 'home'; 
</script>

<div class="min-h-screen bg-slate-50 text-slate-900 flex flex-col">
  <Header bind:currentRoute={currentRoute} />

  <div class="flex-grow">
    {#if currentRoute === 'home'}
      <Home bind:currentRoute={currentRoute} />
    {:else if currentRoute === 'about'}
      <About />
    {:else if currentRoute === 'contact'}
      <Contact />
      
    {:else if currentRoute === 'builder'}
      <Builder />
    {:else if currentRoute === 'library'}
      <Library />
    {/if}
    
  </div>
</div>

