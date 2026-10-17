<script>
  import { onMount } from 'svelte';
  import Header from './Header.svelte';
  import Home from './Home.svelte';
  import About from './About.svelte';
  import Contact from './Contact.svelte';
  import Builder from './Builder.svelte';
  import Library from './Library.svelte';
  import DotPage from './DotPage.svelte';

  let currentRoute = $state('home');

  // Library -> Builder hand-off (open a library structure in post-treatment)
  let builderHandoff = $state(null);
  function openInBuilder(handoff) {
    builderHandoff = handoff;
    currentRoute = 'builder';
  }

  // A dot's own page has a link: #/qd/<record id>/<tab>
  let dotId = $state(null);
  let dotTab = $state('thermo');
  function readHash() {
    const m = location.hash.match(/^#\/qd\/([^/]+)(?:\/([\w-]+))?/);
    if (m) {
      dotId = decodeURIComponent(m[1]);
      dotTab = m[2] || 'thermo';
      currentRoute = 'dot';
    } else if (currentRoute === 'dot') {
      currentRoute = 'library';
    }
  }
  const openDot = (id, tab = 'thermo') => { location.hash = `#/qd/${encodeURIComponent(id)}/${tab}`; };
  onMount(() => {
    readHash();
    window.addEventListener('hashchange', readHash);
    return () => window.removeEventListener('hashchange', readHash);
  });
  // leaving the dot page from the header drops its link
  $effect(() => {
    if (currentRoute !== 'dot' && location.hash.startsWith('#/qd/')) history.replaceState(null, '', location.pathname);
  });
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
      <Builder handoff={builderHandoff} onHandoffConsumed={() => (builderHandoff = null)} />
    {:else if currentRoute === 'library'}
      <Library onOpenInBuilder={openInBuilder} onOpenDot={openDot} />
    {:else if currentRoute === 'dot' && dotId}
      <DotPage id={dotId} tab={dotTab} onTab={(t) => openDot(dotId, t)}
               onBack={() => { history.replaceState(null, '', location.pathname); currentRoute = 'library'; }} />
    {/if}
  </div>
</div>
