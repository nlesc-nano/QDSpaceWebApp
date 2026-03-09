import { defineConfig } from 'vite'
import { svelte } from '@sveltejs/vite-plugin-svelte'

export default defineConfig({
  plugins: [svelte()],
  server: {
    proxy: {
      // Routes any frontend fetch('/api/...') directly to FastAPI
      '/api': 'http://127.0.0.1:8000' 
    }
  }
})

