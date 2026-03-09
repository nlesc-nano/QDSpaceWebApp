# Quantum Dot Space

An integrated computational ecosystem for building, passivating, and analyzing nanocrystals with atomistic precision. [cite_start]Developed by the **Infante Group** at **BCMaterials**. [cite: 1028, 1050]

## 🚀 Key Features

* [cite_start]**Quantum Dot Library**: A high-fidelity database containing DFT-relaxed semiconductor structures across ABX3, II-VI, III-V, and IV-VI systems. [cite: 573, 577, 579, 784]
* [cite_start]**Interactive Properties**: Real-time visualization of PDOS (Projected Density of States), COOP, and Excited State analyses (Spin-Free and SOC). [cite: 704, 758, 865]
* [cite_start]**Quantum Dot Builder**: A specialized interface for creating core/shell architectures from .cif files with customizable radius and aspect ratios. [cite: 908, 933, 968]
* [cite_start]**miniCAT Passivation Engine**: Automated surface engineering logic that attaches anionic and cationic ligands to ensure global charge neutrality. [cite: 706, 767, 800, 982]
* [cite_start]**Advanced 3D Visualization**: Powered by matterviz, supporting both static structures and optimized Molecular Dynamics (MD) trajectories. [cite: 705, 735, 1077, 1078]

## 🛠️ Project Layout

[cite_start]The frontend is a **Svelte 5** application utilizing **Runes** ($state, $derived) for high-performance reactivity. [cite: 700, 705, 911, 1015]

### Core Components
* [cite_start]**src/Library.svelte**: Database explorer with multi-parameter filtering for material, size, and functional type. [cite: 698, 720, 783]
* [cite_start]**src/Builder.svelte**: The nanocrystal assembly interface and passivation controller. [cite: 906, 912, 1014]
* [cite_start]**src/Viewer.svelte**: A specialized 3D canvas that toggles between single-frame parsing and trajectory streaming. [cite: 1078, 1079, 1082]
* [cite_start]**src/Home.svelte**: Dashboard providing entry points to the ecosystem. [cite: 1034, 1037, 1038]
* [cite_start]**src/Header.svelte**: Navigation and routing state management. [cite: 876, 877, 1013]

### Scientific Logic
* [cite_start]**Stoichiometry Parser**: Dynamically calculates atomic counts and elemental ratios directly from XYZ data strings. [cite: 710, 711, 914]
* [cite_start]**Charge Balancer**: Computes total system charge based on common oxidation states (e.g., Cd2+, Pb2+, Se2-) to identify non-neutral structures. [cite: 698, 713, 906, 918]
* [cite_start]**MD Optimization**: Trajectories are downsampled to a target of ~150 frames to maintain smooth browser performance. [cite: 741, 742, 1082]

## 💻 Development Setup

### Recommended Environment
* [cite_start]**IDE**: VS Code with the Svelte Official Extension. [cite: 559, 566]
* [cite_start]**Compiler Settings**: checkJs is enabled in jsconfig.json to provide advanced typechecking. [cite: 550, 567]

### Scripts
Run the following commands in your terminal:

[cite_start]`npm run dev` - Starts the Vite development server (Proxied to :8000 for API). [cite: 557, 558]

[cite_start]`npm run build` - Compiles the project for production. [cite: 572]

[cite_start]`npm run preview` - Locally preview the production build. [cite: 572]

## 📋 Technical Considerations

* [cite_start]**API Integration**: The frontend proxies /api requests to a FastAPI backend for ligand attachment tasks. [cite: 558, 771]
* [cite_start]**Memory Management**: For massive MD files, the viewer utilizes Blob URLs and performance-mode rendering to prevent browser crashes. [cite: 733, 747, 778, 1081]
* [cite_start]**Styling**: Built with Tailwind CSS using a custom brand palette (Deep Indigo and Coral/Rose). [cite: 553, 554, 555]

---
[cite_start]*Supported by Ikerbasque, the Basque Government, and the European Union.* [cite: 1074, 1075]*

## 📄 License

This project is licensed under the **MIT License**. See [`LICENSE`](LICENSE) for details.
