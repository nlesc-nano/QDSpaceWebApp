# Quantum Dot Space

An integrated computational ecosystem for building, passivating, and analyzing nanocrystals with atomistic precision. This project bridges the gap between theoretical DFT modeling and experimental reality through a **FastAPI backend** (miniCAT engine) and a **Svelte 5 frontend**.

## 🛠️ Installation & Setup

This project requires a dual-environment configuration to function correctly.

### 1. Backend (Python Environment)
The backend executes heavy computational logic, including the miniCAT ligand attachment engine and structure generation.

* **Requirements**: Python 3.10+ and the dependencies listed in `requirements.txt`.
* **Setup**:
    ```bash
    pip install -r requirements.txt
    ```
* **Operation**: The frontend is configured via a proxy to expect a FastAPI server running on `http://127.0.0.1:8000` to handle `/api` requests.

### 2. Frontend (Svelte + Vite)
The frontend provides the interactive 3D environment, stoichiometry analysis, and database filtering logic.

* **Setup**:
    ```bash
    npm install
    ```
* **Development**:
    ```bash
    npm run dev
    ```
* **Production Build**:
    ```bash
    npm run build
    ```

---

## 🚀 Key Features

### Quantum Dot Builder
Construct complex core/shell architectures from scratch.
* **Geometry Control**: Define radius, aspect ratios (Sphere, Platelet, Rod), and Miller index facets (e.g., 100, 111).
* **miniCAT Passivation**: Real-time ligand attachment via SMILES strings to achieve global charge neutrality.
* **System Breakdown**: Dynamic reporting of total atoms, shell composition, and ligand counts.

### Quantum Dot Library
A curated repository of pre-computed, DFT-relaxed structures.
* **Intelligent Filtering**: Sort through ABX3, II-VI, III-V, and IV-VI systems by material, size, or functional type.
* **Interactive Properties**: Integrated view of PDOS, COOP, and Excited State (Spin-Free/SOC) plots.
* **Charge Balancer**: Automated detection of non-neutral structures using a built-in oxidation state database.

### 3D Visualization (MatterViz)
* **Static Structures**: High-performance rendering of Jmol-colored atomistic models.
* **MD Trajectories**: Optimized playback of Molecular Dynamics simulations, downsampled to ~150 frames to ensure smooth 60fps browser performance.

---

## 📁 Technical Architecture

* **Framework**: Svelte 5 utilizing **Runes** ($state, $derived, $effect) for high-speed reactivity.
* **Engine**: `matterviz` for WebGL-based 3D rendering of static and trajectory data.
* **Styling**: Tailwind CSS with custom branding including "Deep Indigo" and "Coral/Rose" palettes.
* **Proxying**: Configured via `vite.config.js` to route `/api` calls directly to the FastAPI server.

---
*Developed by the **Infante Group** at **BCMaterials**. Supported by Ikerbasque and the European Union.*
