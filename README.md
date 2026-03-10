
# Quantum Dot Space

An interactive web platform for **building, exploring, and analyzing quantum dots**.

Quantum Dot Space provides an integrated environment for browsing large databases of atomistic nanocrystal structures, visualizing them interactively in the browser, and constructing new quantum dots with automated surface passivation and stoichiometric control.

The platform is designed to accelerate research in **computational nanoscience and quantum dot discovery** by providing an intuitive interface for exploring large structure spaces and associated electronic properties.

---

## Overview

Quantum Dot Space combines two main components:

1. **Quantum Dot Library**  
   A searchable database of precomputed quantum dots including structural and electronic properties.

2. **Quantum Dot Builder**  
   A tool for constructing new nanocrystals with customizable size, morphology, facets, and ligand passivation.

The application runs entirely in the browser and provides real-time 3D visualization and interactive property analysis.

---

## Features

- Interactive **3D visualization** of quantum dot structures
- Exploration of large **quantum dot databases**
- Real-time filtering by material, size, functionalization, and run type
- Interactive **electronic structure analysis**
- Visualization of:
  - fuzzy band structures
  - projected density of states (PDOS)
  - crystal orbital overlap population (COOP)
  - excited states
- Automated **quantum dot builder** with:
  - size control
  - facet selection
  - shape control (sphere, rod, platelet, custom)
  - ligand passivation
  - charge neutrality enforcement
- Download of generated structures

---

# Interface

## Landing Page

The main interface provides access to the two core components of the platform.

![Home](docs/images/home.png)

---

## Quantum Dot Library

The library allows users to browse a large database of precomputed nanocrystals and inspect their structural and electronic properties.

Users can filter structures by:

- material family
- size
- functionalization
- simulation run type

Selected structures are displayed in a fully interactive 3D viewer.

![Library](docs/images/library.png)

---

## Electronic Structure Analysis

The platform provides interactive visualization of electronic properties including band structure, PDOS, and bonding analysis.

Available analyses include:

- Fuzzy band structures
- Projected density of states (PDOS)
- Crystal orbital overlap population (COOP)

![Electronic Properties](docs/images/electronic_properties.png)

---

## Photophysics Analysis

Excited-state properties can also be explored interactively, including exciton characteristics and energy decomposition.

![Photophysics](docs/images/photophysics.png)

---

## Quantum Dot Builder

The builder allows users to generate new nanocrystal structures directly in the browser.

Users can:

- upload a crystal structure (`.cif`)
- control nanocrystal size
- specify facet energies
- adjust aspect ratio
- automatically passivate dangling bonds
- ensure global charge neutrality

The generated structure is visualized in real time.

![Builder](docs/images/builder.png)

---

# Installation

Clone the repository:

```bash
git clone https://github.com/nlesc-nano/QDSpaceWebApp.git
cd QDSpaceWebApp
```

Install frontend dependencies:

```bash
cd qd-frontend
npm install
```

Run the development server:

```bash
npm run dev
```

The application will start at:

```
http://localhost:5173
```

---

# Project Structure

```
QDSpaceWebApp
│
├── qd-frontend
│   ├── src
│   │   ├── components
│   │   ├── routes
│   │   └── utilities
│   │
│   ├── public
│   │   ├── ABX3
│   │   ├── II-VI
│   │   └── III-V
│   │
│   ├── package.json
│   └── vite.config.js
│
├── scripts
└── README.md
```

The `public` directory contains large collections of quantum dot structures used by the library.

---

# Technology Stack

The web application is built using modern frontend technologies:

- **Svelte**
- **Vite**
- **JavaScript / TypeScript**
- **WebGL-based molecular visualization**

The platform is designed for high-performance rendering of atomistic structures directly in the browser.

---

# Data

The repository includes collections of quantum dot structures grouped by material class:

| Dataset | Description |
|-------|-------------|
| ABX3 | Perovskite quantum dots |
| II-VI | II–VI semiconductor nanocrystals |
| III-V | III–V semiconductor nanocrystals |

Each structure includes atomic coordinates and associated metadata used for filtering and visualization.

---

# Development

Start the development server:

```bash
npm run dev
```

Build production version:

```bash
npm run build
```

Preview production build:

```bash
npm run preview
```

---

# Contributing

Contributions are welcome. Possible areas for improvement include:

- improved visualization tools
- additional materials datasets
- new electronic structure analysis modules
- performance optimizations

---

# License

This project is released under the MIT License.
=======
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
