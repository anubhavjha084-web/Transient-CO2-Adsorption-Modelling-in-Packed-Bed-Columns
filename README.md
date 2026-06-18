# Transient CO₂ Adsorption Modelling in Packed Bed Columns

A MATLAB-based simulation tool for modelling transient CO₂ adsorption dynamics in packed bed columns using finite difference methods and ODE solvers.

---

## 📋 Overview

This project simulates the **transient behaviour of CO₂ adsorption** in a packed bed column filled with a solid adsorbent. The model captures gas-phase transport (advection + axial dispersion) coupled with solid-phase adsorption/desorption kinetics.

It includes:
- A **command-line MATLAB script** for quick simulations
- A **GUI-based MATLAB App** for interactive parameter exploration
- A **standalone executable** for users without a MATLAB license

---

## 🔬 Mathematical Model

The simulation solves a system of coupled partial differential equations (PDEs) converted to ODEs via the **method of lines** (finite differences in space, ODE integration in time):

### Gas-Phase Mass Balance

$$\frac{\partial C}{\partial t} = \frac{1}{\varepsilon} \left( D_{ax} \frac{\partial^2 C}{\partial z^2} - u \frac{\partial C}{\partial z} \right) - \frac{(1-\varepsilon)}{\varepsilon} \rho_s \frac{\partial q}{\partial t}$$

### Adsorption Kinetics (Langmuir-type)

$$\frac{\partial q}{\partial t} = k_a \cdot C \cdot (q_{max} - q) - k_d \cdot q$$

### Numerical Methods
- **Spatial discretisation**: Upwind finite differences (1st-order advection) + central differences (2nd-order dispersion)
- **Time integration**: MATLAB `ode15s` (stiff solver)
- **Boundary conditions**: Dirichlet inlet (`C = C_in`), open outlet (`∂C/∂z = 0`)

---

## 📁 Project Structure

| File | Description |
|------|-------------|
| `Term.m` | Main MATLAB script — runs the simulation with default parameters and generates plots |
| `packedBedODE.m` | ODE function defining the packed bed mass balance equations |
| `AdsorptionModelApp.m` | MATLAB App (GUI) for interactive simulation with adjustable parameters |
| `Adsorption_Model_App.exe` | Standalone compiled executable (no MATLAB required) |
| `MyAppInstaller_web.exe` | Web installer for the standalone application |
| `TERM PROJECT.pdf` | Detailed project report and documentation |

---

## 🚀 Getting Started

### Prerequisites

- **MATLAB R2020a** or later (for `.m` scripts)
- No MATLAB required for the standalone `.exe`

### Running the Script

```matlab
% Open MATLAB and navigate to the project directory
cd path/to/TERM_PROJECT

% Run the main simulation
run('Term.m')
```

### Running the GUI App

```matlab
% Launch the interactive app
AdsorptionModelApp
```

### Running the Standalone App

Simply double-click `Adsorption_Model_App.exe` or run `MyAppInstaller_web.exe` to install it.

---

## ⚙️ Default Parameters

| Parameter | Symbol | Default Value | Unit |
|-----------|--------|---------------|------|
| Column length | `L` | 1.0 | m |
| Spatial points | `n` | 200 | — |
| Superficial velocity | `u` | 0.01 | m/s |
| Axial dispersion | `D_ax` | 1×10⁻⁶ | m²/s |
| Bed porosity | `ε` | 0.35 | — |
| Adsorbent density | `ρ_s` | 50 | kg/m³ |
| Max adsorption capacity | `q_max` | 1.0 | mol/kg |
| Adsorption rate constant | `k_a` | 0.01 | 1/s |
| Desorption rate constant | `k_d` | 0.005 | 1/s |
| Inlet concentration | `C_in` | 1.0 | mol/m³ |
| Simulation time | `t_span` | 5000 | s |

---

## 📊 Output Plots

The simulation generates four diagnostic plots:

1. **CO₂ Concentration Profile** — Spatial distribution of gas-phase CO₂ at final time
2. **Breakthrough Curve** — Outlet concentration vs. time (semi-log scale)
3. **Adsorbed Phase Profile** — Spatial distribution of adsorbed CO₂ at final time
4. **Outlet Adsorption vs. Time** — Adsorbed amount at column outlet over time

---

## 🛠️ Customisation

You can modify the simulation parameters either:
- **In the script**: Edit the variables in `Term.m` directly
- **In the GUI**: Use the interactive fields in `AdsorptionModelApp` to change parameters and re-run

---

## 📄 License

This project is part of an academic term project. Please refer to the project report (`TERM PROJECT.pdf`) for full details and references.

---

## 👤 Author

**Anubhav Jha**
- GitHub: [@anubhavjha084-web](https://github.com/anubhavjha084-web)
