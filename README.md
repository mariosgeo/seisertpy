# seisertpy: 3D ERT Forward Simulation & Inversion in Horseshoe Tunnels

`seisertpy` provides tools and workflows for **3D Electrical Resistivity Tomography (ERT)** modeling, mesh generation, and inversion in underground tunnel environments using [PyGIMLi](https://www.pygimli.org/).

This project focuses on simulating and resolving subsurface anomalies (e.g., fault zones, water-bearing cavities, or conductive intrusions) beneath the flat floor of a **horseshoe-shaped tunnel** using a multi-line floor electrode array.

---

## 🚀 Key Features

- **Horseshoe Tunnel Geometry (PLC)**:
  - 250° circular roof arch ($r = 3.0\text{ m}$) with a flat bottom floor.
  - Refinement shell surrounding the excavation zone ($r = 4.5\text{ m}$) for accurate potential field computation.
  - Interior treated as an unmeshed void (`addHoleMarker`) to prevent high numerical contrast instabilities and reduce memory footprint.
- **Robust Electrode Node Integration**:
  - Longitudinal multi-line floor array (5 lines $\times$ 21 electrodes = 105 electrodes).
  - 1 mm safety offset into the host rock domain to prevent singularity and boundary node collisions with the interior void.
- **2D-to-3D Extrusion Meshing**:
  - High-quality 2D unstructured triangular cross-section extruded along $Z$ into 3D prismatic elements with strict slice alignment at electrode positions.
  - Anomaly localization in 3D space ($X \in [-1.5, 1.5]\text{ m}$, $Y \in [-2.72, -5.72]\text{ m}$, $Z \in [13.5, 16.5]\text{ m}$).
- **Forward Simulation & Inversion**:
  - Line-by-line Gradient (`'gr'`) array measurement sequence.
  - Forward numerical simulation with Gaussian noise addition ($3\%$).
  - VTK export for 3D rendering and inspection in **ParaView** / **PyVista**.

---

## 📂 Repository Structure

```text
seisertpy/
├── make_mesh_v2.ipynb                # Main Jupyter Notebook (mesh, simulation & inversion)
├── final_tunnel_test.py              # Standalone Python script for full tunnel ERT simulation & inversion
├── tunnel_test_6_greek.py            # Mesh construction and geometry definition script with detailed annotations
├── test2.py                          # Lightweight 3D mesh and inversion test script (low-RAM)
├── horseshoe_tunnel_floor_only.vtk   # Exported 3D tunnel prism mesh (VTK format)
├── simulated_ert_floor_only.dat      # Synthetic apparent resistivity dataset (PyGIMLi format)
├── .gitignore                        # Standard Git exclusions for Python / Jupyter
└── README.md                         # Documentation and usage guide
```

---

## ⚙️ Installation & Requirements

The project relies on `pygimli` and standard scientific Python libraries. We recommend using a dedicated Conda environment:

```bash
# Create conda environment with PyGIMLi from conda-forge
conda create -n pg2 -c conda-forge pygimli numpy matplotlib jupyter

# Activate the environment
conda activate pg2
```

---

## 🔍 Detailed Code Workflow (`make_mesh_v2.ipynb`)

### 1. Geometric Parameterization & Cross-Section
- **Tunnel Cross-Section**: Defined with a $250^\circ$ circular arch ($r = 3.0\text{ m}$) and a flat floor at $Y \approx -1.721\text{ m}$.
- **Surrounding Domain**: Modeled as a $30\text{ m} \times 30\text{ m}$ rock mass.
- **Target Anomaly**: A $3\text{ m} \times 3\text{ m}$ conductive body ($10\ \Omega\cdot\text{m}$) embedded $1\text{ m}$ below the tunnel floor within a background rock domain ($1000\ \Omega\cdot\text{m}$).

### 2. Electrode Alignment & Boundary Safe Offset
- Electrodes are arranged in **5 parallel lines** on the tunnel floor across $X \in [-1.0, 1.0]\text{ m}$ and along $Z \in [5.0, 25.0]\text{ m}$.
- A small offset ($\Delta Y = -1\text{ mm}$) shifts the electrodes into the solid domain (`safe_floor_y = floor_y - 0.001`), ensuring that electrode nodes belong strictly to the solid mesh boundary and do not lie inside the void.

### 3. 2D Meshing & 3D Extrusion
- `pygimli.meshtools.createMesh` creates the 2D triangular cross-section.
- `pygimli.meshtools.extrudeMesh` extrudes the 2D mesh into 3D prisms using coordinate slices aligned with the electrode $Z$-positions and anomaly boundaries.
- The 3D mesh is exported to `horseshoe_tunnel_floor_only.vtk`.

### 4. ERT Forward Simulation
- A Gradient array sequence (`schemeName='gr'`) is generated for each floor line, producing 420 apparent resistivity data points.
- `ert.simulate(mesh, scheme, res=rhomap, noiseLevel=0.03)` calculates synthetic potential responses and adds $3\%$ relative noise.
- Validated synthetic data are stored in `simulated_ert_floor_only.dat`.

### 5. 3D ERT Inversion & Visualization
- The synthetic dataset is inverted using `pygimli.physics.ert.ERTManager` with smoothness regularization (`lam=20`).
- Results and reconstructed resistivity distributions are exported to VTK format for 3D rendering in **ParaView**.

---

## 📊 Visualization in ParaView

To visualize the generated mesh and inversion results:
1. Open [ParaView](https://www.paraview.org/).
2. Load `horseshoe_tunnel_floor_only.vtk`.
3. Apply a **Threshold** filter (by `Marker`) or **Clip** / **Slice** filter along the $Y$ or $Z$ plane to inspect the tunnel void and subsurface anomaly.
