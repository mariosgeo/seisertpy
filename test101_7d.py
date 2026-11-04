import os
from pathlib import Path
import pygimli as pg
import pygimli.meshtools as mt
import numpy as np
import pygimli.physics.ert as ert
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# ============================================================
# 1. ΟΡΙΣΜΟΣ ΦΑΚΕΛΟΥ & VTK FOLDER
# ============================================================
current_dir = Path(__file__).parent
os.chdir(current_dir)

output_dir = current_dir / "simulation_3d_forward_2d_inversion_outputs"
output_dir.mkdir(parents=True, exist_ok=True)

summary_file_path = output_dir / "inversion_summary_3d_fwd_2d_inv.txt"
with open(summary_file_path, "w", encoding="utf-8") as f:
    f.write("Line_Number\tX_Position(m)\tChi2_Misfit\tRMS_Error\n")
    f.write("-" * 55 + "\n")

# ==========================================
# 2. ΠΑΡΑΜΕΤΡΟΙ (ΓΡΑΜΜΗ ΣΤΟ x = -1.0 m)
# ==========================================
tunnel_radius = 3.0
tunnel_length = 70.0

floor_x_positions = [+0.5]  # Αλλαγή από 0.0 σε -1.0
n_floor_electrodes = 21
floor_z_positions = np.linspace(25, 45, n_floor_electrodes)

# ==========================================
# 3. ΓΕΩΜΕΤΡΙΑ 3D ΤΟΥΝΕΛ & 3D MESH (ΣΥΜΜΕΤΡΙΚΟ)
# ==========================================
deg = 250
start_angle_deg = 90 - (deg / 2)
end_angle_deg = 90 + (deg / 2)
start_rad = np.radians(start_angle_deg)
end_rad = np.radians(end_angle_deg)
angles = np.linspace(start_rad, end_rad)

arch_points = [[tunnel_radius * np.cos(a), tunnel_radius * np.sin(a)] for a in angles]
floor_y = tunnel_radius * np.sin(start_rad)

ref_radius = 4.5
ref_arch_points = [[ref_radius * np.cos(a), ref_radius * np.sin(a)] for a in angles]

world_2d = mt.createWorld(start=[-30, -30], end=[30, 30], isClosed=True)
tunnel_2d = mt.createPolygon(arch_points, isClosed=True)
ref_shell_2d = mt.createPolygon(ref_arch_points, isClosed=True)

anom_top = floor_y - 1.0
anom_bot = floor_y - 4.0
anomaly_2d = mt.createPolygon([[-1.5, anom_top], [1.5, anom_top], 
                              [1.5, anom_bot], [-1.5, anom_bot]], isClosed=True)

plc_2d = world_2d + ref_shell_2d + tunnel_2d + anomaly_2d
plc_2d.addHoleMarker(pg.RVector3(0, 0, 0))

shift = 0.001
safe_floor_y = floor_y - shift

for x in floor_x_positions:
    plc_2d.createNode([x, safe_floor_y])

plc_2d.addRegionMarker([0, -10], marker=1, area=100.0)      # Background rock
plc_2d.addRegionMarker([0, 3.75], marker=2, area=20.0)      # Refined rock shell
plc_2d.addRegionMarker([0, anom_top - 1.5], marker=3, area=0.5) # Anomaly
plc_2d.addHoleMarker([0, 1.0])

print("Δημιουργία Συμμετρικού 3D Mesh...")
mesh_2d = mt.createMesh(plc_2d, quality=32)

z_slices = np.linspace(0, tunnel_length, 35)
for z in floor_z_positions:
    if z not in z_slices:
        z_slices = np.append(z_slices, z)

anomaly_z_min, anomaly_z_max = 33.5, 36.5
z_slices = np.append(z_slices, [anomaly_z_min, anomaly_z_max])
z_slices = np.unique(np.sort(z_slices))

mesh_3d = mt.extrudeMesh(mesh_2d, a=z_slices)

for cell in mesh_3d.cells():
    if cell.marker() == 3:
        z = cell.center().z()
        if z < anomaly_z_min or z > anomaly_z_max:
            cell.setMarker(1)

# --- Αποθήκευση του 3D Mesh σε VTK ---
vtk_mesh_3d_path = output_dir / "mesh_3d.vtk"
mesh_3d.exportVTK(str(vtk_mesh_3d_path))
print(f"Το 3D πλέγμα αποθηκεύτηκε στο: {vtk_mesh_3d_path}")

# ============================================================
# 4. ΚΑΤΑΣΚΕΥΗ SCHEME ΚΑΙ 3D FORWARD ΠΡΟΣΟΜΟΙΩΣΗ
# ============================================================
all_lines_3d = []
for x in floor_x_positions:
    line = [[x, safe_floor_y, z] for z in floor_z_positions]
    all_lines_3d.append(line)

scheme = pg.DataContainerERT()
flat_sensors = [pos for line in all_lines_3d for pos in line]
for pos in flat_sensors:
    scheme.createSensor(pos)

data_idx = 0
offset = 0
line_info = []

for line_sensors in all_lines_3d:
    ls = ert.createData(elecs=line_sensors, schemeName='dd')
    start_idx = data_idx
    for i in range(ls.size()):
        scheme.createFourPointData(data_idx, 
                                   int(ls('a')[i]) + offset,
                                   int(ls('b')[i]) + offset,
                                   int(ls('m')[i]) + offset,
                                   int(ls('n')[i]) + offset)
        data_idx += 1
    end_idx = data_idx
    line_info.append((line_sensors, start_idx, end_idx))
    offset += len(line_sensors)

rhomap = [
    [1, 1000.0],  # Background rock
    [2, 1000.0],  # Refined rock shell
    [3, 10.0],    # Conductive anomaly
    [0, 1000.0]   # Boundary 
]

print("Εκτέλεση 3D Forward ERT Simulation...")
full_data = ert.simulate(
    mesh_3d,
    scheme=scheme,
    res=rhomap,
    noiseLevel=0.03,
    noiseAbs=1e-4,
    seed=42
)

# =========================================================================== 
# 5. ΑΠΟΜΟΝΩΣΗ ΔΕΔΟΜΕΝΩΝ ΚΑΙ 2D ΑΝΤΙΣΤΡΟΦΗ
# ===========================================================================
c_min_val = 100.0
c_max_val = 1000.0

line_sensors_3d, start_idx, end_idx = line_info[0]
x_val = line_sensors_3d[0][0]

print(f"\n==========================================")
print(f"ΕΚΤΕΛΕΣΗ 2D ΑΝΤΙΣΤΡΟΦΗΣ ΓΙΑ ΤΗ ΓΡΑΜΜΗ (x = {x_val:.2f} m)")
print(f"==========================================")

line_sensors_2d = [[s[2], s[1]] for s in line_sensors_3d]
sub_data = ert.createData(elecs=line_sensors_2d, schemeName='dd')

sub_data['rhoa'] = full_data['rhoa'][start_idx:end_idx]

sub_data.markInvalid((sub_data("rhoa") <= 0) | (sub_data("rhoa") > 2000))
sub_data.removeInvalid()
sub_data["err"] = 0.03 

mgr = ert.ERTManager(sub_data)

inv_res = mgr.invert(sub_data, lam=20, paraDepth=6.0, paraMaxCellSize=0.3, verbose=True)

chi2_val = mgr.inv.chi2()            
pred_rhoa = mgr.inv.response     
rms_val = np.sqrt(np.mean(((sub_data['rhoa'] - pred_rhoa) / sub_data['rhoa']) ** 2)) * 100

print(f"--> 2D Inversion | Chi2: {chi2_val:.4f} | RMS Error: {rms_val:.2f}%")

with open(summary_file_path, "a", encoding="utf-8") as f:
    f.write(f"1\t\t{x_val:.2f}\t\t{chi2_val:.4f}\t\t{rms_val:.2f}%\n")

vtk_filename = output_dir / f"inversion_2d_from_3fwd_line_x_{x_val:.2f}.vtk"
mgr.paraDomain.addData("Resistivity", inv_res)
mgr.paraDomain.exportVTK(str(vtk_filename))

ax, cbar = mgr.showResult(
    cMin=c_min_val,
    cMax=c_max_val,
    logScale=True,
    cMap="jet",
    colorBar=True,
    coverage=1
)

if cbar is not None:
    cbar.set_label(r"Resistivity ($\Omega$m)")
else:
    sm = plt.cm.ScalarMappable(cmap="jet", norm=LogNorm(vmin=c_min_val, vmax=c_max_val))
    sm.set_array([])
    fig = ax.get_figure()
    cbar = fig.colorbar(sm, ax=ax, orientation="vertical", pad=0.02)
    cbar.set_label(r"Resistivity ($\Omega$m)")

ax.set_xlim(25.0, 45.0)

current_yticks = ax.get_yticks()
new_labels = [f"{abs(y - safe_floor_y):.1f}" if y <= safe_floor_y else "" for y in current_yticks]
ax.set_yticklabels(new_labels)

ax.set_xlabel("Distance along tunnel (z in m)")
ax.set_ylabel("Depth from floor (m)")
ax.set_title(f"2D ERT Inversion (3D Forward) - Line (x = {x_val:.2f} m - $\chi^2$={chi2_val:.2f}, RMS={rms_val:.1f}%)")

fig = ax.get_figure()
png_filename = output_dir / f"inversion_2d_line_x_{x_val:.2f}.png"
pdf_filename = output_dir / f"inversion_2d_line_x_{x_val:.2f}.pdf"

fig.savefig(output_dir / "inversion_2d_line_auto.png", dpi=300, bbox_inches="tight")
plt.savefig(png_filename, dpi=300, bbox_inches="tight")
plt.savefig(pdf_filename, dpi=300, bbox_inches="tight")
plt.close(fig)
mgr.saveResult(str(output_dir / '2d_inversion_result'))

print(f"\nΗ 2D αντιστροφή ολοκληρώθηκε με επιτυχία!")