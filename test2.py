import pygimli as pg
import pygimli.meshtools as mt
import pygimli.physics.ert as ert
import numpy as np


# 1. 2D Διατομή με λιγότερα elements
x_bounds = [-15.0, 15.0]
y_bounds = [-15.0, 0.0]

world_2d = mt.createWorld(start=[x_bounds[0], y_bounds[0]], 
                          end=[x_bounds[1], y_bounds[1]], 
                          worldMarker=True)

# Προσθήκη Ανωμαλίας
anomaly_2d = mt.createPolygon([[-1.5, -4.0], [1.5, -4.0], 
                               [1.5, -7.0], [-1.5, -7.0]], isClosed=True)

plc_2d = world_2d + anomaly_2d
plc_2d.addRegionMarker([0.0, -10.0], marker=1, area=5.0)  # Αραίο πλέγμα για χαμηλή RAM
plc_2d.addRegionMarker([0.0, -5.5], marker=2, area=1.0)

mesh_2d = mt.createMesh(plc_2d, quality=31.0)

# 2. Extrusion με ΜΟΝΟ 10 Slices
z_slices = np.linspace(0, 40.0, 11)
mesh = mt.extrudeMesh(mesh_2d, a=z_slices)

# Ορισμός ανωμαλίας στο Z: [15m - 25m]
for cell in mesh.cells():
    if cell.marker() == 2:
        z_center = cell.center()[2]
        if z_center < 15.0 or z_center > 25.0:
            cell.setMarker(1)

print(f"--> Low-Memory Mesh Created! Cells: {mesh.cellCount()} (Πολύ ελαφρύ)")

# ==========================================
# 2. ΗΛΕΚΤΡΟΔΙΑ & ΠΡΩΤΟΚΟΛΛΟ
# ==========================================
floor_x = np.linspace(-1.0, 1.0, 3)    # 3 γραμμές
floor_z = np.linspace(10.0, 30.0, 11)  # 11 θέσεις στο Z

sensors_3d = []
for z in floor_z:
    for x in floor_x:
        sensors_3d.append([x, 0.0, z])

scheme = ert.createData(elecs=sensors_3d, schemeName='dd') # Dipole-Dipole

# ==========================================
# 3. FORWARD MODELING & INVERSION
# ==========================================

rhomap = {1: 1000.0, 2: 10.0}
res = np.array([rhomap.get(cell.marker(), 1000.0) for cell in mesh.cells()])

# Forward Manager
fop = ert.ERTModelling()
fop.setMesh(mesh)
fop.setData(scheme)

voltages = fop.response(res)
res_homog = np.full(mesh.cellCount(), 1000.0)
voltages_homog = fop.response(res_homog)

rhoa = (voltages / voltages_homog) * 1000.0

data = pg.DataContainerERT(scheme)
data["rhoa"] = rhoa
data["err"] = ert.estimateError(data, relativeError=0.03, absoluteUError=1e-4)

data.markInvalid(data("rhoa") <= 0)
data.removeInvalid()

# 3D Inversion
mgr = ert.ERTManager()
inv_mesh = pg.Mesh(mesh)
for cell in inv_mesh.cells():
    cell.setMarker(1)

inv_res = mgr.invert(
    data=data,
    mesh=inv_mesh,
    lam=20,
    maxIter=5,
    verbose=True
)

print("Αποθήκευση VTK")
inv_mesh.exportVTK("tunnel_inverted_light.vtk")
print(" Τέλειωσε!")