import pygimli as pg
import pygimli.meshtools as mt
import numpy as np
import pygimli.physics.ert as ert


# ==========================================
# 1. ΠΑΡΑΜΕΤΡΟΙ
# ==========================================

# Ορισμός παραμέτρων
tunnel_radius = 3.0       # ακτίνα τούνελ σε μέτρα
tunnel_length = 50.0      # μήκος τούνελ σε μέτρα


# Ορισμός παραμέτρων ERT
n_floor_lines = 5         # αριθμός γραμμών ERT

# βάζω μέσα στο τούνελ τις 5 γραμμές μου
# με απόσταση 0.5 η μία από την άλλη, από -1 ως 1
floor_x_positions = np.linspace(-1.0, 1.0, n_floor_lines)

n_floor_electrodes = 21   # αριθμός ηλεκτροδίων σε κάθε γραμμή

# οριοθετούμε τις θέσεις των ηλεκτροδίων στο πάτωμα
# έχουμε 21 ηλεκτρόδια με 1 m απόσταση μεταξύ τους
# η γραμμή ξεκινάει από το 15 m και τελειώνει στο 35 m
floor_z_positions = np.linspace(15, 35, n_floor_electrodes)

# ==========================================
# 2. ΓΕΩΜΕΤΡΙΑ ΤΟΥ ΤΟΥΝΕΛ
# ==========================================

# Δημιουργούμε συμμετρικά το πρώτο και το τελευταίο
# σημείο από το ημικύκλιο σχήμα του τούνελ

# συνολικό ανάπτυγμα τόξου του τούνελ
deg = 250

# -35 μοίρες το πρώτο σημείο
start_angle_deg = 90 - (deg / 2)

# 215 μοίρες το τελευταίο σημείο
end_angle_deg = 90 + (deg / 2)

# μετατροπή των σημείων σε ακτίνια
start_rad = np.radians(start_angle_deg)
end_rad = np.radians(end_angle_deg)


# Υπολογισμός σημείων του ημικυκλίου της σήραγγας
angles = np.linspace(start_rad, end_rad)
# μετατρέπει τις γωνίες σε καρτεσιανές συντεταγμένες
# x = R cosθ
# y = R sinθ
arch_pts = [[tunnel_radius * np.cos(a), tunnel_radius * np.sin(a)] for a in angles] 

# Ορίζουμε το ύψος του δαπέδου
floor_y = tunnel_radius * np.sin(start_rad)
# δημιουργούμε ένα εξτρά ημικύκλιο μέσα στο οποίο θα μπει το πλέγμα 
ref_radius = 4.5
ref_arch_pts = [[ref_radius * np.cos(a), ref_radius * np.sin(a)] for a in angles]

# Δημιουργία εξωτερικού χώρου/βράχου γύρω από το τούνελ ---
# Δημιουργία τετραγώνου 60*60 με γωνίες χ, y 
domain_2d = mt.createPolygon([[-30, -30], [30, -30], [30, 30], [-30, 30]], isClosed=True)

ref_shell_2d = mt.createPolygon(ref_arch_pts, isClosed=True)
# Δημιουργία του ημικυκλίου της σήραγγας ενώνοντας τα σημεία που δημιουργήσαμε πριν
# ενώνει το πρώτο με το τελευταίο σημείο με ευθεία γραμμή και έτσι δημιουργεί το πάτωμα της σήραγγας
tunnel_2d = mt.createPolygon(arch_pts, isClosed=True)

# δημιουργώ ανωμαλία στο υπέδαφος
# με κέντρο X=0, ξεκινάει 1m κάτω από το πάτωμα floor_y, ως 3m
anom_top = floor_y - 1.0
anom_bot = floor_y - 4.0
anomaly_2d = mt.createPolygon([[-1.5, anom_top], [1.5, anom_top], 
                               [1.5, anom_bot], [-1.5, anom_bot]], isClosed=True)

# δημιουργία του PLC από την ένωση του κύβου, των δύο κυλίνδρων και της ανωμαλίας
plc_2d = domain_2d + ref_shell_2d + tunnel_2d + anomaly_2d



# --- 1mm Safe Shift ---
# Μετατοπίζουμε τα ηλεκτρόδια 1mm μέσα στο πέτρωμα (προς τα κάτω στον άξονα Y)
# ώστε οι κόμβοι να μην βρίσκονται ακριβώς πάνω στο ανοιχτό όριο του κενού (hole).
shift = 0.001
safe_floor_y = floor_y - shift

# ==========================================
# 3. ΠΛΕΓΜΑ
# ==========================================

# δημιουργώ κόμβους στο πάτωμα (μετατοπισμένους κατά 1mm προς τα μέσα)
for x in floor_x_positions:
    plc_2d.createNode([x, safe_floor_y])

# διαχωρίζω τα σημεία ωστε να ορίσω σε αυτά την ποιότητα του πλέγματος που θέλω 
# 1.εξωτερικός χώρος-κύβος (max area 5.0)
plc_2d.addRegionMarker([0, 25], marker=1, area=5.0)

# 2.ενδιάμεσος χώρος-πλέγμα (max area 0.05)
# είναι ο χώρος από την οροφή του τούνελ(3) ως την οροφή του δευτερου ημικυκλίου που φτιάξαμε (4.5)
plc_2d.addRegionMarker([0, 3.75], marker=2, area=0.05)

# 3. ανωμαλία υπεδάφους (max area 0.2)
plc_2d.addRegionMarker([0, anom_top - 1.5], marker=3, area=0.2)

# Ορίζουμε το εσωτερικό του τούνελ ως κενό/τρύπα (σε 2D συντεταγμένες).
# Έτσι δεν θα γεμίσει με πλέγμα αυτή η περιοχή.
plc_2d.addHoleMarker([0, 1.0])

# Δημιουργία του 2D πλέγματος
print("Generating 2D cross-section with anomaly...")
mesh_2d = mt.createMesh(plc_2d, quality=34.0)

plc_2d.exportVTK("plc_2d.vtk")  
mesh_2d.exportVTK("mesh_2d.vtk")

# δημιουργώ πιο μικρές θέσεις γύρω από τα ηλεκτρόδια ώστε να έχω μεγαλύτερη ανάλυση γύρω από αυτά 
refined_z_slices = []
for z in floor_z_positions:
    refined_z_slices.extend([z - 0.2, z, z + 0.2])

# οριοθετούμε ακριβώς την ανωμαλία στον Ζ άξονα γιατί αλλιώς θα βγει 3*3*50
# θα είναι 3x3x3m κύβος στο κέντρο του τούνελ (Z = 25)
anomaly_z_start = 23.5
anomaly_z_end = 26.5

# δημιουργώ 25 βασικές θέσεις κατά μήκος του τούνελ 
base_z_slices = np.linspace(0, tunnel_length, 25)
# ενώνω τις βασικές με τις πιο μικρές θέσεις που έφτιαξα
all_z_targets = np.concatenate((base_z_slices, refined_z_slices, [anomaly_z_start, anomaly_z_end]))

# σορτάρω τις θέσεις αυτές και σβήνω τυχόν κοινές θέσεις. βεβαιώνομαι ότι όλες οι θέσεις είναι εντός του τούνελ (από 0 ως 50)
z_slices = np.unique(np.sort(all_z_targets))
z_slices = z_slices[(z_slices >= 0) & (z_slices <= tunnel_length)]

# εξώθηση του πλέγματος
print("Extruding 3D mesh... This might take a moment.")
mesh_3d = mt.extrudeMesh(mesh_2d, a=z_slices)

# οριοθέτηση της ανωμαλίας
# προσέχουμε που βρίσκεται το κέντρο των κελιών της ώστε να μην βγαίνει έξω από τα όρια που έχω ορίσει
print("Localizing the anomaly volume...")
for cell in mesh_3d.cells():
    if cell.marker() == 3:
        z_center = cell.center()[2]
        # αν το όριο της είναι έξω από τα όρια να μετατραπεί σε βράχο
        if z_center < anomaly_z_start or z_center > anomaly_z_end:
            cell.setMarker(1)

print(f"3D Mesh successfully generated!")
print(f"Nodes: {mesh_3d.nodeCount()}")
print(f"Cells: {mesh_3d.cellCount()}")

# Export σε VTK
mesh_3d.exportVTK("final_tunnel.vtk")
print("Mesh exported to 'final_tunnel.vtk'")

# ==========================================
# 4. FORWARD MODEL
# ==========================================
print("\n3. Generating Gradient array protocol (line by line)...")

# Δημιουργία των 5 γραμμών ηλεκτροδίων (21 ηλεκτρόδια ανά γραμμή)
all_lines = []
for x in floor_x_positions:
    line = [[x, safe_floor_y, z] for z in floor_z_positions]
    all_lines.append(line)

# Κατασκευή DataContainerERT
scheme = pg.DataContainerERT()
flat_sensors = [pos for line in all_lines for pos in line]
for pos in flat_sensors:
    scheme.createSensor(pos)

# Κατασκευή μετρήσεων Gradient ανά γραμμή
data_idx = 0
offset = 0
for line_sensors in all_lines:
    ls = ert.createData(elecs=line_sensors, schemeName='gr')
    for i in range(ls.size()):
        scheme.createFourPointData(data_idx, 
                                   int(ls('a')[i]) + offset,
                                   int(ls('b')[i]) + offset,
                                   int(ls('m')[i]) + offset,
                                   int(ls('n')[i]) + offset)
        data_idx += 1
    offset += len(line_sensors)

print(f"Total Sensors: {scheme.sensorCount()}")
print(f"Total Measurements: {scheme.size()}")

rhomap = [
    [1, 1000.0],  # Background rock
    [2, 1000.0],  # Refined rock shell
    [3, 10.0],    # Conductive anomaly
    [0, 1000.0]   # Boundary
]

print("Running Forward ERT Simulation...")
data = ert.simulate(
    mesh_3d,
    scheme=scheme,
    res=rhomap,
    noiseLevel=0.03,
    noiseAbs=1e-4,
    seed=42
)

data.markInvalid(data("rhoa") <= 0)
data.removeInvalid()
print(f"Simulation success! Usable data points: {data.size()}")

# ==========================================
# 5. INVERSION
# ==========================================
print("Preparing Inversion Mesh...")
inv_mesh = pg.Mesh(mesh_3d)

for cell in inv_mesh.cells():
    if cell.marker() == 3:
        cell.setMarker(2)

print("Starting 3D Inversion... Monitor your computer's RAM usage!")
mgr = ert.ERTManager()
inv_res = mgr.invert(data, mesh=inv_mesh, lam=20, verbose=True)

inv_mesh.exportVTK("final_inverted_tunnel.vtk")
print("All done!")