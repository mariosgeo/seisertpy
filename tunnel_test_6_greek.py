import pygimli as pg
import pygimli.meshtools as mt
import numpy as np

# Ορισμός παραμέτρων
tunnel_radius = 3.0 # ακτίνα τούνελ σε μέτρα
tunnel_length = 50.0 # μήκος τούνελ σε μέτρα 

# Ορισμός παραμέτρων ert
n_floor_lines = 5 # αριθμός γραμμών ert
floor_x_positions = np.linspace(-1.0, 1.0, n_floor_lines) #βάζω μέσα στο τούνελ τις 5 γραμμές μου με απόσταση 0.5 η μία από την άλλη οπότε από -1 ως 1
n_floor_electrodes = 21 # αριθμός ηλεκτροδίων σε κάθε γραμμή
floor_z_positions = np.linspace(15, 35, n_floor_electrodes)# οριοθετούμε τις θέσεις των ηλεκτροδίων στο πάτωμα, έχουμε 21 ηλεκτρόδια με 1μετρο απόσταση μεταξύ τους άρα σύνολο 20 μέτρα
# για να μπουν στο κέντρο του τούνελ υπολογίζουμε οτι η γραμμή πρέπει να ξεκινάει από το 15 μέτρο και να τελειώνει στο 35

# Δημιουργουμε συμμετρικά το πρώτο και το τελευταίο σημείο από το ημικύκλιο σχήμα του τούνελ --- 
#συνολικό ανάπτυγμα τόξου του τούνελ
start_rad = np.radians(90 - (250 / 2))# -35 μοίρες το πρώτο σημείο
end_rad = np.radians(90 + (250 / 2))# 215 μοίρες το τελευταίο σημείο 

# Υπολογισμός σημείων του ημικυκλίου της σήραγγας
angles = np.linspace(start_rad, end_rad) #δημιουργεί γωνίες από -35-215 μοίρες, πρακτικά το άνω μισό της σήραγγας
arch_pts = [[tunnel_radius * np.cos(a), tunnel_radius * np.sin(a)] for a in angles]  # μετατρέπει τις γωνίες σε καρτεσιανές συντεταγμένες από τους τύπους x=Rcosθ y=Rsinθ 

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

#δημιουργώ κόμβους στο πάτωμα 
for x in floor_x_positions:
    plc_2d.createNode([x, floor_y])

# διαχωρίζω τα σημεία ωστε να ορίσω σε αυτά την ποιότητα του πλέγματος που θέλω 
# 1.εξωτερικός χώρος-κύβος (max area 5.0)
plc_2d.addRegionMarker([0, 25], marker=1, area=5.0)

# 2.ενδιάμεσος χώρος-πλέγμα (max area 0.05)
# είναι ο χώρος από την οροφή του τούνελ(3) ως την οροφή του δευτερου ημικυκλίου που φτιάξαμε (4.5)
plc_2d.addRegionMarker([0, 3.75], marker=2, area=0.05)

# 3. ανωμαλία υπεδάφους (max area 0.2)
plc_2d.addRegionMarker([0, anom_top - 1.5], marker=3, area=0.2)

# Ορίζουμε το εσωτερικό του τούνελ ωσ κενό/τρύπα. 
# Έτσι δεν θα γεμίσει με πλέγμα αυτή η περιοχή .
#plc_2d.addHoleMarker([0, 1.0])# δηλώνει ότι το κέντρο  του κυλίνδρου βρίσκεται μέσα σε κενό (άρα όλο)
plc_2d.addHoleMarker([0, 0, tunnel_length / 2])

# Δημιουργία του 2D πλέγματος
print("Generating 2D cross-section with anomaly...")
mesh_2d = mt.createMesh(plc_2d, quality=34.0)#δεν προσθέτουμε το area=1.0 όπως στους προηγούμενους κώδικες οπότε το κάθε marker που ορίσαμε π

# δημιουργώ πιο μικρές θέσεις γύρω από τα ηλεκτρόδια ώστε να έχω μεγαλύτερη ανάλυση γύρω από αυτά 
refined_z_slices = []
for z in floor_z_positions:
    refined_z_slices.extend([z - 0.2, z, z + 0.2])

# οριοθετούμε ακριβώς την ανωμαλία στον Ζ άξονα γιατί αλλιώς θα βγει 3*3*50
# θα είναι 3x3x3m κύβος στο κέντρο του τούνελ (Z = 25)
anomaly_z_start = 23.5
anomaly_z_end = 26.5

#δημιουργώ 25 βασικές θέσεις κατά μήκος του τούνελ 
base_z_slices = np.linspace(0, tunnel_length, 25)
#ενώνω τις βασικές με τις πιο μικρές θέσεις που έφτιαξα
all_z_targets = np.concatenate((base_z_slices, refined_z_slices, [anomaly_z_start, anomaly_z_end]))

# σορτάρω τις θέσεις αυτές και σβήνω τυχόν κοινές θέσεις. βεβαιώνομαι ότι όλες οι θέσεις είνια εντός του τούνελ (από 0 ως 50)
z_slices = np.unique(np.sort(all_z_targets))
z_slices = z_slices[(z_slices >= 0) & (z_slices <= tunnel_length)]

#εξώθηση του πλέγματος
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

# προσθήκη σενσορα σε κάθε θέση ηλεκτροδίου που φτιάξαμε πριν 
sensors_3d = []

for z in floor_z_positions:
    for x in floor_x_positions:
        sensors_3d.append([x, floor_y, z])

# Export σε VTK
#mesh_3d.exportVTK(r"D:\pygimli-diplomatiki\tunnel_mesh\tunnel_test_6_greek.vtk")
mesh_3d.exportVTK(r"/export/home/mnotia/final_inverted_tunnel.vtk")
print("All done! Results saved to final_inverted_tunnel.vtk")