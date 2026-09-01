# STEP0: Download the large original file `neuronalindividuals_new.mat` and
# Place into the folder `Data/Matlab_originals_M3424F`, together with `behav.mat

# STEP1: This script produces a .csv files with labels to ease the extraction of .mat data into R
# First, run IN R the following
# py_require("h5py")
# import("h5py")
# Then, run in Python
import h5py, csv, os

mice = ["M3424F"]   # fill in your mouse folder names
############################## MAKE SURE TO EXTRACT AND PLACE THIS FOLDER FROM THE DATA REPO HERE:
base = "../Data/Matlab_originals"

rows = []
missing = []

for mouse in mice:
    mouse_clean = mouse.strip("/")  # guard against accidental leading/trailing slash
    behav_path  = os.path.join(base, mouse_clean, "behav.mat")
    neuron_path = os.path.join(base, mouse_clean, "neuronIndividuals_new.mat")

    if not os.path.exists(behav_path):
        missing.append(behav_path)
        continue
    if not os.path.exists(neuron_path):
        missing.append(neuron_path)
        continue

    with h5py.File(behav_path, "r") as f:
        behav_order = [f[f["behavIndividuals"][i,0]].name.split("/")[-1] for i in range(6)]

    with h5py.File(neuron_path, "r") as f:
        neuron_order = [f[f["neuronIndividuals_new"][i,0]].name.split("/")[-1] for i in range(6)]

    for i in range(6):
        rows.append([mouse_clean, i+1, behav_order[i], neuron_order[i]])

with open("../Data/Matlab_originals/order_lookup_M3424_F.csv", "w", newline="") as out:
    w = csv.writer(out)
    w.writerow(["mouse", "trial", "behav_group", "neuron_group"])
    w.writerows(rows)

print(f"wrote order_lookup_M3424_F.csv with {len(rows)} rows")
if missing:
    print(f"\n{len(missing)} file(s) not found, skipped:")
    for m in missing:
        print(" -", m)

