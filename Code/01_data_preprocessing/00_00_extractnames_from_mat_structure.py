# R> reticulate::py_require("h5py")
# resolve_trial_order.py
import h5py, csv, os

mice = ["M3411", "M3412",  "M3421F", "M3422F", "M3424F", "M3425F"]   # fill in your mouse folder names
  ############################## MAKE SURE TO EXTRACT AND PLACE THIS FOLDER FROM THE DATA REPO HERE:
base = "Fig_1_3_tri_cir_sqr/"

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

with open("Data_extraction/shapes/trial_order_lookup.csv", "w", newline="") as out:
    w = csv.writer(out)
    w.writerow(["mouse", "trial", "behav_group", "neuron_group"])
    w.writerows(rows)

print(f"wrote trial_order_lookup.csv with {len(rows)} rows")
if missing:
    print(f"\n{len(missing)} file(s) not found, skipped:")
    for m in missing:
        print(" -", m)
