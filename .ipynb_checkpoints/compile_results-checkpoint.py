import numpy as np
import json, os, sys, pickle
from pathlib import Path

C6 = 2 * np.pi * 862690

top_dir = sys.argv[1]
data_dir = top_dir + "/phase_diagram_data"

f_names = os.listdir(data_dir)

rbs_over_a = []
deltas_over_omega = []
energies = []
ent_entropies = []
densities = []

for f_name in f_names:

  with open(Path(data_dir, f_name), 'r') as f:
    data_dict = json.load(f)

  omega = data_dict["params"]["ham_config"]["omega"]
    
  Rb = (C6 / omega) ** (1 / 6)

  rbs_over_a.append(Rb / data_dict["params"]["ham_config"]["a"])
  deltas_over_omega.append(data_dict["params"]["ham_config"]["delta"] / omega)
  
  energies.append(data_dict["results"]["energy"])
  ent_entropies.append(data_dict["results"]["entanglement_entropy"])
  densities.append(data_dict["results"]["rydberg_density"])

rb_over_a_pts = len(set(rbs_over_a))
delta_over_omega_pts = len(set(deltas_over_omega))

rbs_over_a = np.array(rbs_over_a)
deltas_over_omega = np.array(deltas_over_omega)
energies = np.array(energies)
ent_entropies = np.array(ent_entropies)
densities = np.array(densities)

index_array = np.argsort(rbs_over_a)

rbs_over_a = np.take_along_axis(rbs_over_a, index_array, axis=0)
deltas_over_omega = np.take_along_axis(deltas_over_omega, index_array, axis=0)
energies= np.take_along_axis(energies, index_array, axis=0)
ent_entropies = np.take_along_axis(ent_entropies, index_array, axis=0)

for ii in range(densities.shape[-1]):
    densities[:, ii] = np.take_along_axis(densities[:, ii], index_array, axis=0)

rbs_over_a = np.reshape(rbs_over_a, (rb_over_a_pts, delta_over_omega_pts))
deltas_over_omega = np.reshape(deltas_over_omega, (rb_over_a_pts, delta_over_omega_pts))
energies = np.reshape(energies, (rb_over_a_pts, delta_over_omega_pts))
ent_entropies = np.reshape(ent_entropies, (rb_over_a_pts, delta_over_omega_pts))
densities = np.reshape(densities, (rb_over_a_pts, delta_over_omega_pts, -1))

for ii in range(rb_over_a_pts):
    index_array = np.argsort(deltas_over_omega[ii])

    deltas_over_omega[ii] = np.take_along_axis(deltas_over_omega[ii], index_array, axis=0)
    energies[ii] = np.take_along_axis(energies[ii], index_array, axis=0)

    for jj in range(densities.shape[-1]):
        densities[ii, :, jj] = np.take_along_axis(densities[ii, :, jj], index_array, axis=0)

delta_over_omega_ax = deltas_over_omega[0]
rb_over_a_ax = rbs_over_a[:, 0]

compiled_results = {"rb_over_a_ax": rb_over_a_ax, "delta_over_omega_ax": delta_over_omega_ax, "energies": energies, "entanglement_entropies": ent_entropies, "densities": densities}

out_f_name = top_dir + "/compiled_results.pkl"

with open(out_f_name, 'wb') as io:
   pickle.dump(compiled_results, io)