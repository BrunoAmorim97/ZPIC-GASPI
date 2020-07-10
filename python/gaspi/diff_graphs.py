import zdf
import sys
import numpy as np
import matplotlib.pyplot as plt

NUM_DIMS = 3
NUM_DATA = 2

dims_name = ["x", "y", "z"]

folder = "/home/bruno/zpic-out/"

emf_folder_gaspi = "EMF-gaspi/"
current_folder_gaspi = "CURRENT-gaspi/"

emf_folder = "EMF/"
current_folder = "CURRENT/"

iteration = sys.argv[1].zfill(6)

emf_files = ["B1-" + iteration + ".zdf", "B2-" + iteration + ".zdf", "B3-" + iteration + ".zdf"]
current_files = ["J1-" + iteration + ".zdf", "J2-" + iteration + ".zdf", "J3-" + iteration + ".zdf"]

data_names = [["Current x", "Current y", "Current z"], ["EMF B x","EMF B y", "EMF B z"]]

NUM_FILES = len(emf_files) + len(current_files)

print(f"Reading files: ", end="")

for file_name in current_files + emf_files:
	print(file_name + " ", end="")
print()

file_paths = []
file_paths_gaspi = []

for file_name in current_files:
	file_paths.append(folder + current_folder + file_name)
	file_paths_gaspi.append(folder + current_folder_gaspi + file_name)

for file_name in emf_files:
	file_paths.append(folder + emf_folder + file_name)
	file_paths_gaspi.append(folder + emf_folder_gaspi + file_name)

data = []
data_gaspi = []

for file_path in file_paths:
	data.append(zdf.read(file_path)[0])

for file_path in file_paths_gaspi:
	data_gaspi.append(zdf.read(file_path)[0])


for i in range(NUM_FILES):
	if (data[i].shape != data_gaspi[i].shape):
		print("Data files have diferent shapes, aborting")
		exit()

num_cells = data[0].size

diffs = []

for data_i in range(NUM_FILES):
	diffs.append(np.absolute(np.array(data[data_i]) - np.array(data_gaspi[data_i])))

nx = 128
extent_size = (0, nx, 0, nx)

fig, axs = plt.subplots(nrows=NUM_DATA, ncols=NUM_DIMS)
fig.suptitle("Absolute error")

for row in range(NUM_DATA):
	for col in range(NUM_DIMS):
		diff = diffs[col + row*NUM_DIMS]
		ax = axs[row, col]
		ax.title.set_text( "%s (max= %.6f)" % (data_names[row][col], np.max(diff)) )
		pcm = ax.pcolormesh(diff)
		fig.colorbar(pcm, ax=axs[row, col])

plt.show()