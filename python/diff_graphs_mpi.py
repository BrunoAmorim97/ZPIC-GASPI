import zdf
import sys
import numpy as np
import matplotlib.pyplot as plt

NUM_DIMS = 3
NUM_DATA = 2



# directory = "/home/bruno/zpic-out-saved/weibel-512x512-36-500/"
# extent_size = (0, 51.2, 0, 51.2)

directory = "/home/bruno/zpic-out-saved/lwfa-2000x256-8-1450/"
extent_size = (0, 40.0, 0, 51.2)

emf_folder = "EMF/"

folders = ["serial/", "gaspi/", "mpi/"]
NUM_FILES = len(folders)


iteration = str( eval(sys.argv[1]) ).zfill(6)

emf_file = "B2-" + iteration + ".zdf"

data_names = ["Magnetic field GASPI error", "Magnetic field MPI error"]

file_paths = []

for folder in folders:
	file_paths.append(directory + folder + emf_folder + emf_file)

data = []

for file_path in file_paths:
	data.append(zdf.read(file_path)[0])

for i in range(NUM_FILES):
	if (data[0].shape != data[i].shape):
		print("Data files have diferent shapes, aborting")
		exit()



# GASPI diff
plt.subplot(121)
plt.imshow( np.absolute(np.array(data[0]) - np.array(data[1])), interpolation = 'bilinear', origin = 'lower',
			extent = extent_size, aspect = 'auto', cmap = 'jet')

plt.colorbar().set_label('Absolute difference')
plt.title(data_names[0])

# MPI diff
plt.subplot(122)
plt.imshow( np.absolute(np.array(data[0]) - np.array(data[2])), interpolation = 'bilinear', origin = 'lower',
			extent = extent_size, aspect = 'auto', cmap = 'jet')

plt.colorbar().set_label('Absolute difference')
plt.title(data_names[1])



plt.show()