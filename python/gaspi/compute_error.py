import zdf
import sys
import numpy as np

NUM_DIMS = 3

out_folder = "/home/bruno/zpic-out/"

serial_folder_name = "serial/"
gaspi_folder_name = "gaspi/"


emf_folder = "EMF/"
current_folder = "CURRENT/"

iteration = sys.argv[1].zfill(6)

current_files = ["J1-" + iteration + ".zdf", "J2-" + iteration + ".zdf", "J3-" + iteration + ".zdf"]
emf_b_files = ["B1-" + iteration + ".zdf", "B2-" + iteration + ".zdf", "B3-" + iteration + ".zdf"]
emf_e_files = ["E1-" + iteration + ".zdf", "E2-" + iteration + ".zdf", "E3-" + iteration + ".zdf"]

data_names = ["Current x", "Current y", "Current z", "EMF B x","EMF B y", "EMF B z", "EMF E x","EMF E y", "EMF E z"]

NUM_FILES = len(emf_b_files) + len(emf_e_files) + len(current_files)
NUM_FILES_TYPE = int(NUM_FILES/3)

print(f"Reading files: ", end="")

for file_name in current_files + emf_b_files + emf_e_files:
	print(file_name + " ", end="")
print()

data = []
data_gaspi = []

for file_name in current_files:
	data.append( zdf.read(out_folder + serial_folder_name + current_folder + file_name)[0] )
	data_gaspi.append( zdf.read(out_folder + gaspi_folder_name + current_folder + file_name)[0] )

for file_name in emf_b_files + emf_e_files:
	data.append( zdf.read(out_folder + serial_folder_name + emf_folder + file_name)[0] )
	data_gaspi.append( zdf.read(out_folder + gaspi_folder_name + emf_folder + file_name)[0] )

for i in range(NUM_FILES):
	if (data[i].shape != data_gaspi[i].shape):
		print("ERROR: Data files have diferent shapes, aborting")
		exit()

# print("Serial\n")
# print("Curr")
# for line in data[2]:
# 	for row in line:
# 		print("%.6f " % (row), end = "")
# 	print()

# print("EMF B")
# for line in data[5]:
# 	for row in line:
# 		print("%.6f " % (row), end = "")
# 	print()

# print("\n\nGASPI\n")
# print("Curr")
# for line in data_gaspi[2]:
# 	for row in line:
# 		print("%.6f " % (row), end = "")
# 	print()

# print("EMF B")
# for line in data_gaspi[5]:
# 	for row in line:
# 		print("%.6f " % (row), end = "")
# 	print()

num_cells = data[0].size

error_sum = []
max_error = []

for data_i in range(NUM_FILES):
	abs_diff_matrix = np.absolute(np.array(data[data_i]) - np.array(data_gaspi[data_i]))
	error_sum.append(np.sum(abs_diff_matrix))
	max_error.append(np.max(abs_diff_matrix))



print(f"In {num_cells} cells [{data[0].shape[0]},{data[0].shape[1]}] we got:")
for data_i in range(NUM_FILES):
	if data_i != 0 and data_i % 3 == 0:
		print()
	print("%s error\t Sum: %.10f, Average: %.10f, Max: %.10f" % (data_names[data_i], error_sum[data_i], error_sum[data_i]/num_cells, max_error[data_i]) )

print()

print("Current Total Error Sum = %.10f, Total Average: %.10f, Max %.10f"% (sum(error_sum[0:NUM_FILES_TYPE]), sum([error_sum[i]/num_cells for i in range(NUM_FILES_TYPE)]) / NUM_FILES_TYPE, max(max_error[0:NUM_FILES_TYPE])) )

print("EMF B   Total Error Sum = %.10f, Total Average: %.10f, Max %.10f"% (sum(error_sum[NUM_FILES_TYPE:NUM_FILES_TYPE*2]), sum([error_sum[i]/num_cells for i in range(NUM_FILES_TYPE, NUM_FILES_TYPE*2)]) / NUM_FILES_TYPE, max(max_error[NUM_FILES_TYPE:NUM_FILES_TYPE*2])) )

print("EMF E   Total Error Sum = %.10f, Total Average: %.10f, Max %.10f"% (sum(error_sum[NUM_FILES_TYPE*2:NUM_FILES]), sum([error_sum[i]/num_cells for i in range(NUM_FILES_TYPE*2, NUM_FILES)]) / NUM_FILES_TYPE, max(max_error[NUM_FILES_TYPE*2:NUM_FILES])) )


