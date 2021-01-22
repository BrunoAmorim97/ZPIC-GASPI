import numpy as np
import matplotlib.pyplot as plt


# weibel-512x512-1024-500
seq_time = 30040
# num_procs = [1, 14, 28, 42, 56, 70, 84, 98, 112, 126, 140, 154, 168, 182, 196, 210, 224]
# times = [30040, 2205.5, 1113.5, 746, 579, 456, 405.5, 351, 309, 272, 253, 240, 216, 208, 197, 191, 177.5]
# times_row = [30040, 2229.5, 1153, 797.5, 615, 505.5, 419.5, 378, 322, 297, 260, 254.5, 240, 198, 197, 195, 191]
# times_mpi = [30040, 2221.5, 1138, 771.5, 600, 474, 421.5, 359, 314, 285, 272.5, 250, 236.5, 216, 206, 199, 184]
# times_mpi_row = [30040, 2277.5, 1185, 816, 631, 519, 431, 388, 331, 306, 267, 261, 247, 201, 201, 201, 197]

# times_omp_1_14 = [30040, 2603.5, 1349, 927, 684, 550, 459, 401.5, 347, 303.5, 276, 259.5, 230.5, 221, 205.5, 186.5, 174.5]
# times_omp_2_7 = [30040, 2612.5, 1215.5, 812, 614.5, 494, 410, 365, 308, 278, 248, 230, 210, 197, 185, 167, 158]

# lwfa-2000x512-64-4000
seq_time = 19498
num_procs = [1, 28, 56, 84, 112, 140, 168, 196, 224]
times = [19498, 1053, 641, 438, 370.5, 330.5, 302.5, 273.5, 241]
times_mpi = [19498, 1195, 744.5, 616, 435.5, 371.5, 316.5, 276.5, 251]
times_row = [19498, 758, 406, 287, 211.5, 172, 168, 132, 130]
times_mpi_row = [19498, 782, 419.5, 295.5, 216, 174, 170.5, 131, 131]


# seq_time = 19498
# times_row = [453]
# num_procs_row = [47]



num_procs_str = list(map(str, num_procs))

speedups = []
speedups_mpi = []
speedups_row = []
speedups_mpi_row = []

speedups_omp_1_14 = []
speedups_omp_2_7 = []


for i in range(len(times)):

	speedups.append( seq_time/times[i] )
	speedups_mpi.append( seq_time/times_mpi[i] )
	speedups_row.append( seq_time/times_row[i] )
	speedups_mpi_row.append( seq_time/times_mpi_row[i] )

	# speedups_omp_1_14.append( seq_time/times_omp_1_14[i] )
	# speedups_omp_2_7.append( seq_time/times_omp_2_7[i] )

# speedups_row.append(seq_time/times_row[0])

# print("Speedups:\n" + str(speedups))
# print("Speedups MPI:\n" + str(speedups_2))

# plt.title("Weibel Instability (14 procs per node)")
# plt.title("Weibel Instability (14 threads per node)")

plt.title("Laser Wakefield Acceleration (14 procs per node)")



#speedup
plt.plot(num_procs, num_procs, "--", label="1:1")
plt.plot(num_procs, speedups, "o-", label="GASPI")
plt.plot(num_procs, speedups_mpi, ".--", label="MPI")

plt.plot(num_procs, speedups_row, "o--", label="GASPI row-div")
plt.plot(num_procs, speedups_mpi_row, ".--", label="MPI row-div")

# plt.plot(num_procs, speedups_omp_1_14, ".--", label="GASPI 1 proc/node 14 thrd/proc")
# plt.plot(num_procs, speedups_omp_2_7, ".--", label="GASPI 2 proc/node 7 thrd/proc")

# for x,y in zip(num_procs, speedups):
# 	label = "(%.1f)" % (y)
# 	plt.annotate(label, # this is the text
# 				 (x,y), # this is the point to label
# 				 textcoords="offset points", # how to position the text
# 				 xytext=(10,-15), # distance from text to points (x,y)
# 				 ha='center') # horizontal alignment can be left, right or center

# for x,y in zip(num_procs_row, speedups_row):
# 	label = "(%.1f)" % (y)
# 	plt.annotate(label, # this is the text
# 				 (x,y), # this is the point to label
# 				 textcoords="offset points", # how to position the text
# 				 xytext=(20,-3), # distance from text to points (x,y)
# 				 ha='center') # horizontal alignment can be left, right or center

plt.ylabel('Speedup')
# plt.yscale("log")
# plt.xscale("log")
plt.yticks(num_procs)
# num_procs.insert(2, 47)
plt.xticks(num_procs)
plt.xlabel('Num Processes')
# plt.xlabel('Num Threads')
plt.legend()
plt.grid(True)

plt.show()


