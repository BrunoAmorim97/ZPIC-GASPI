import numpy as np
import matplotlib.pyplot as plt


seq_time = 30040
times = [30040, 2205.5, 1113.5, 746, 579, 456, 405.5, 351, 309, 272, 253, 240, 216, 208, 197, 191, 177.5]
times_2 = [30040, 2221.5, 1138, 771.5, 600, 474, 421.5, 359, 314, 285, 272.5, 250, 236.5, 216, 206, 199, 184]
num_procs = [1, 14, 28, 42, 56, 70, 84, 98, 112, 126, 140, 154, 168, 182, 196, 210, 224]

num_procs_str = list(map(str, num_procs))

speedups = []
speedups_2 = []
efficiency = []

for i in range(len(times)):

	speedups.append( seq_time/times[i] )
	speedups_2.append( seq_time/times_2[i] )



print("Speedups:\n" + str(speedups))
print("Speedups MPI:\n" + str(speedups_2))

plt.title("Weibel (14 procs per node)")
# plt.title("Laser Wakefield Acceleration")



#speedup
plt.plot(num_procs_str, num_procs, "--", label="1:1")
plt.plot(num_procs_str, speedups, "o-", label="GASPI")
plt.plot(num_procs_str, speedups_2, ".--", label="MPI")

# zip joins x and y coordinates in pairs
""" for x,y in zip(num_procs_str, speedups):

	label = "(%.1f)" % (y)

	plt.annotate(label, # this is the text
				 (x,y), # this is the point to label
				 textcoords="offset points", # how to position the text
				 xytext=(7,-15), # distance from text to points (x,y)
				 ha='center') # horizontal alignment can be left, right or center """



plt.ylabel('Speedup')
# plt.yscale("log")
# plt.xscale("log")
plt.yticks(num_procs)
plt.xlabel('Num Processes')
plt.legend()
plt.grid(True)

plt.show()


