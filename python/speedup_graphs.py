import numpy as np
import matplotlib.pyplot as plt


seq_time = 796
times = [792, 416, 290, 255, 137, 134]
num_procs = [1, 2, 3, 4, 7, 8]

num_procs_str = list(map(str, num_procs))

speedups = []
efficiency = []

for i in range(len(times)):

	speedup = seq_time/times[i]
	speedups.append( speedup )
	efficiency.append( speedup/num_procs[i] * 100)


print("Speedups:\n" + str(speedups))
print("Efficiencies:\n" + str(efficiency))

# plt.title("Weibel [512,512]")
plt.title("Laser Wakefield Acceleration")



#speedup
# plt.subplot(1, 3, 2)
plt.plot(num_procs, speedups, "o-", label="GASPI")

# zip joins x and y coordinates in pairs
for x,y in zip(num_procs, speedups):

	label = "(%.1f)" % (y)

	plt.annotate(label, # this is the text
				 (x,y), # this is the point to label
				 textcoords="offset points", # how to position the text
				 xytext=(7,-12), # distance from text to points (x,y)
				 ha='center') # horizontal alignment can be left, right or center



plt.plot(num_procs, num_procs, "--", label="1:1")
plt.ylabel('Speedup')
# plt.yscale("log")
# plt.xscale("log")
plt.yticks([])
plt.xlabel('Num Procs')
plt.legend()

plt.show()


