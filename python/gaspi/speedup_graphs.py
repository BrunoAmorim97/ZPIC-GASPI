import numpy as np
import matplotlib.pyplot as plt


times = [9688, 4871, 2461, 1226, 630, 322]
num_procs = [1, 2, 4, 8, 16, 32]

num_procs_str = list(map(str, num_procs))

speedups = []
efficiency = []

for i in range(len(times)):

	speedup = times[0]/times[i]
	speedups.append( speedup )
	efficiency.append( speedup/num_procs[i] * 100)


print("Speedups:\n" + str(speedups))
print("Efficiencies:\n" + str(efficiency))

plt.title("Weibel [512,512]")



#speedup
# plt.subplot(1, 3, 2)
plt.plot(num_procs_str, speedups, "o-", label="GASPI")

# zip joins x and y coordinates in pairs
for x,y in zip(num_procs_str, speedups):

	label = "(%.1f)" % (y)

	plt.annotate(label, # this is the text
				 (x,y), # this is the point to label
				 textcoords="offset points", # how to position the text
				 xytext=(10,-15), # distance from text to points (x,y)
				 ha='center') # horizontal alignment can be left, right or center



plt.plot(num_procs_str, num_procs, "--", label="1:1")
plt.ylabel('Speedup')
plt.yscale("log")
plt.yticks([])
plt.xlabel('Num Procs')
plt.legend()

plt.show()


