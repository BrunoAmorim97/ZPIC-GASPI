import numpy as np
import matplotlib.pyplot as plt


times = [602.2, 281.0, 124.3, 55.3, 28.1, 14.65, 5.3]
num_procs = [1, 2, 4, 8, 16, 32, 64]

num_procs_str = list(map(str, num_procs))

speedups = []
efficiency = []

for i in range(len(times)):

	speedup = times[0]/times[i]
	speedups.append( speedup )
	efficiency.append( speedup/num_procs[i] * 100)


print(num_procs)
print(times)
print(speedups)
print(efficiency)


plt.title("Weibel [128,128]")

#times
# plt.subplot(1, 3, 1)
# plt.plot(num_procs, times, 'o-')
# plt.ylabel('Execution Times (s)')
# plt.xlabel('Num Procs')
# plt.yscale('log')


#speedup
# plt.subplot(1, 3, 2)
plt.plot(num_procs_str, speedups, "o-")

# zip joins x and y coordinates in pairs
for x,y in zip(num_procs_str, speedups):

	label = "(%.1f)" % (y)

	plt.annotate(label, # this is the text
				 (x,y), # this is the point to label
				 textcoords="offset points", # how to position the text
				 xytext=(0,10), # distance from text to points (x,y)
				 ha='center') # horizontal alignment can be left, right or center

plt.ylabel('Speedup')
plt.xlabel('Num Procs')

#efficiency
# plt.subplot(1, 3, 3)
# plt.plot(num_procs, efficiency, 'o-')
# plt.ylabel('Efficiency (%)')
# plt.xlabel('Num Procs')


plt.show()


