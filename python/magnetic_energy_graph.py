import zdf
import numpy as np
import matplotlib.pyplot as plt

folder = "/home/bruno/zpic-out/"

emf_folder_gaspi = "gaspi/EMF/"
emf_folder = "serial/EMF/"

# nx = (512, 512)
nx = (2000, 256)

norm = 0.5 * nx[0] * nx[1]


diag_freq = 50
num_iter = 1450

t_max = 35.0

Bperp = np.zeros(int(num_iter/diag_freq))
Bperp_gaspi = np.zeros(int(num_iter/diag_freq))

Bperp_index = 0

for iteration in range(0, num_iter, diag_freq):
	
	iteration_str = str(iteration).zfill(6)
	files = ("B1-" + iteration_str + ".zdf", "B2-" + iteration_str + ".zdf")

	file_paths = (folder + emf_folder + files[0], folder + emf_folder + files[1])
	file_paths_gaspi = (folder + emf_folder_gaspi + files[0], folder + emf_folder_gaspi + files[1])
	
	data = (zdf.read(file_paths[0])[0] , zdf.read(file_paths[1])[0])
	data_gaspi = (zdf.read(file_paths_gaspi[0])[0] , zdf.read(file_paths_gaspi[1])[0])

	Bperp[Bperp_index] = np.sum(data[0]**2 + data[1]**2) * norm
	Bperp_gaspi[Bperp_index] = np.sum(data_gaspi[0]**2 + data_gaspi[1]**2) * norm

	Bperp_index += 1


plt.plot(np.linspace(0, t_max, num = int(num_iter/diag_freq)), Bperp, marker="x")
plt.plot(np.linspace(0, t_max, num = int(num_iter/diag_freq)), Bperp_gaspi, linewidth=2, linestyle='dashed')

plt.yscale('log')
plt.grid(True)
plt.legend(['Serial', 'GASPI'], loc='upper right')
plt.xlabel("$t$ [$1/\omega_n$]")
plt.ylabel("$B_{\perp}$ energy [$m_e c^2$]")
plt.title("Magnetic field energy")
plt.show()