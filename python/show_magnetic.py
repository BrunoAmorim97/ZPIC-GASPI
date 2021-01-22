import zdf
import numpy as np
import matplotlib.pyplot as plt
import sys

iteration = str( eval(sys.argv[1]) ).zfill(6)
print(f"Reading files for iteration {iteration}...")


folder = "/home/bruno/zpic-out-saved/zpic-out-weibel-512x512-36-500/serial/EMF/"
# folder = "/home/bruno/zpic-out-saved/lwfa-2000x256-8-1450/serial/EMF/"

extent_size = ( 0, 51.2, 0, 51.2 )
# extent_size = ( 0, 40.0, 0, 51.2 )

range = [[0, 512], [0, 512]]
# range = [[0, 2000], [0, 256]]

files = [folder + "B1-" + iteration + ".zdf", folder + "B2-" + iteration + ".zdf", folder + "B3-" + iteration + ".zdf"]
data = [None] * len(files)

(data[0], info) = zdf.read(files[0])
(data[1], info) = zdf.read(files[1])
(data[2], info) = zdf.read(files[2])


Bperp = np.sqrt( np.array(data[0])**2 + np.array(data[1])**2 )

plt.imshow( Bperp, interpolation = 'bilinear', origin = 'lower',
		extent = ( range[0][0], range[0][1], range[1][0], range[1][1] ),
		aspect = 'auto', cmap = 'jet')


plt.colorbar().set_label('Magnetic Field')
plt.xlabel("$x_1$")
plt.ylabel("$x_2$")
plt.title("Magnetic Field")

plt.show()