import zdf
import numpy as np
import matplotlib.pyplot as plt

# current_file = "current-weibel-small.zdf"

# emf_x_file = "emf-x-weibel-small.zdf"
# emf_y_file = "emf-y-weibel-small.zdf"

current_file = "/home/bruno/zpic-out/CURRENT-gaspi/J3-000500.zdf"

emf_x_file = "/home/bruno/zpic-out/EMF-gaspi/B1-000500.zdf"
emf_y_file = "/home/bruno/zpic-out/EMF-gaspi/B2-000500.zdf"

box_size = 51.2
extent_size = ( 0, box_size, 0, box_size )

(current_data, info) = zdf.read(current_file)
(emf_x_data , info) = zdf.read(emf_x_file)
(emf_y_data , info) = zdf.read(emf_y_file)

plt.figure(figsize=(11, 4))

#Current
plt.subplot(121)
plt.imshow( current_data, interpolation = 'bilinear', origin = 'lower',
			extent = extent_size, aspect = 'auto', cmap = 'Spectral')

plt.colorbar().set_label('Electric Current')
plt.xlabel("$x_1$")
plt.ylabel("$x_2$")
plt.title("Electric Current\n")

# EMF
plt.subplot(122)

Bperp = np.sqrt( emf_x_data**2 + emf_y_data**2 )

plt.imshow( Bperp, interpolation = 'bilinear', origin = 'lower',
		  extent = extent_size, aspect = 'auto', cmap = 'jet')

plt.colorbar().set_label('Magnetic Field')
plt.xlabel("$x_1$")
plt.ylabel("$x_2$")
plt.title("Magnetic Field\n")

plt.show()
# print(data)