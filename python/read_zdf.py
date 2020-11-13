import zdf
import numpy as np
import matplotlib.pyplot as plt
import sys

iteration = str( eval(sys.argv[1]) ).zfill(6)
print(f"Reading files for iteration {iteration}...")


current_file_gaspi = "/home/bruno/zpic-out/gaspi/CURRENT/J2-" + iteration + ".zdf"
emf_b_x_file_gaspi = "/home/bruno/zpic-out/gaspi/EMF/B2-" + iteration + ".zdf"

current_file_serial = "/home/bruno/zpic-out/serial/CURRENT/J2-" + iteration + ".zdf"
emf_b_x_file_serial = "/home/bruno/zpic-out/serial/EMF/B2-" + iteration + ".zdf"

extent_size = ( 0, 40.0, 0, 51.2 )

(current_data_gaspi, info) = zdf.read(current_file_gaspi)
(emf_b_x_data_gaspi , info) = zdf.read(emf_b_x_file_gaspi)

(current_data_serial, info) = zdf.read(current_file_serial)
(emf_b_x_data_serial , info) = zdf.read(emf_b_x_file_serial)


# Current serial
plt.subplot(231)
plt.imshow( current_data_serial, interpolation = 'bilinear', origin = 'lower',
			extent = extent_size, aspect = 'auto', cmap = 'jet')

plt.colorbar().set_label('Electric Current')
plt.title("Electric Current serial\n")

# Current GASPI
plt.subplot(232)
plt.imshow( current_data_gaspi, interpolation = 'bilinear', origin = 'lower',
			extent = extent_size, aspect = 'auto', cmap = 'jet')

plt.colorbar().set_label('Electric Current')
plt.title("Electric Current GASPI\n")

# Current diff
plt.subplot(233)
plt.imshow( np.absolute(np.array(current_data_serial) - np.array(current_data_gaspi)), interpolation = 'bilinear', origin = 'lower',
			extent = extent_size, aspect = 'auto', cmap = 'jet')

plt.colorbar().set_label('Electric Current')
plt.title("Electric Current Diff\n")

# EMF serial
plt.subplot(234)

# Bperp = np.sqrt( emf_b_x_data**2 + emf_y_data**2 )

plt.imshow( emf_b_x_data_serial, interpolation = 'bilinear', origin = 'lower',
		  extent = extent_size, aspect = 'auto', cmap = 'jet')

plt.colorbar().set_label('Magnetic Field')
plt.title("Magnetic Field serial\n")

# EMF GASPI
plt.subplot(235)

# Bperp = np.sqrt( emf_b_x_data**2 + emf_y_data**2 )

plt.imshow( emf_b_x_data_gaspi, interpolation = 'bilinear', origin = 'lower',
		  extent = extent_size, aspect = 'auto', cmap = 'jet')

plt.colorbar().set_label('Magnetic Field')
plt.title("Magnetic Field GASPI\n")

# EMF diff
plt.subplot(236)
plt.imshow( np.absolute(np.array(emf_b_x_data_serial) - np.array(emf_b_x_data_gaspi)), interpolation = 'bilinear', origin = 'lower',
			extent = extent_size, aspect = 'auto', cmap = 'jet')

plt.colorbar().set_label('Magnetic Field')
plt.title("Magnetic Field Diff\n")



plt.show()