import scipy.io
import numpy as np

# Load initial MATLAB file
y_vs_e = np.genfromtxt('WavelengthVsEmissivity.csv', delimiter=',', skip_header=1, missing_values=.2)

y = y_vs_e[:,0] * 1e-6 # All wavelengths from optical to lwir spectrum in
e = y_vs_e[:, 1] # All emissivities corresponding to y
y_vs_e = np.column_stack((y, e))

