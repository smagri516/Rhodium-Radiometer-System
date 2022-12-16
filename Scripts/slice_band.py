import numpy as np

def defineWorkingRegion(y, ymin, ymax):
    """
    Slice wavelength array for working region for a detector
    :param y: whole wavelength array to slice
    :param ymin: minimum wavelength (in m)
    :param ymax: max wavelength
    :return: 2d array of wavelengths in working region with indexes from original
    """

    working_indicies = np.where((y[:] >= ymin) & (y[:] <= ymax))
    y_working = y[working_indicies]
    return np.column_stack((working_indicies[0], y_working))