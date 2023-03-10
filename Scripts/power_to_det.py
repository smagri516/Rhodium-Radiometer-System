import math
from scipy.integrate import simps


def system_specs(Rdet, Rspot, Fno, zML_red):
    """
    By changing the system inputs, the output will give Rlens, zLD_red, f1
    :param Rdet:
    :param Rspot:
    :param Fno:
    :param zML_red:
    :return:
    """
    m = Rdet / Rspot
    zLD_red = zML_red * m  # reduced thickness of lens to detector

    f1 = 1 / ((1 / zLD_red) - (1 / zML_red))  # Calculate focal length of lens
    Dlens = f1 / Fno  # Calculate diameter of lens using F/#
    Rlens = Dlens / 2

    return Rlens, f1, zLD_red


def radianceToIncidentPower(Ltotal, rhoa, rhob, aPrime, tWindow, Rdet, Rlens, zLD_red):
    """
    Given the system params and the total radiance, calculate the incident
    power on the detector
    :param Ltotal: Total radiance
    :param rhoa: reflection of first window (0-1)
    :param rhob: reflection of second window (0-1)
    :param aPrime: absorption coeff of window over det
    :param tWindow: thickness of window
    :param Rdet: Radius of detector
    :param Rlens: Radius of Lens
    :param zLD_red: reduced thickness from lens to detector
    :return: float power in Watts incident upon detector
    """

    Tau = ((1 - rhoa) * (1 - rhob) * math.exp(-aPrime * tWindow)) / (1 - rhoa * rhob * math.exp(-2 * aPrime * tWindow))
    Coeff = Tau * (math.pi * Rlens ** 2) * (math.pi * Rdet ** 2) / (zLD_red ** 2)
    Phi = Coeff * Ltotal # Total power incident on the detector

    return Phi


def systemCoeff(rhoa, rhob, aPrime, tWindow, Rdet, Rlens, zLD_red):
    """
    Given the system params and the total radiance, calculate the incident
    power on the detector
    :param rhoa: reflection of first window (0-1)
    :param rhob: reflection of second window (0-1)
    :param aPrime: absorption coeff of window over det
    :param tWindow: thickness of window
    :param Rdet: Radius of detector
    :param Rlens: Radius of Lens
    :param zLD_red: reduced thickness from lens to detector
    :return: float power in Watts incident upon detector
    """

    Tau = ((1 - rhoa) * (1 - rhob) * math.exp(-aPrime * tWindow)) / (1 - rhoa * rhob * math.exp(-2 * aPrime * tWindow))
    Coeff = Tau * 1e-2 *(math.pi * Rlens ** 2) * (math.pi * Rdet ** 2) / (zLD_red ** 2)
    return Coeff

def calculateSignal(y_working, sysCoeff):
    """
    Calculate the signal from numerical integration of
    :param y_working:  take 2d spectral radiance vs R(y) and integrate. Then multiply by system coeff
    :return:
    """
    RxL = y_working[:, 2] * y_working[:, 3] # Ly * Ry
    signal = sysCoeff * simps(RxL, y_working[:,1])
    return signal


def calc_totalNoise(noises):
    """
    Calculate total noises from a list of noise values
    :param noises: [noise1, noise2, ...]
    :return: total noise value
    """

    i_noise = 0
    for noise in noises:
        i_noise += noise**2
    i_noise = math.sqrt(i_noise)
    return i_noise

def calculate_1f(idc, fTop, fBot, k=10):
    """
    Calculate 1/f noise from signal
    :param idc:
    :return:
    """

    return k * idc*1e-1 * math.log(fTop/fBot)

def calc_ijohnson(T, band, R, kb=1.38e-23):
    """Calculate ijohnson"""

    return math.sqrt(4 * kb * T * band / R)


