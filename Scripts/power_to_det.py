import math

def system_specs(Rdet, Rspot, Fno, zML_red):
    """
    By changing the system inputs, the output will give Rlens, zLD_red, f1
    :param Rdet:
    :param Rspot:
    :param Fno:
    :param zML_red:
    :return:
    """
    m = Rdet/Rspot
    zLD_red = zML_red * m # reduced thickness of lens to detector

    f1 = 1 / ((1/zLD_red) - (1/zML_red)) # Calculate focal length of lens
    Dlens = f1/Fno # Calculate diameter of lens using F/#
    Rlens = Dlens/2

    return Rlens, f1, zLD_red








