__all__ = [
    "stredni_napeti", "amplituda_napeti", "rozkmit_napeti", "soucinitel_asymetrie",
    "amplituda_z_R", "stredni_z_R", "frekvence", "Kt", "pomerny_gradient",
    "bezrozmerny_gradient",
]


def stredni_napeti(sigma_h, sigma_d):
    """(2.1)"""
    return (sigma_h + sigma_d) / 2


def amplituda_napeti(sigma_h, sigma_d):
    """(2.2)"""
    return (sigma_h - sigma_d) / 2


def rozkmit_napeti(sigma_h, sigma_d):
    """(2.3)"""
    return sigma_h - sigma_d


def soucinitel_asymetrie(sigma_h, sigma_d):
    """(2.4)"""
    return sigma_d / sigma_h


def amplituda_z_R(sigma_h, R):
    """(2.5)"""
    return sigma_h / 2 * (1 - R)


def stredni_z_R(sigma_h, R):
    """(2.5)"""
    return sigma_h / 2 * (1 + R)


def frekvence(T):
    """(2.6)"""
    return 1 / T


def Kt(sigma_max, sigma_n):
    """(2.7)"""
    return sigma_max / sigma_n


def pomerny_gradient(dsigma_dx, sigma_max):
    """(2.8)"""
    return abs(dsigma_dx) / sigma_max


def bezrozmerny_gradient(chi, r):
    """(2.9)"""
    return chi * r
