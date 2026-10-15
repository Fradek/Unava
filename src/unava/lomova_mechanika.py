import math

__all__ = [
    "sigma_z_rovinna_deformace", "kappa", "K_I", "K_II", "K_III", "pole_mod_I",
    "pole_mod_II", "pole_mod_III", "hlavni_napeti_rovinna", "plasticka_zona_I",
    "plasticka_zona_II", "plasticka_zona_III", "rp_I", "rp_II", "rp_III",
    "mala_plasticka_zona", "platnost_KIC", "kriterium_K", "KIC_barson_rolfe",
    "KIC_sailors_corten", "energie_bez_trhliny", "E_efektivni", "griffith_WT",
    "plocha_trhliny", "griffith_kriticke_napeti", "griffith_kriticka_delka",
    "hnaci_sila", "hnaci_sila_poddajnost", "G_I", "G_II", "G_III", "sih_Sc",
    "sih_Smin_II", "ctod", "ctod_elasticko_plasticky", "J_z_K", "J_z_ctod",
    "ramberg_osgood_hrr",
]


def sigma_z_rovinna_deformace(sigma_x, sigma_y, mu):
    """(6.1)"""
    return mu * (sigma_x + sigma_y)


def kappa(mu, rovinna_deformace=True):
    """(6.6)"""
    return 3 - 4 * mu if rovinna_deformace else (3 - mu) / (1 + mu)


def K_I(sigma, a, Y=1.0):
    """(6.7)"""
    return sigma * math.sqrt(math.pi * a) * Y


def K_II(tau, a, Y=1.0):
    """(6.12)"""
    return tau * math.sqrt(math.pi * a) * Y


def K_III(tau, a, Y=1.0):
    """(6.15)"""
    return tau * math.sqrt(math.pi * a) * Y


def pole_mod_I(KI, r, theta, G=None, mu=None, rovinna_deformace=True):
    """(6.5)  vrací (σx, σy, τxy[, u, v])"""
    f = KI / math.sqrt(2 * math.pi * r)
    c, s = math.cos(theta / 2), math.sin(theta / 2)
    s3 = math.sin(3 * theta / 2)
    vysledek = (f * c * (1 - s * s3), f * c * (1 + s * s3), f * c * s * math.cos(3 * theta / 2))
    if G is None:
        return vysledek
    k = kappa(mu, rovinna_deformace)
    g = KI / (2 * G) * math.sqrt(r / (2 * math.pi))
    return vysledek + (g * c * (k - 1 + 2 * s ** 2), g * s * (k + 1 - 2 * c ** 2))


def pole_mod_II(KII, r, theta, G=None, mu=None, rovinna_deformace=True):
    """(6.11)  vrací (σx, σy, τxy[, u, v])"""
    f = KII / math.sqrt(2 * math.pi * r)
    c, s = math.cos(theta / 2), math.sin(theta / 2)
    c3, s3 = math.cos(3 * theta / 2), math.sin(3 * theta / 2)
    vysledek = (-f * s * (2 + c * c3), f * s * c * c3, f * c * (1 - s * s3))
    if G is None:
        return vysledek
    k = kappa(mu, rovinna_deformace)
    g = KII / (2 * G) * math.sqrt(r / (2 * math.pi))
    return vysledek + (g * s * (k + 1 + 2 * c ** 2), -g * c * (k - 1 - 2 * s ** 2))


def pole_mod_III(KIII, r, theta, G=None):
    """(6.14)  vrací (τxz, τyz[, w])"""
    f = KIII / math.sqrt(2 * math.pi * r)
    vysledek = (-f * math.sin(theta / 2), f * math.cos(theta / 2))
    if G is None:
        return vysledek
    return vysledek + (KIII / G * math.sqrt(2 * r / math.pi) * math.sin(theta / 2),)


def hlavni_napeti_rovinna(sigma_x, sigma_y, tau_xy):
    """(6.16)"""
    s = (sigma_x + sigma_y) / 2
    r = math.sqrt(((sigma_x - sigma_y) / 2) ** 2 + tau_xy ** 2)
    return s + r, s - r


def plasticka_zona_I(KI, sigma_Y, theta, mu=None):
    """(6.20)  mu=None → rovinná napjatost, jinak rovinná deformace"""
    c2, s2 = math.cos(theta / 2) ** 2, math.sin(theta / 2) ** 2
    if mu is None:
        return KI ** 2 / (2 * math.pi * sigma_Y ** 2) * c2 * (1 + 3 * s2)
    return KI ** 2 / (2 * math.pi * sigma_Y ** 2) * c2 * ((1 - 2 * mu) ** 2 + 3 * s2)


def plasticka_zona_II(KII, sigma_Y, theta, mu=None):
    """(6.22)  mu=None → rovinná napjatost, jinak rovinná deformace"""
    f = 3 * KII ** 2 / (4 * math.pi * sigma_Y ** 2)
    if mu is None:
        return f * ((7 - math.cos(theta)) / 3 - 1.5 * math.sin(theta) ** 2)
    return f * (2 + (1 - 2 * mu) ** 2 * (1 - math.cos(theta)) / 3 - 1.5 * math.sin(theta) ** 2)


def plasticka_zona_III(KIII, sigma_Y):
    """(6.26)  τY = σY/√3"""
    tau_Y = sigma_Y / math.sqrt(3)
    return KIII ** 2 / (2 * math.pi * tau_Y ** 2)


def rp_I(KI, sigma_Y):
    """(6.21)"""
    return (KI / sigma_Y) ** 2 / math.pi


def rp_II(KII, sigma_Y):
    """(6.23)"""
    return 3 * (KII / sigma_Y) ** 2 / math.pi


def rp_III(KIII, sigma_Y):
    """(6.27)"""
    return (KIII / (sigma_Y / math.sqrt(3))) ** 2 / math.pi


def mala_plasticka_zona(sigma, sigma_Y):
    """(6.28)"""
    return sigma / sigma_Y < 0.3


def platnost_KIC(B, a, W, KIC, sigma_Y):
    """(6.30)"""
    mez = 2.5 * (KIC / sigma_Y) ** 2
    return B >= mez and a >= mez and (W - a) >= mez


def kriterium_K(KI, KIC, k=2.0):
    """(6.31)"""
    return KI < KIC / k


def KIC_barson_rolfe(E, KV):
    """(6.33)"""
    return 0.0149 * math.sqrt(E) * KV ** 0.75


def KIC_sailors_corten(E, KV):
    """(6.34)"""
    return 0.0321 * math.sqrt(E * KV)


def energie_bez_trhliny(sigma, E, V):
    """(6.35)"""
    return sigma ** 2 / (2 * E) * V


def E_efektivni(E, mu=None):
    """(6.40)  mu=None → rovinná napjatost"""
    return E if mu is None else E / (1 - mu ** 2)


def griffith_WT(sigma, a, B, E_ef):
    """(6.39)"""
    return math.pi * B * sigma ** 2 * a ** 2 / E_ef


def plocha_trhliny(a, B, pruchozi=True):
    """(6.42)"""
    return 4 * a * B if pruchozi else 2 * a * B


def griffith_kriticke_napeti(E_ef, gamma, a):
    """(6.43)"""
    return math.sqrt(2 * E_ef * gamma / (math.pi * a))


def griffith_kriticka_delka(E_ef, gamma, sigma):
    """(6.44)"""
    return 2 * E_ef * gamma / (math.pi * sigma ** 2)


def hnaci_sila(sigma, a, E_ef):
    """(6.45)"""
    return math.pi * sigma ** 2 * a / E_ef


def hnaci_sila_poddajnost(F, B, dC_da, konstantni_posuv=False):
    """(6.48), (6.49)"""
    G = F ** 2 / (2 * B) * dC_da
    return -G if konstantni_posuv else G


def G_I(KI, E, mu=None):
    """(6.50)  mu=None → rovinná napjatost"""
    return KI ** 2 / E if mu is None else (1 - mu ** 2) * KI ** 2 / E


def G_II(KII, E, mu):
    """(6.51)"""
    return (1 - mu ** 2) * KII ** 2 / E


def G_III(KIII, E, mu):
    """(6.51)"""
    return (1 + mu) * KIII ** 2 / E


def sih_Sc(KIC, mu, Gs):
    """(6.53)"""
    return KIC ** 2 * (1 - 2 * mu) / (4 * math.pi * Gs)


def sih_Smin_II(tau, a, mu, Gs):
    """(6.54)"""
    return tau ** 2 * a / (12 * Gs) * (2 * (1 - mu) - mu ** 2)


def ctod(V, W, a, r):
    """(6.58)"""
    return r * (W - a) * V / (r * (W - a) + a)


def ctod_elasticko_plasticky(KI, Ye, E_ef, Vp, W, a, rp):
    """(6.59)"""
    return KI ** 2 / (Ye * E_ef) + ctod(Vp, W, a, rp)


def J_z_K(KI, E_ef):
    """(6.62)"""
    return KI ** 2 / E_ef


def J_z_ctod(sigma_Y, delta):
    """(6.63)"""
    return sigma_Y * delta


def ramberg_osgood_hrr(sigma, sigma_0, E, alpha, n):
    """(6.64)  vrací ε"""
    eps_0 = sigma_0 / E
    return eps_0 * (sigma / sigma_0 + alpha * (sigma / sigma_0) ** n)
