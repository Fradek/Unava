import math

__all__ = [
    "basquin_napeti", "basquin_cykly", "basquin_regrese", "modifikovana_mez_unavy",
    "eta_D", "ohybove_napeti", "gradient_ohyb", "vg_heywood", "vg_nemec", "psi",
    "sigma_F", "haigh_dovolena_amplituda", "haigh_ekvivalentni_amplituda",
    "smith_horni_dolni", "gerber_weyrauch_horni_dolni", "Kf", "citlivost_vrubu",
    "Kf_z_citlivosti", "Kf_z_gradientu", "gough_pollard", "gough_pollard_krehke",
    "gough_zobecnene", "liu_zenner", "lee_empiricke", "alpha_np", "sigma_eq_tah_krut",
    "sigma_eq_rankine", "sigma_eq_guest", "alpha_hmh", "hmh_amplituda", "sigma_eq_hmh",
    "hydrostaticke_napeti", "tau_oct", "sigma_eq_sines", "sigma_eq_lee", "dang_van",
    "fdp", "findley", "findley_basquin", "mcdiarmid", "biaxialita",
    "skalovane_normalove", "ETA_L", "ETA_R",
]


def basquin_napeti(sigma_f, b, Nf):
    """(3.1)  σa = σf·(2Nf)^b"""
    return sigma_f * (2 * Nf) ** b


def basquin_cykly(sigma_a, sigma_f, b):
    """(3.1)  Nf = ½·(σa/σf)^(1/b)"""
    if sigma_a <= 0:
        return math.inf
    return 0.5 * (sigma_a / sigma_f) ** (1 / b)


def basquin_regrese(sigma_a, Nf):
    """(3.2)  vrací (σf, b) z lineární regrese log σa = log σf + b·log 2Nf"""
    x = [math.log10(2 * n) for n in Nf]
    y = [math.log10(s) for s in sigma_a]
    xs, ys = sum(x) / len(x), sum(y) / len(y)
    b = sum((xi - xs) * (yi - ys) for xi, yi in zip(x, y)) / sum((xi - xs) ** 2 for xi in x)
    return 10 ** (ys - b * xs), b


def modifikovana_mez_unavy(sigma_c, *soucinitele):
    """(3.3)  σc× = σc·ηD·ηL·vg·ηR"""
    for eta in soucinitele:
        sigma_c *= eta
    return sigma_c


def eta_D(D, d):
    """(3.4)"""
    return 1 - math.sqrt(0.02 * math.log(D / d))


ETA_L = {"tah": 0.9, "tah s ohybem": 0.7, "ohyb": 1.0, "krut ocel": 0.58, "krut litina": 0.8}
"""Tab. 3.1"""

ETA_R = {0.5: 1.0, 0.9: 0.897, 0.95: 0.868, 0.99: 0.814, 0.999: 0.753,
         0.9999: 0.702, 0.99999: 0.659, 0.999999: 0.620}
"""Tab. 3.2"""


def ohybove_napeti(Mo, d, Jy):
    """(3.5)"""
    return Mo * d / (2 * Jy)


def gradient_ohyb(d):
    """(3.6)"""
    return 2 / d


def vg_heywood(c, d):
    """(3.7)"""
    return 1 + math.sqrt(2 * c / d)


def vg_nemec(sigma_max, sigma_s, c=None, d=None):
    """(3.8), s c a d (3.9)"""
    pomer = (sigma_max - sigma_s) / sigma_s
    if (c is None) != (d is None):
        raise ValueError("pro rovnici (3.9) je nutné zadat c i d")
    if c is not None:
        pomer *= 2 * c / d
    return 1 + math.sqrt(pomer)


def psi(sigma_c, sigma_hc):
    """(3.10)"""
    return 2 * sigma_c / sigma_hc - 1


def sigma_F(sigma_c, sigma_hc):
    """(3.10)"""
    return sigma_c / psi(sigma_c, sigma_hc)


def _mocnina(x, exponent):
    return math.copysign(abs(x) ** exponent, x)


def _haigh_clen(sigma_m, metoda, sigma_u=None, sigma_Y=None, n=None, k=None):
    if metoda == "goodman":
        return sigma_m / sigma_u
    if metoda == "soderberg":
        return sigma_m / sigma_Y
    if metoda == "gerber":
        return (sigma_m / sigma_u) ** 2
    if metoda == "asme":
        return (sigma_m / sigma_u) ** 2
    if metoda == "bagci":
        return (sigma_m / sigma_Y) ** 4
    if metoda == "clemson":
        return _mocnina(sigma_m / sigma_Y, n)
    if metoda == "sekercioglu":
        return _mocnina(sigma_m / sigma_u, k)
    raise ValueError(f"neznámá metoda: {metoda!r}")


def haigh_dovolena_amplituda(sigma_m, sigma_c, metoda, sigma_u=None, sigma_Y=None,
                             n=None, k=None, tlak_bez_prinosu=True):
    """(3.11)–(3.17)  mezní σa pro dané σm"""
    if tlak_bez_prinosu and sigma_m < 0:
        sigma_m = 0
    clen = _haigh_clen(sigma_m, metoda, sigma_u, sigma_Y, n, k)
    if clen >= 1:
        return 0.0
    if metoda == "asme":
        return sigma_c * math.sqrt(1 - clen)
    return sigma_c * (1 - clen)


def haigh_ekvivalentni_amplituda(sigma_a, sigma_m, metoda, sigma_u=None, sigma_Y=None,
                                 n=None, k=None, tlak_bez_prinosu=False):
    """(3.11)–(3.17)  ekvivalentní amplituda při σm = 0 (posun S-N křivky, kap. 3.3.2)"""
    if tlak_bez_prinosu and sigma_m < 0:
        sigma_m = 0
    clen = _haigh_clen(sigma_m, metoda, sigma_u, sigma_Y, n, k)
    if clen >= 1:
        return math.inf
    if metoda == "asme":
        return sigma_a / math.sqrt(1 - clen)
    return sigma_a / (1 - clen)


def smith_horni_dolni(sigma_m, sigma_c, sigma_u):
    """(3.18)  vrací (σh, σd) na mezi Smithova diagramu"""
    sigma_a = sigma_c * (1 - sigma_m / sigma_u)
    return sigma_m + sigma_a, sigma_m - sigma_a


def gerber_weyrauch_horni_dolni(sigma_m, sigma_c, sigma_u):
    """(3.19)  vrací (σh, σd)"""
    sigma_a = sigma_c * (1 - (sigma_m / sigma_u) ** 2)
    return sigma_m + sigma_a, sigma_m - sigma_a


def Kf(sigma_c, sigma_c_vrub):
    """(3.20)"""
    return sigma_c / sigma_c_vrub


def citlivost_vrubu(Kf, Kt):
    """(3.22)"""
    return (Kf - 1) / (Kt - 1)


def Kf_z_citlivosti(q, Kt):
    """(3.22)"""
    return 1 + q * (Kt - 1)


def Kf_z_gradientu(Kt, C_SS, chi):
    """(3.24)"""
    return Kt / (1 + math.sqrt(C_SS * chi))


def gough_pollard(sigma_a, tau_a, sigma_c, tau_c):
    """(3.25)  levá strana, porušení při ≥ 1"""
    return (sigma_a / sigma_c) ** 2 + (tau_a / tau_c) ** 2


def gough_pollard_krehke(sigma_a, tau_a, sigma_c, tau_c):
    """(3.26)"""
    return sigma_a / sigma_c + (tau_a / tau_c) ** 2


def gough_zobecnene(sigma_a, tau_a, sigma_c, tau_c):
    """(3.27)"""
    k = sigma_c / tau_c
    return (sigma_a / sigma_c) ** k + (tau_a / tau_c) ** 2


def liu_zenner(sigma_a, tau_a, sigma_c, tau_c):
    """(3.28)"""
    k = sigma_c / tau_c
    return (k - 1) * (sigma_a / sigma_c) ** 2 + (2 - k) * (tau_a / tau_c) ** 2


def lee_empiricke(sigma_a, tau_a, sigma_c, tau_c, alpha_np, phi, zeta=2.0):
    """(3.29)  phi v radiánech"""
    return ((sigma_a / sigma_c) ** (zeta * (1 + alpha_np * math.sin(phi)))
            + (tau_a / tau_c) ** (zeta * (1 + alpha_np * math.cos(phi))))


def alpha_np(sigma_hmh_90, sigma_hmh_0):
    """(3.30)"""
    return sigma_hmh_90 / sigma_hmh_0


def sigma_eq_tah_krut(sigma_a, tau_a, sigma_c, tau_c):
    """(3.31)"""
    k = sigma_c / tau_c
    return math.sqrt(sigma_a ** 2 + k ** 2 * tau_a ** 2)


def sigma_eq_rankine(sigma_1a):
    """(3.32)"""
    return sigma_1a


def sigma_eq_guest(sigma_1a, sigma_2a):
    """(3.33)"""
    return sigma_1a - sigma_2a


def alpha_hmh(sigma_c_R_minus1, sigma_c_R0):
    """(3.35)"""
    return (sigma_c_R_minus1 - sigma_c_R0) / sigma_c_R0


def hmh_amplituda(s1, s2, s3):
    """(3.36)"""
    return math.sqrt(0.5 * ((s1 - s2) ** 2 + (s2 - s3) ** 2 + (s1 - s3) ** 2))


def sigma_eq_hmh(sigma_a, sigma_mx, sigma_my, sigma_mz, alpha):
    """(3.34)"""
    return sigma_a + alpha * (sigma_mx + sigma_my + sigma_mz)


def hydrostaticke_napeti(s1, s2, s3):
    """(3.38)"""
    return (s1 + s2 + s3) / 3


def tau_oct(s1, s2, s3):
    """(3.39)"""
    return math.sqrt((s1 - s2) ** 2 + (s2 - s3) ** 2 + (s1 - s3) ** 2) / 3


def sigma_eq_sines(tau_oct_a, alpha_oct, s1m, s2m, s3m):
    """(3.37)"""
    return tau_oct_a + alpha_oct * (s1m + s2m + s3m)


def sigma_eq_lee(sigma_hmh_a, alpha_NP, f_NP):
    """(3.40)"""
    return sigma_hmh_a * (1 + alpha_NP * f_NP)


def dang_van(tau_mic, sigma_hp, a):
    """(3.41)  levá strana, porovnat s b"""
    return tau_mic + a * sigma_hp


def fdp(alpha, beta, sigma, tau):
    """(3.42)"""
    return alpha * sigma + beta * tau


def findley(tau_a, sigma_n_max, k_F=0.2):
    """(3.44)"""
    return tau_a + k_F * sigma_n_max


def findley_basquin(tau_f, b, Nf):
    """(3.45)"""
    return tau_f * (2 * Nf) ** b


def mcdiarmid(tau_a, sigma_n_max, tau_c_AB, sigma_u):
    """(3.46)"""
    return tau_a + tau_c_AB / (2 * sigma_u) * sigma_n_max


def biaxialita(s1, s3):
    """(3.48)"""
    return s3 / s1 if abs(s1) > abs(s3) else s1 / s3


def skalovane_normalove(k, V):
    """(3.47)"""
    return 1 + (1 - k) * V
