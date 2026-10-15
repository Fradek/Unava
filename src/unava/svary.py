import math

__all__ = [
    "fat_C", "fat_cykly", "iiw_omezeni", "delta_K_svar", "paris_rychlost",
    "paris_cykly", "korekce_tloustky", "km_efektivni", "km_kombinovane",
]


def fat_C(FAT, m=3):
    """konstanta C z FAT třídy (únavová pevnost při 2·10⁶ cyklech)"""
    return FAT ** m * 2e6


def fat_cykly(delta_sigma, FAT, m=3, mez_cyklu=1e7, za_mezi="vodorovne", m2=22):
    """(7.1)  N = C/Δσ^m; za_mezi = 'vodorovne' (CAFL) nebo 'sklon' (mF = m2)"""
    if delta_sigma <= 0:
        return math.inf
    C = fat_C(FAT, m)
    N = C / delta_sigma ** m
    if N <= mez_cyklu:
        return N
    if za_mezi == "vodorovne":
        return math.inf
    delta_sigma_mez = (C / mez_cyklu) ** (1 / m)
    return mez_cyklu * (delta_sigma_mez / delta_sigma) ** m2


def iiw_omezeni(delta_sigma_n, sigma_n_max, sigma_Y):
    """(7.2), (7.3)"""
    return delta_sigma_n <= 1.5 * sigma_Y and sigma_n_max <= sigma_Y


def delta_K_svar(a, d_mem, Ym, Mkm, d_b, Yb, Mkb):
    """(7.4)"""
    return math.sqrt(math.pi * a) * (d_mem * Ym * Mkm + d_b * Yb * Mkb)


def paris_rychlost(C, delta_K, m):
    """(7.5)"""
    return C * delta_K ** m


def paris_cykly(C, m, Y, delta_sigma, a0, af):
    """(7.6)–(7.8)"""
    if not 0 < a0 < af:
        raise ValueError("musí platit 0 < a0 < af")
    zaklad = C * (Y * delta_sigma * math.sqrt(math.pi)) ** m
    if m == 2:
        return math.log(af / a0) / zaklad
    e = (2 - m) / 2
    return 2 / ((2 - m) * zaklad) * (af ** e - a0 ** e)


def korekce_tloustky(t_eff, n, t_ref=25):
    """(7.9)"""
    return (t_ref / t_eff) ** n


def km_efektivni(km_vypoctene, km_zahrnute):
    """(7.10)"""
    return km_vypoctene / km_zahrnute


def km_kombinovane(km_osove, km_uhlove):
    """(7.11)"""
    return 1 + (km_osove - 1) + (km_uhlove - 1)
