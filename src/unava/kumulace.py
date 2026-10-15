import math
from .zaklady import amplituda_napeti, stredni_napeti
from .sn import basquin_cykly

__all__ = [
    "poskozeni", "palmgren_miner", "pocet_bloku", "svum_ekvivalentni", "svum_poskozeni",
    "haibach_cykly", "corten_dolan_cykly_bloku", "corten_dolan_pocet_bloku",
    "corten_dolan_N_l", "corten_dolan_D", "extremy", "preusporadej_od_maxima",
    "stekajici_dest", "reverzni_amplituda", "morrow_cykly", "zivotnost_bloku",
]


def poskozeni(n, N):
    """(4.3)"""
    return n / N


def palmgren_miner(n, N):
    """(4.5)  D = Σ ni/Ni"""
    return sum(ni / Ni for ni, Ni in zip(n, N))


def pocet_bloku(n, N):
    """(4.6)  B = 1 / Σ ni/Ni"""
    return 1 / palmgren_miner(n, N)


def svum_ekvivalentni(sigma_a, sigma_m, k_svum):
    """(4.7)"""
    return sigma_a + k_svum * sigma_m


def svum_poskozeni(sigma_eq, sigma_c, m_svum):
    """(4.8)"""
    return (sigma_eq / sigma_c) ** m_svum


def haibach_cykly(sigma_a, C, m, sigma_c):
    """(4.9), (4.11), (4.12)  pod mezí únavy sklon 2m"""
    if sigma_a >= sigma_c:
        return C / sigma_a ** m
    C2 = C * sigma_c ** m
    return C2 / sigma_a ** (2 * m)


def corten_dolan_cykly_bloku(N1, sigma_a, n, d):
    """(4.13)  Nb; sigma_a[0] je nejvyšší stupeň"""
    if sigma_a[0] != max(sigma_a):
        raise ValueError("sigma_a[0] musí být nejvyšší stupeň zatížení (k němu patří N1)")
    nb = sum(n)
    jmenovatel = sum(ni / nb * (si / sigma_a[0]) ** d for ni, si in zip(n, sigma_a))
    return N1 / jmenovatel


def corten_dolan_pocet_bloku(N1, sigma_a, n, d):
    """(4.13), (4.14)  B = Nb / nb"""
    return corten_dolan_cykly_bloku(N1, sigma_a, n, d) / sum(n)


def corten_dolan_N_l(N1, sigma_a1, sigma_ai, d):
    """(4.19)"""
    return N1 * (sigma_a1 / sigma_ai) ** d


def corten_dolan_D(n, N, alpha):
    """(4.20)"""
    return sum((ni / Ni) ** alpha for ni, Ni in zip(n, N))


def extremy(hodnoty, nazvy=None):
    """redukce historie na vrcholy (Peak – Valley)"""
    if nazvy is None:
        nazvy = list(range(len(hodnoty)))
    h, j = [], []
    for hodnota, nazev in zip(hodnoty, nazvy):
        if h and hodnota == h[-1]:
            continue
        if len(h) >= 2 and (h[-1] - h[-2]) * (hodnota - h[-1]) > 0:
            h[-1], j[-1] = hodnota, nazev
        else:
            h.append(hodnota)
            j.append(nazev)
    return h, j


def preusporadej_od_maxima(hodnoty, nazvy=None):
    """kap. 4.2.1  blok začínající a končící absolutním extrémem"""
    if nazvy is None:
        nazvy = list(range(len(hodnoty)))
    i = max(range(len(hodnoty)), key=lambda k: abs(hodnoty[k]))
    return extremy(list(hodnoty[i:]) + list(hodnoty[:i + 1]),
                   list(nazvy[i:]) + list(nazvy[:i + 1]))


def stekajici_dest(hodnoty, nazvy=None, opakovany_blok=True):
    """(4.21)  třibodová metoda (ASTM E1049); vrací [(od, do, σh, σd, počet)]"""
    if opakovany_blok:
        hodnoty, nazvy = preusporadej_od_maxima(hodnoty, nazvy)
    else:
        hodnoty, nazvy = extremy(hodnoty, nazvy)
    if len(hodnoty) < 2:
        return []
    zasobnik, cykly = [], []
    for bod in zip(hodnoty, nazvy):
        zasobnik.append(bod)
        while len(zasobnik) >= 3:
            Y = abs(zasobnik[-2][0] - zasobnik[-3][0])
            X = abs(zasobnik[-1][0] - zasobnik[-2][0])
            if X < Y or (opakovany_blok and len(zasobnik) == 3):
                break
            (h1, j1), (h2, j2) = zasobnik[-3], zasobnik[-2]
            if len(zasobnik) == 3:
                cykly.append((j1, j2, max(h1, h2), min(h1, h2), 0.5))
                zasobnik.pop(0)
            else:
                cykly.append((j1, j2, max(h1, h2), min(h1, h2), 1.0))
                del zasobnik[-3:-1]
    if opakovany_blok:
        (h1, j1), (h2, j2) = zasobnik[0], zasobnik[1]
        cykly.append((j1, j2, max(h1, h2), min(h1, h2), 1.0))
    else:
        for (h1, j1), (h2, j2) in zip(zasobnik, zasobnik[1:]):
            cykly.append((j1, j2, max(h1, h2), min(h1, h2), 0.5))
    return cykly


def reverzni_amplituda(sigma_a, sigma_m, sigma_f):
    """(4.25)  Morrow"""
    if sigma_m >= sigma_f:
        raise ValueError("střední napětí musí být menší než sigma_f")
    return sigma_a / (1 - sigma_m / sigma_f)


def morrow_cykly(sigma_ar, sigma_f, b):
    """(4.24)  σar = σf·(2Nf)^b"""
    return basquin_cykly(sigma_ar, sigma_f, b)


def zivotnost_bloku(hodnoty, sigma_f, b, nazvy=None, opakovany_blok=True):
    """(4.22)–(4.26)  vrací (B, tabulka cyklů)"""
    tabulka, D = [], 0
    for od, do, sh, sd, pocet in stekajici_dest(hodnoty, nazvy, opakovany_blok):
        sa = amplituda_napeti(sh, sd)
        sm = stredni_napeti(sh, sd)
        sar = reverzni_amplituda(sa, sm, sigma_f)
        Nf = morrow_cykly(sar, sigma_f, b)
        D += pocet / Nf
        tabulka.append({"cyklus": f"{od}–{do}", "pocet": pocet, "sigma_a": sa,
                        "sigma_m": sm, "sigma_ar": sar, "Nf": Nf})
    return (1 / D if D > 0 else math.inf), tabulka
