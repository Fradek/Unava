import math

from ._numerika import _koren, _pulcykly

__all__ = [
    "celkova_deformace", "ramberg_osgood_napeti", "ramberg_osgood_deformace",
    "ramberg_osgood_napeti_z_celkove", "mcb_elasticka", "mcb_plasticka",
    "mcb_deformace", "mcb_cykly", "n_z_mcb", "K_z_mcb", "tranzitni_pulcykly",
    "morrow_deformace", "morrow_deformace_cykly", "swt_parametr", "swt_prava_strana",
    "swt_cykly", "neuber", "glinka", "hmh_amplituda_napeti", "hmh_amplituda_deformace",
    "efektivni_poisson", "smykova_zivotnost", "brown_miller_cykly", "S_kbm",
    "kandil_brown_miller", "wang_brown_A_B", "wang_brown_cykly", "fatemi_socie",
    "fatemi_socie_cykly", "jiang_dY", "jiang_dD", "swt_energeticke_cykly",
]


def celkova_deformace(eps_el, eps_pl):
    """(5.1)"""
    return eps_el + eps_pl


def ramberg_osgood_napeti(eps_a_pl, K, n):
    """(5.2)"""
    return K * eps_a_pl ** n


def ramberg_osgood_deformace(sigma_a, E, K, n):
    """(5.3)"""
    return sigma_a / E + (sigma_a / K) ** (1 / n)


def ramberg_osgood_napeti_z_celkove(eps_a, E, K, n):
    """(5.3)  inverze"""
    if eps_a <= 0:
        raise ValueError("amplituda deformace musí být kladná")
    return _koren(lambda s: ramberg_osgood_deformace(s, E, K, n) - eps_a, 1e-12, E * eps_a)


def mcb_elasticka(sigma_f, b, E, Nf):
    """(5.4)"""
    return sigma_f / E * (2 * Nf) ** b


def mcb_plasticka(eps_f, c, Nf):
    """(5.5)"""
    return eps_f * (2 * Nf) ** c


def mcb_deformace(sigma_f, b, eps_f, c, E, Nf):
    """(5.6)"""
    return mcb_elasticka(sigma_f, b, E, Nf) + mcb_plasticka(eps_f, c, Nf)


def mcb_cykly(eps_a, sigma_f, b, eps_f, c, E):
    """(5.6)  inverze → Nf"""
    return _pulcykly(lambda x: sigma_f / E * x ** b + eps_f * x ** c, eps_a) / 2


def n_z_mcb(b, c):
    """(5.9)"""
    return b / c


def K_z_mcb(sigma_f, eps_f, b, c):
    """(5.10)"""
    return sigma_f * eps_f ** (-b / c)


def tranzitni_pulcykly(sigma_f, b, eps_f, c, E):
    """2NT z rovnosti (5.4) = (5.5)"""
    return (eps_f * E / sigma_f) ** (1 / (b - c))


def morrow_deformace(sigma_f, b, eps_f, c, E, Nf, sigma_m):
    """(5.11)"""
    return (sigma_f - sigma_m) / E * (2 * Nf) ** b + eps_f * (2 * Nf) ** c


def morrow_deformace_cykly(eps_a, sigma_m, sigma_f, b, eps_f, c, E):
    """(5.11)  inverze → Nf"""
    if sigma_m >= sigma_f:
        raise ValueError("střední napětí musí být menší než sigma_f")
    return _pulcykly(lambda x: (sigma_f - sigma_m) / E * x ** b + eps_f * x ** c, eps_a) / 2


def swt_parametr(sigma_max, eps_a):
    """(5.12)"""
    return sigma_max * eps_a


def swt_prava_strana(sigma_f, b, eps_f, c, E, Nf):
    """(5.13)"""
    return sigma_f ** 2 / E * (2 * Nf) ** (2 * b) + sigma_f * eps_f * (2 * Nf) ** (b + c)


def swt_cykly(sigma_max, eps_a, sigma_f, b, eps_f, c, E):
    """(5.13)  inverze → Nf; σmax ≤ 0 → nekonečno"""
    if sigma_max <= 0:
        return math.inf
    return _pulcykly(lambda x: sigma_f ** 2 / E * x ** (2 * b)
                     + sigma_f * eps_f * x ** (b + c), sigma_max * eps_a) / 2


def neuber(K, sigma_EL, E, K_cykl, n_cykl):
    """(5.15)/(5.16)  K = Kt nebo Kf, σEL nominální elastické napětí; vrací (σ, ε)"""
    if sigma_EL == 0:
        return 0.0, 0.0
    znamenko = math.copysign(1, sigma_EL)
    s_el = abs(K * sigma_EL)
    cil = s_el ** 2 / E
    sigma = _koren(lambda s: s * ramberg_osgood_deformace(s, E, K_cykl, n_cykl) - cil,
                   1e-12 * s_el, s_el)
    return znamenko * sigma, znamenko * cil / sigma


def glinka(K, sigma_EL, E, K_cykl, n_cykl):
    """(5.17)–(5.19)  vrací (σ, ε)"""
    U_EL = (K * sigma_EL) ** 2 / (2 * E)

    def U(s):
        return s ** 2 / (2 * E) + s * (s / K_cykl) ** (1 / n_cykl) / (1 + n_cykl)

    if sigma_EL == 0:
        return 0.0, 0.0
    znamenko = math.copysign(1, sigma_EL)
    s_el = abs(K * sigma_EL)
    sigma = _koren(lambda s: U(s) - U_EL, 1e-12 * s_el, s_el)
    return znamenko * sigma, znamenko * ramberg_osgood_deformace(sigma, E, K_cykl, n_cykl)


def hmh_amplituda_napeti(s1a, s2a, s3a):
    """(5.21)"""
    return math.sqrt((s1a - s2a) ** 2 + (s2a - s3a) ** 2 + (s3a - s1a) ** 2) / math.sqrt(2)


def hmh_amplituda_deformace(e1a, e2a, e3a, mu_ef):
    """(5.22)"""
    return (math.sqrt((e1a - e2a) ** 2 + (e2a - e3a) ** 2 + (e3a - e1a) ** 2)
            / (math.sqrt(2) * (1 + mu_ef)))


def efektivni_poisson(eps_el, eps_pl, mu_el=0.3, mu_pl=0.5):
    """(5.23)"""
    return (mu_el * eps_el + mu_pl * eps_pl) / (eps_el + eps_pl)


def smykova_zivotnost(tau_f, b_g, gamma_f, c_g, G, Nf):
    """pravá strana (5.24) a (5.31)"""
    return tau_f / G * (2 * Nf) ** b_g + gamma_f * (2 * Nf) ** c_g


def brown_miller_cykly(dgamma_max, tau_f, b_g, gamma_f, c_g, G):
    """(5.24)  inverze → Nf"""
    return _pulcykly(lambda x: tau_f / G * x ** b_g + gamma_f * x ** c_g, dgamma_max / 2) / 2


def S_kbm(tau_c, sigma_c, mu):
    """(5.27)"""
    return 2 * (tau_c / sigma_c - 0.5 * (1 + mu) / (1 - mu))


def kandil_brown_miller(dgamma_max, deps, S, alpha):
    """(5.26)  levá strana"""
    return ((dgamma_max / 2) ** alpha + S * deps ** alpha) ** (1 / alpha)


def wang_brown_A_B(S):
    """(5.29), (5.30)"""
    return 1.3 + 0.7 * S, 1.5 + 0.5 * S


def wang_brown_cykly(dgamma_max, deps, sigma_m, S, sigma_f, b, eps_f, c, E):
    """(5.28)  inverze → Nf"""
    A, B = wang_brown_A_B(S)
    leva = dgamma_max / 2 + S * deps
    return _pulcykly(lambda x: A * (sigma_f - 2 * sigma_m) / E * x ** b
                     + B * eps_f * x ** c, leva) / 2


def fatemi_socie(dgamma, sigma_max, sigma_Y, k_FS):
    """(5.31)  levá strana"""
    return dgamma / 2 * (1 + k_FS * sigma_max / sigma_Y)


def fatemi_socie_cykly(dgamma, sigma_max, sigma_Y, k_FS, tau_f, b_g, gamma_f, c_g, G):
    """(5.31)  inverze → Nf"""
    return _pulcykly(lambda x: tau_f / G * x ** b_g + gamma_f * x ** c_g,
                     fatemi_socie(dgamma, sigma_max, sigma_Y, k_FS)) / 2


def jiang_dY(a, sigma, deps_p, tau, dgamma_p):
    """(5.33)"""
    return a * sigma * deps_p + (1 - a) / 2 * tau * dgamma_p


def jiang_dD(sigma_mr, sigma_c, m, sigma, sigma_f, dY):
    """(5.32)"""
    return max(sigma_mr - sigma_c, 0) ** m * (1 + sigma / sigma_f) * dY


def swt_energeticke_cykly(sigma_max, deps_1, sigma_f, b, eps_f, c, E):
    """(5.34)  inverze → Nf"""
    return _pulcykly(lambda x: sigma_f * x ** b * (sigma_f / E * x ** b + eps_f * x ** c),
                     sigma_max * deps_1 / 2) / 2
