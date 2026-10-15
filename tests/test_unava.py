import math

import pytest

import unava

HISTORIE = [170, -270, 390, 130, 260, -160, 130, -420, 180, 60, 320, -90]
NAZVY = list("ABCDEFGHIJKL")

E = 206000
SIGMA_F, B, EPS_F, C = 742.3, -0.0658, 0.0603, -0.3085
K_CYKL, N_CYKL = 1370, 0.216


def test_zakladni_pojmy():
    assert unava.stredni_napeti(300, -100) == 100
    assert unava.amplituda_napeti(300, -100) == 200
    assert unava.rozkmit_napeti(300, -100) == 400
    assert unava.soucinitel_asymetrie(300, -150) == -0.5
    assert unava.amplituda_z_R(300, -0.5) == pytest.approx(225)
    assert unava.stredni_z_R(300, -0.5) == pytest.approx(75)


def test_basquin_inverze():
    Nf = unava.basquin_cykly(400, 1020, -0.08432)
    assert unava.basquin_napeti(1020, -0.08432, Nf) == pytest.approx(400)
    assert unava.basquin_cykly(0, 1020, -0.08432) == math.inf


def test_basquin_regrese():
    Nf = [1e3, 1e4, 1e5, 1e6]
    sa = [unava.basquin_napeti(900, -0.09, n) for n in Nf]
    sigma_f, b = unava.basquin_regrese(sa, Nf)
    assert sigma_f == pytest.approx(900)
    assert b == pytest.approx(-0.09)


def test_haigh():
    assert unava.haigh_dovolena_amplituda(100, 300, "goodman", sigma_u=600) == pytest.approx(250)
    assert unava.haigh_ekvivalentni_amplituda(250, 100, "goodman", sigma_u=600) == pytest.approx(300)
    assert unava.haigh_dovolena_amplituda(-100, 300, "goodman", sigma_u=600) == 300
    assert unava.haigh_dovolena_amplituda(700, 300, "goodman", sigma_u=600) == 0
    a = unava.haigh_dovolena_amplituda(-100, 300, "clemson", sigma_Y=500, n=1.5,
                                       tlak_bez_prinosu=False)
    assert isinstance(a, float) and a > 300
    with pytest.raises(ValueError):
        unava.haigh_dovolena_amplituda(100, 300, "neexistuje", sigma_u=600)


def test_haibach_spojitost():
    C, m, sc = 1e12, 3, 150
    assert unava.haibach_cykly(sc, C, m, sc) == pytest.approx(
        unava.haibach_cykly(sc * (1 - 1e-12), C, m, sc), rel=1e-9)


def test_corten_dolan():
    with pytest.raises(ValueError):
        unava.corten_dolan_cykly_bloku(1e4, [200, 300], [10, 10], 2.5)
    assert unava.corten_dolan_pocet_bloku(1e4, [300], [1], 2.5) == pytest.approx(1e4)


def test_extremy_plato():
    assert unava.extremy([0, 5, 5, 0, 3])[0] == [0, 5, 0, 3]
    assert unava.extremy([0, 2, 5, 3, 1, 4])[0] == [0, 5, 1, 4]


def test_stekajici_dest_skripta():
    cykly = unava.stekajici_dest(HISTORIE, NAZVY)
    nalezene = {(od, do): (sh - sd, (sh + sd) / 2) for od, do, sh, sd, _ in cykly}
    assert nalezene == {
        ("I", "J"): (120, 120), ("L", "A"): (260, 40), ("K", "B"): (590, 25),
        ("D", "E"): (130, 195), ("F", "G"): (290, -15), ("H", "C"): (810, -15),
    }
    assert all(p == 1.0 for *_, p in cykly)


def test_stekajici_dest_astm():
    cykly = unava.stekajici_dest([-2, 1, -3, 5, -1, 3, -4, 4, -2], opakovany_blok=False)
    vysledek = sorted((sh - sd, p) for _, _, sh, sd, p in cykly)
    assert vysledek == [(3, 0.5), (4, 0.5), (4, 1.0), (6, 0.5), (8, 0.5), (8, 0.5), (9, 0.5)]


def test_stekajici_dest_kratka_historie():
    assert unava.stekajici_dest([5, 5, 5]) == []


def test_zivotnost_bloku_skripta():
    B_bloku, tabulka = unava.zivotnost_bloku(HISTORIE, 1020, -0.08432)
    assert round(B_bloku) == 32787
    assert len(tabulka) == 6


def test_zivotnost_bloku_jeden_pruchod():
    B_bloku, _ = unava.zivotnost_bloku(HISTORIE, 1020, -0.08432, opakovany_blok=False)
    assert round(B_bloku) == 46975


def test_reverzni_amplituda_mez():
    with pytest.raises(ValueError):
        unava.reverzni_amplituda(100, 1020, 1020)


def test_mcb_inverze():
    for eps_a in (0.02, 0.005, 0.002):
        Nf = unava.mcb_cykly(eps_a, SIGMA_F, B, EPS_F, C, E)
        assert unava.mcb_deformace(SIGMA_F, B, EPS_F, C, E, Nf) == pytest.approx(eps_a)


def test_mcb_mimo_rozsah():
    assert unava.mcb_cykly(1e-6, SIGMA_F, B, EPS_F, C, E) == math.inf
    with pytest.raises(ValueError):
        unava.mcb_cykly(5, SIGMA_F, B, EPS_F, C, E)


def test_tranzitni_bod():
    N2t = unava.tranzitni_pulcykly(SIGMA_F, B, EPS_F, C, E)
    assert unava.mcb_elasticka(SIGMA_F, B, E, N2t / 2) == pytest.approx(
        unava.mcb_plasticka(EPS_F, C, N2t / 2))


def test_kompatibilita():
    n = unava.n_z_mcb(B, C)
    K = unava.K_z_mcb(SIGMA_F, EPS_F, B, C)
    Nf = 1000
    sa = unava.basquin_napeti(SIGMA_F, B, Nf)
    assert unava.ramberg_osgood_napeti(unava.mcb_plasticka(EPS_F, C, Nf), K, n) == pytest.approx(sa)


def test_ramberg_osgood_inverze():
    eps = unava.ramberg_osgood_deformace(400, E, K_CYKL, N_CYKL)
    assert unava.ramberg_osgood_napeti_z_celkove(eps, E, K_CYKL, N_CYKL) == pytest.approx(400)


def test_morrow_a_swt_inverze():
    Nf = unava.morrow_deformace_cykly(0.004, 50, SIGMA_F, B, EPS_F, C, E)
    assert unava.morrow_deformace(SIGMA_F, B, EPS_F, C, E, Nf, 50) == pytest.approx(0.004)
    Nf = unava.swt_cykly(400, 0.004, SIGMA_F, B, EPS_F, C, E)
    assert unava.swt_prava_strana(SIGMA_F, B, EPS_F, C, E, Nf) == pytest.approx(1.6)
    assert unava.swt_cykly(-10, 0.004, SIGMA_F, B, EPS_F, C, E) == math.inf


def test_neuber_glinka():
    s_n, e_n = unava.neuber(2.5, 200, E, K_CYKL, N_CYKL)
    assert s_n * e_n == pytest.approx((2.5 * 200) ** 2 / E)
    assert e_n == pytest.approx(unava.ramberg_osgood_deformace(s_n, E, K_CYKL, N_CYKL))
    s_g, _ = unava.glinka(2.5, 200, E, K_CYKL, N_CYKL)
    assert s_g < s_n < 2.5 * 200
    s_t, e_t = unava.neuber(2.5, -200, E, K_CYKL, N_CYKL)
    assert (s_t, e_t) == pytest.approx((-s_n, -e_n))


def test_lomova_mechanika():
    KI = unava.K_I(100, 0.01)
    assert KI == pytest.approx(100 * math.sqrt(math.pi * 0.01))
    E_ef = unava.E_efektivni(E, 0.3)
    assert unava.G_I(KI, E, 0.3) == pytest.approx(unava.J_z_K(KI, E_ef))
    assert unava.hnaci_sila(100, 0.01, E) == pytest.approx(unava.J_z_K(KI, E))
    gamma = 1e-3
    sc = unava.griffith_kriticke_napeti(E_ef, gamma, 0.01)
    assert unava.griffith_kriticka_delka(E_ef, gamma, sc) == pytest.approx(0.01)
    assert unava.plasticka_zona_I(KI, 500, 0) == pytest.approx(unava.rp_I(KI, 500) / 2)
    sx, sy, txy = unava.pole_mod_I(KI, 1e-3, 0)
    assert sx == pytest.approx(sy) and txy == pytest.approx(0)


def test_svary():
    assert unava.fat_cykly(90, 90) == pytest.approx(2e6)
    assert unava.fat_cykly(10, 90) == math.inf
    assert unava.fat_cykly(10, 90, za_mezi="sklon") > 1e7
    assert unava.fat_cykly(0, 90) == math.inf
    N2 = unava.paris_cykly(1e-12, 2, 1, 100, 1e-3, 2e-2)
    N2b = unava.paris_cykly(1e-12, 2.000001, 1, 100, 1e-3, 2e-2)
    assert N2 == pytest.approx(N2b, rel=1e-4)
    with pytest.raises(ValueError):
        unava.paris_cykly(1e-12, 3, 1, 100, 2e-2, 1e-3)


def test_jmenny_prostor():
    assert "math" not in dir(unava)
    assert unava.__version__
