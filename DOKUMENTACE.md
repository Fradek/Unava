# Dokumentace knihovny `unava`

Knihovna obsahuje vzorce ze skript **Životnost a spolehlivost** (Kořínek, Tošková, Fusek; VŠB-TUO, 2025) pro předmět Únava materiálu. Každá funkce má v docstringu číslo rovnice ze skript, takže se dá snadno dohledat, odkud vzorec pochází.

Verze: 0.1.0 · Python ≥ 3.8 · bez závislostí (jen standardní `math`)

---

## Obsah

1. [Instalace](#1-instalace)
2. [Rychlý start](#2-rychlý-start)
3. [Obecné konvence](#3-obecné-konvence)
4. [Značení parametrů](#4-značení-parametrů)
5. [Kapitola 2 – Základní pojmy (`unava.zaklady`)](#5-kapitola-2--základní-pojmy-unavazaklady)
6. [Kapitola 3 – Vysokocyklová únava S-N (`unava.sn`)](#6-kapitola-3--vysokocyklová-únava-s-n-unavasn)
7. [Kapitola 4 – Kumulace poškození (`unava.kumulace`)](#7-kapitola-4--kumulace-poškození-unavakumulace)
8. [Kapitola 5 – Nízkocyklová únava E-N (`unava.en`)](#8-kapitola-5--nízkocyklová-únava-e-n-unavaen)
9. [Kapitola 6 – Lomová mechanika (`unava.lomova_mechanika`)](#9-kapitola-6--lomová-mechanika-unavalomova_mechanika)
10. [Kapitola 7 – Svarové spoje (`unava.svary`)](#10-kapitola-7--svarové-spoje-unavasvary)
11. [Opravy oproti skriptům](#11-opravy-oproti-skriptům)
12. [Řešené příklady](#12-řešené-příklady)
13. [Testy](#13-testy)

---

## 1. Instalace

Ve složce, kde leží `pyproject.toml`:

```
pip install .
```

Pro vývoj (změny v souborech se projeví bez přeinstalování):

```
pip install -e .
```

Kontrola:

```python
import unava
print(unava.__version__)
```

---

## 2. Rychlý start

```python
import unava

napeti = [170, -270, 390, 130, 260, -160, 130, -420, 180, 60, 320, -90]
B, _ = unava.zivotnost_bloku(napeti, sigma_f=1020, b=-0.08432)
print(f"B = {B:.0f}")          # B = 32787
```

Všechny funkce lze volat přímo z `unava`, nebo po kapitolách:

```python
from unava import kumulace, en
kumulace.palmgren_miner([1, 2], [1e5, 1e6])
en.mcb_cykly(0.005, 742.3, -0.0658, 0.0603, -0.3085, 206000)
```

---

## 3. Obecné konvence

**Jednotky.** Knihovna jednotky nepřevádí. Stačí je dodržovat konzistentně, typicky napětí a moduly v MPa a délky v mm. Pro lomovou mechaniku je obvyklé dosazovat délku trhliny v metrech; K pak vyjde v MPa·√m.

**Cykly vs. půlcykly.** Funkce s `cykly` v názvu vrací vždy **počet cyklů N_f**. Převod z počtu půlcyklů 2N_f dělají interně. Jedinou výjimkou je `tranzitni_pulcykly`, která vrací 2N_T.

**Exponenty.** `b` a `c` se zadávají záporně, tak jak jsou ve skriptech (např. `b = -0.08432`).

**Úhly.** Všechny úhly jsou v **radiánech** (`math.radians(90)`).

**Kritéria typu „levá strana = 1“.** Funkce vrací hodnotu levé strany. K porušení dochází, když je výsledek ≥ 1, případně když přesáhne mez uvedenou u funkce.

**Numerické řešení.** Implicitní rovnice (MCB, SWT, Neuber, Glinka, inverze Ramberg-Osgooda apod.) řeší knihovna bisekcí s přesností na úrovni strojové přesnosti. Životnost hledá v rozsahu 2N_f od 10⁻³ do 10³⁰. Pokud je zatížení tak malé, že by životnost přesáhla horní mez, vrací funkce `math.inf`. Pokud je tak velké, že by byla pod dolní mezí, vyhodí `ValueError`.

**Nekonečná životnost.** Hodnotu `math.inf` vrací i v dalších případech: nulová amplituda (`basquin_cykly`, `morrow_cykly`), SWT při σmax ≤ 0, FAT křivka pod mezí únavy a `zivotnost_bloku` pro historii bez cyklů.

**Chybové stavy.** `ValueError` vyhodí knihovna při vstupech, pro které vzorec nemá smysl:
- σm ≥ σ′f u Morrowovy korekce (`reverzni_amplituda`, `morrow_deformace_cykly`),
- neznámá `metoda` u Haighova diagramu,
- u `vg_nemec` zadané jen jedno z `c`, `d`,
- `sigma_a[0]` není nejvyšší stupeň u Corten-Dolana,
- neplatí 0 < a0 < af u `paris_cykly`,
- nekladná amplituda deformace u `ramberg_osgood_napeti_z_celkove`.

**Nevyužité parametry.** Volitelné parametry s hodnotou `None` se zadávají jen tehdy, když je zvolená varianta vzorce potřebuje (např. `sigma_u` pro Goodmana, `sigma_Y` pro Soderberga).

---

## 4. Značení parametrů

| Parametr | Symbol | Význam |
|---|---|---|
| `sigma_h`, `sigma_d` | σh, σd | horní a dolní napětí cyklu |
| `sigma_a`, `sigma_m` | σa, σm | amplituda a střední napětí |
| `sigma_ar` | σar | reverzní (ekvivalentní) amplituda napětí |
| `sigma_max` | σmax | maximální napětí |
| `sigma_n` | σn | nominální napětí |
| `sigma_c` | σc | mez únavy (R = −1) |
| `sigma_hc` | σhc | mez únavy při míjivém cyklu (R = 0) |
| `sigma_u` | σu | mez pevnosti |
| `sigma_Y` | σY | mez kluzu |
| `sigma_f` | σ′f | součinitel únavové pevnosti |
| `b` | b | exponent únavové pevnosti |
| `eps_f` | ε′f | součinitel únavové tažnosti |
| `c` | c | exponent únavové tažnosti |
| `K`, `n` / `K_cykl`, `n_cykl` | K′, n′ | součinitel a exponent cyklického zpevnění |
| `E`, `G` | E, G | modul pružnosti v tahu a ve smyku |
| `mu` | µ | Poissonovo číslo |
| `Nf` | Nf | počet cyklů do porušení |
| `R` | R | součinitel asymetrie cyklu |
| `Kt`, `Kf` | Kt, Kf | teoretický součinitel koncentrace napětí, únavový součinitel vrubu |
| `tau_a`, `tau_c`, `tau_f` | τa, τc, τ′f | smyková amplituda, mez únavy ve smyku, součinitel únavové pevnosti ve smyku |
| `gamma_f`, `b_g`, `c_g` | γ′f, bγ, cγ | smykové parametry E-N křivky |
| `a` | a | délka trhliny (polovina u průchozí trhliny délky 2a) |
| `Y` | Y | korekční (tvarový) faktor |
| `KI`, `KII`, `KIII` | K | součinitele intenzity napětí pro módy I–III |
| `E_ef` | E* | efektivní modul (rovinná napjatost/deformace), viz `E_efektivni` |

---

## 5. Kapitola 2 – Základní pojmy (`unava.zaklady`)

| Funkce | Rov. | Vzorec | Vrací |
|---|---|---|---|
| `stredni_napeti(sigma_h, sigma_d)` | 2.1 | σm = (σh + σd)/2 | σm |
| `amplituda_napeti(sigma_h, sigma_d)` | 2.2 | σa = (σh − σd)/2 | σa |
| `rozkmit_napeti(sigma_h, sigma_d)` | 2.3 | Δσ = σh − σd | Δσ |
| `soucinitel_asymetrie(sigma_h, sigma_d)` | 2.4 | R = σd/σh | R |
| `amplituda_z_R(sigma_h, R)` | 2.5 | σa = σh/2·(1 − R) | σa |
| `stredni_z_R(sigma_h, R)` | 2.5 | σm = σh/2·(1 + R) | σm |
| `frekvence(T)` | 2.6 | f = 1/T | f |
| `Kt(sigma_max, sigma_n)` | 2.7 | Kt = σmax/σn | Kt |
| `pomerny_gradient(dsigma_dx, sigma_max)` | 2.8 | χσ = \|dσ/dx\|/σmax | χσ |
| `bezrozmerny_gradient(chi, r)` | 2.9 | Cσ = χσ·r | Cσ |

---

## 6. Kapitola 3 – Vysokocyklová únava S-N (`unava.sn`)

### 6.1 Wöhlerova křivka

| Funkce | Rov. | Vzorec | Vrací |
|---|---|---|---|
| `basquin_napeti(sigma_f, b, Nf)` | 3.1 | σa = σf·(2Nf)^b | σa |
| `basquin_cykly(sigma_a, sigma_f, b)` | 3.1 | Nf = ½·(σa/σf)^(1/b) | Nf |
| `basquin_regrese(sigma_a, Nf)` | 3.2 | lineární regrese v log-log | `(sigma_f, b)` |

`basquin_regrese` bere seznamy naměřených amplitud a počtů cyklů.

### 6.2 Modifikace S-N křivky

| Funkce | Rov. | Vzorec | Vrací |
|---|---|---|---|
| `modifikovana_mez_unavy(sigma_c, *soucinitele)` | 3.3 | σc× = σc·ηD·ηL·vg·ηR | σc× |
| `eta_D(D, d)` | 3.4 | ηD = 1 − √(0,02·ln(D/d)) | ηD |
| `ETA_L` | Tab. 3.1 | slovník: `"tah"`, `"tah s ohybem"`, `"ohyb"`, `"krut ocel"`, `"krut litina"` | ηL |
| `ETA_R` | Tab. 3.2 | slovník: pravděpodobnost → ηR (0.5 … 0.999999) | ηR |
| `ohybove_napeti(Mo, d, Jy)` | 3.5 | σo = Mo·d/(2Jy) | σo |
| `gradient_ohyb(d)` | 3.6 | χσ = 2/d | χσ |
| `vg_heywood(c, d)` | 3.7 | vg = 1 + √(2c/d) | vg |
| `vg_nemec(sigma_max, sigma_s, c=None, d=None)` | 3.8, 3.9 | vg = 1 + √((σmax − σs)/σs); s `c`, `d` navíc ·2c/d | vg |

`modifikovana_mez_unavy` přijme libovolný počet součinitelů:

```python
unava.modifikovana_mez_unavy(300, unava.eta_D(50, 10), unava.ETA_L["tah"], unava.ETA_R[0.99])
```

### 6.3 Vliv středního napětí (Haighův, Smithův, Gerber-Weyrauchův diagram)

| Funkce | Rov. | Popis | Vrací |
|---|---|---|---|
| `psi(sigma_c, sigma_hc)` | 3.10 | ψ = 2σc/σhc − 1 | ψ |
| `sigma_F(sigma_c, sigma_hc)` | 3.10 | σF = σc/ψ | σF |
| `haigh_dovolena_amplituda(sigma_m, sigma_c, metoda, sigma_u=None, sigma_Y=None, n=None, k=None, tlak_bez_prinosu=True)` | 3.11–3.17 | mezní amplituda pro dané σm | σa |
| `haigh_ekvivalentni_amplituda(sigma_a, sigma_m, metoda, sigma_u=None, sigma_Y=None, n=None, k=None, tlak_bez_prinosu=False)` | 3.11–3.17 | přepočet na σm = 0 (kap. 3.3.2) | σar |
| `smith_horni_dolni(sigma_m, sigma_c, sigma_u)` | 3.18 | mez Smithova diagramu | `(σh, σd)` |
| `gerber_weyrauch_horni_dolni(sigma_m, sigma_c, sigma_u)` | 3.19 | mez Gerber-Weyrauchova diagramu | `(σh, σd)` |

Parametr `metoda` a potřebné vstupy:

| `metoda` | Rov. | Aproximace | Potřebuje |
|---|---|---|---|
| `"goodman"` | 3.11 | σa/σc + σm/σu = 1 | `sigma_u` |
| `"soderberg"` | 3.12 | σa/σc + σm/σY = 1 | `sigma_Y` |
| `"gerber"` | 3.13 | σa/σc + (σm/σu)² = 1 | `sigma_u` |
| `"asme"` | 3.14 | (σa/σc)² + (σm/σu)² = 1 | `sigma_u` |
| `"bagci"` | 3.15 | σa/σc + (σm/σY)⁴ = 1 | `sigma_Y` |
| `"clemson"` | 3.16 | σa/σc + (σm/σY)ⁿ = 1 | `sigma_Y`, `n` |
| `"sekercioglu"` | 3.17 | σa/σc + (σm/σu)ᵏ = 1 | `sigma_u`, `k` (ocel ≈ 0,8; Al ≈ 0,45) |

Parametr `tlak_bez_prinosu=True` (výchozí u `haigh_dovolena_amplituda`) nahradí tlakové σm nulou. Odpovídá to konzervativnímu předpokladu z obr. 3.8, že tlakové střední napětí únavovou pevnost nezvyšuje. U `haigh_ekvivalentni_amplituda` je výchozí `False`, takže tlakové σm životnost prodlužuje, jak ukazuje obr. 3.9.

Chování na okrajích:
- Pokud σm dosáhne meze (σu, resp. σY), vrací `haigh_dovolena_amplituda` 0 a `haigh_ekvivalentni_amplituda` `math.inf`.
- Gerber, ASME a Bagci mají sudou mocninu, takže jsou vůči znaménku σm symetrické a tlak u nich životnost nezlepší.
- U Clemsona a Sekercioglu s neceločíselným exponentem se pro tlakové σm použije −(|σm|/σ)ⁿ, aby výsledek zůstal reálný. Tlak se pak chová podobně jako u Goodmana.

```python
unava.haigh_dovolena_amplituda(100, 300, "goodman", sigma_u=600)          # 250.0
unava.haigh_ekvivalentni_amplituda(250, 100, "goodman", sigma_u=600)      # 300.0
```

### 6.4 Vruby

| Funkce | Rov. | Vzorec | Vrací |
|---|---|---|---|
| `Kf(sigma_c, sigma_c_vrub)` | 3.20 | Kf = σc/σc* | Kf |
| `citlivost_vrubu(Kf, Kt)` | 3.22 | q = (Kf − 1)/(Kt − 1) | q |
| `Kf_z_citlivosti(q, Kt)` | 3.22 | Kf = 1 + q·(Kt − 1) | Kf |
| `Kf_z_gradientu(Kt, C_SS, chi)` | 3.24 | Kf = Kt/(1 + √(CSS·χσ)) | Kf |

### 6.5 Multiaxiální zatěžování – empirické přístupy

Všechny funkce vrací levou stranu kritéria. Porušení nastává při hodnotě ≥ 1.

| Funkce | Rov. | Vzorec |
|---|---|---|
| `gough_pollard(sigma_a, tau_a, sigma_c, tau_c)` | 3.25 | (σa/σc)² + (τa/τc)² |
| `gough_pollard_krehke(sigma_a, tau_a, sigma_c, tau_c)` | 3.26 | σa/σc + (τa/τc)² |
| `gough_zobecnene(sigma_a, tau_a, sigma_c, tau_c)` | 3.27 | (σa/σc)^k + (τa/τc)², k = σc/τc |
| `liu_zenner(sigma_a, tau_a, sigma_c, tau_c)` | 3.28 | (k−1)(σa/σc)² + (2−k)(τa/τc)² |
| `lee_empiricke(sigma_a, tau_a, sigma_c, tau_c, alpha_np, phi, zeta=2.0)` | 3.29 | Leeovo kritérium; ζ = 2 houževnaté, 1,5 křehké; φ v rad |
| `alpha_np(sigma_hmh_90, sigma_hmh_0)` | 3.30 | αnp = σa,HMH(φ=90°)/σa,HMH(φ=0°) |

### 6.6 Multiaxiální zatěžování – ekvivalentní napětí

Výsledek se porovnává s σc pro trvalou pevnost, případně s Basquinem (3.1) pro časovanou pevnost.

| Funkce | Rov. | Vzorec |
|---|---|---|
| `sigma_eq_tah_krut(sigma_a, tau_a, sigma_c, tau_c)` | 3.31 | √(σa² + k²τa²) |
| `sigma_eq_rankine(sigma_1a)` | 3.32 | σ1a |
| `sigma_eq_guest(sigma_1a, sigma_2a)` | 3.33 | σ1a − σ2a (podle skript porovnat s τc) |
| `hmh_amplituda(s1, s2, s3)` | 3.36 | √(½·[(σ1−σ2)² + (σ2−σ3)² + (σ1−σ3)²]) |
| `alpha_hmh(sigma_c_R_minus1, sigma_c_R0)` | 3.35 | (σc,R=−1 − σc,R=0)/σc,R=0 |
| `sigma_eq_hmh(sigma_a, sigma_mx, sigma_my, sigma_mz, alpha)` | 3.34 | σa + α·(σmx + σmy + σmz) |
| `hydrostaticke_napeti(s1, s2, s3)` | 3.38 | (σ1 + σ2 + σ3)/3 |
| `tau_oct(s1, s2, s3)` | 3.39 | ⅓·√Σ(σi − σj)² |
| `sigma_eq_sines(tau_oct_a, alpha_oct, s1m, s2m, s3m)` | 3.37 | τoct,a + αoct·(σ1m + σ2m + σ3m) |
| `sigma_eq_lee(sigma_hmh_a, alpha_NP, f_NP)` | 3.40 | σHMH,a·(1 + αNP·fNP) |
| `dang_van(tau_mic, sigma_hp, a)` | 3.41 | τmic + a·σhp (porovnat s b) |

Pro Sinesovo pravidlo platí τoct,a = `tau_oct(...)/2`.

### 6.7 Metoda kritické roviny

| Funkce | Rov. | Vzorec |
|---|---|---|
| `fdp(alpha, beta, sigma, tau)` | 3.42 | FDP = ασ + βτ |
| `findley(tau_a, sigma_n_max, k_F=0.2)` | 3.44 | τa + kF·σn,max |
| `findley_basquin(tau_f, b, Nf)` | 3.45 | τ′f·(2Nf)^b |
| `mcdiarmid(tau_a, sigma_n_max, tau_c_AB, sigma_u)` | 3.46 | τa + τc(A,B)/(2σu)·σn,max |
| `skalovane_normalove(k, V)` | 3.47 | f = 1 + (1 − k)·V |
| `biaxialita(s1, s3)` | 3.48 | V = σ3/σ1 pro \|σ1\| > \|σ3\|, jinak σ1/σ3 |

---

## 7. Kapitola 4 – Kumulace poškození (`unava.kumulace`)

### 7.1 Hypotézy kumulace

| Funkce | Rov. | Vzorec | Vrací |
|---|---|---|---|
| `poskozeni(n, N)` | 4.3 | Di = ni/Ni | Di |
| `palmgren_miner(n, N)` | 4.5 | D = Σ ni/Ni | D |
| `pocet_bloku(n, N)` | 4.6 | B = 1/Σ(ni/Ni) | B |
| `svum_ekvivalentni(sigma_a, sigma_m, k_svum)` | 4.7 | σeq = σa + kSVÚM·σm | σeq |
| `svum_poskozeni(sigma_eq, sigma_c, m_svum)` | 4.8 | Di = (σeq/σc)^m | Di |
| `haibach_cykly(sigma_a, C, m, sigma_c)` | 4.9–4.12 | N = C/σa^m pro σa ≥ σc, jinak C′/σa^(2m) | N |
| `corten_dolan_cykly_bloku(N1, sigma_a, n, d)` | 4.13 | Nb = N1/Σ vi·(σai/σa1)^d | Nb |
| `corten_dolan_pocet_bloku(N1, sigma_a, n, d)` | 4.13, 4.14 | B = Nb/nb | B |
| `corten_dolan_N_l(N1, sigma_a1, sigma_ai, d)` | 4.19 | N̄l = N1·(σa1/σai)^d | N̄l |
| `corten_dolan_D(n, N, alpha)` | 4.20 | D = Σ(ni/Ni)^α | D |

`n` a `N` jsou seznamy stejné délky. U Corten-Dolana musí být `sigma_a[0]` nejvyšší stupeň zatížení a `N1` počet cyklů do porušení na této hladině. Exponent `d` bývá 0,7–0,9·m.

### 7.2 Metoda stékajícího deště

| Funkce | Popis | Vrací |
|---|---|---|
| `extremy(hodnoty, nazvy=None)` | redukce historie na vrcholy (Peak – Valley) | `(hodnoty, nazvy)` |
| `preusporadej_od_maxima(hodnoty, nazvy=None)` | blok začínající a končící absolutním extrémem (kap. 4.2.1) | `(hodnoty, nazvy)` |
| `stekajici_dest(hodnoty, nazvy=None, opakovany_blok=True)` | třibodová metoda podle ASTM E1049, rov. 4.21 | seznam cyklů |
| `reverzni_amplituda(sigma_a, sigma_m, sigma_f)` | 4.25: σar = σa/(1 − σm/σf) | σar |
| `morrow_cykly(sigma_ar, sigma_f, b)` | 4.24: σar = σf·(2Nf)^b | Nf |
| `zivotnost_bloku(hodnoty, sigma_f, b, nazvy=None, opakovany_blok=True)` | celý postup 4.22–4.26 | `(B, tabulka)` |

`extremy` nejdřív odstraní opakované hodnoty (plató) a body uvnitř monotónních úseků. U plató ponechá název prvního bodu.

`stekajici_dest` vrací seznam n-tic `(od, do, σh, σd, pocet)`. Historie s méně než dvěma různými hodnotami dá prázdný seznam. `od` a `do` jsou názvy vrcholů, a pokud nejsou zadané, jejich indexy. `pocet` je 1.0 pro celý cyklus a 0.5 pro půlcyklus.

Parametr `opakovany_blok`:
- `True` (výchozí) odpovídá postupu ve skriptech: historie se přeuspořádá od absolutního extrému a vycházejí jen celé cykly.
- `False` znamená jeden průchod historií; zbytek se počítá jako půlcykly.

`zivotnost_bloku` vrací B (počet opakování bloku do porušení) a seznam slovníků s klíči `cyklus`, `pocet`, `sigma_a`, `sigma_m`, `sigma_ar`, `Nf`. Celkový počet cyklů do porušení je B × počet cyklů v bloku.

```python
B, tab = unava.zivotnost_bloku(napeti, 1020, -0.08432, nazvy=list("ABCDEFGHIJKL"))
celkem_cyklu = B * sum(c["pocet"] for c in tab)
```

---

## 8. Kapitola 5 – Nízkocyklová únava E-N (`unava.en`)

### 8.1 Cyklická deformační křivka a Manson–Coffin–Basquin

| Funkce | Rov. | Vzorec | Vrací |
|---|---|---|---|
| `celkova_deformace(eps_el, eps_pl)` | 5.1 | ε = εel + εpl | ε |
| `ramberg_osgood_napeti(eps_a_pl, K, n)` | 5.2 | σa = K′·εa,pl^n′ | σa |
| `ramberg_osgood_deformace(sigma_a, E, K, n)` | 5.3 | εa = σa/E + (σa/K′)^(1/n′) | εa |
| `ramberg_osgood_napeti_z_celkove(eps_a, E, K, n)` | 5.3 | numerická inverze | σa |
| `mcb_elasticka(sigma_f, b, E, Nf)` | 5.4 | σ′f/E·(2Nf)^b | εa,el |
| `mcb_plasticka(eps_f, c, Nf)` | 5.5 | ε′f·(2Nf)^c | εa,pl |
| `mcb_deformace(sigma_f, b, eps_f, c, E, Nf)` | 5.6 | součet 5.4 + 5.5 | εa |
| `mcb_cykly(eps_a, sigma_f, b, eps_f, c, E)` | 5.6 | numerická inverze | Nf |
| `n_z_mcb(b, c)` | 5.9 | n′ = b/c | n′ |
| `K_z_mcb(sigma_f, eps_f, b, c)` | 5.10 | K′ = σ′f·ε′f^(−b/c) | K′ |
| `tranzitni_pulcykly(sigma_f, b, eps_f, c, E)` | 5.4 = 5.5 | 2NT = (ε′f·E/σ′f)^(1/(b−c)) | 2NT |

### 8.2 Vliv středního napětí

| Funkce | Rov. | Vzorec | Vrací |
|---|---|---|---|
| `morrow_deformace(sigma_f, b, eps_f, c, E, Nf, sigma_m)` | 5.11 | (σ′f − σm)/E·(2Nf)^b + ε′f·(2Nf)^c | εa |
| `morrow_deformace_cykly(eps_a, sigma_m, sigma_f, b, eps_f, c, E)` | 5.11 | numerická inverze | Nf |
| `swt_parametr(sigma_max, eps_a)` | 5.12 | PSWT = σmax·εa | PSWT |
| `swt_prava_strana(sigma_f, b, eps_f, c, E, Nf)` | 5.13 | σ′f²/E·(2Nf)^(2b) + σ′f·ε′f·(2Nf)^(b+c) | PSWT |
| `swt_cykly(sigma_max, eps_a, sigma_f, b, eps_f, c, E)` | 5.13 | numerická inverze; pro σmax ≤ 0 vrací `inf` | Nf |

### 8.3 Korekce napětí ve vrubu

| Funkce | Rov. | Popis | Vrací |
|---|---|---|---|
| `neuber(K, sigma_EL, E, K_cykl, n_cykl)` | 5.15, 5.16 | K²·σEL·εEL = σ·ε s Ramberg-Osgoodovou křivkou | `(σ, ε)` |
| `glinka(K, sigma_EL, E, K_cykl, n_cykl)` | 5.17–5.19 | rovnost energií, integrál pod R-O křivkou | `(σ, ε)` |

`K` je Kt (rov. 5.15) nebo Kf (rov. 5.16) a `sigma_EL` je nominální elastické napětí. Neuber dává vždy vyšší (konzervativnější) lokální napětí než Glinka. Pro záporné `sigma_EL` (tlak) vrací obě funkce symetrické řešení se záporným σ i ε.

```python
sigma, eps = unava.neuber(2.5, 200, 206000, 1370, 0.216)   # σ ≈ 349 MPa
sigma, eps = unava.glinka(2.5, 200, 206000, 1370, 0.216)   # σ ≈ 326 MPa
```

### 8.4 Multiaxiální zatěžování

| Funkce | Rov. | Popis | Vrací |
|---|---|---|---|
| `hmh_amplituda_napeti(s1a, s2a, s3a)` | 5.21 | ekvivalentní amplituda napětí | σa,HMH |
| `hmh_amplituda_deformace(e1a, e2a, e3a, mu_ef)` | 5.22 | ekvivalentní amplituda deformace | εa,HMH |
| `efektivni_poisson(eps_el, eps_pl, mu_el=0.3, mu_pl=0.5)` | 5.23 | µ̄ = (µel·εel + µpl·εpl)/εa | µ̄ |
| `smykova_zivotnost(tau_f, b_g, gamma_f, c_g, G, Nf)` | 5.24, 5.31 | τ′f/G·(2Nf)^bγ + γ′f·(2Nf)^cγ | pravá strana |
| `brown_miller_cykly(dgamma_max, tau_f, b_g, gamma_f, c_g, G)` | 5.24 | inverze Brown-Millera | Nf |
| `S_kbm(tau_c, sigma_c, mu)` | 5.27 | S = 2·(τc/σc − 0,5·(1+µ)/(1−µ)) | S |
| `kandil_brown_miller(dgamma_max, deps, S, alpha)` | 5.26 | ((Δγmax/2)^α + S·Δε^α)^(1/α) | levá strana |
| `wang_brown_A_B(S)` | 5.29, 5.30 | A = 1,3 + 0,7S; B = 1,5 + 0,5S | `(A, B)` |
| `wang_brown_cykly(dgamma_max, deps, sigma_m, S, sigma_f, b, eps_f, c, E)` | 5.28 | inverze Wang-Browna | Nf |
| `fatemi_socie(dgamma, sigma_max, sigma_Y, k_FS)` | 5.31 | Δγ/2·(1 + kFS·σmax/σY) | levá strana |
| `fatemi_socie_cykly(dgamma, sigma_max, sigma_Y, k_FS, tau_f, b_g, gamma_f, c_g, G)` | 5.31 | inverze Fatemi-Socie | Nf |
| `jiang_dY(a, sigma, deps_p, tau, dgamma_p)` | 5.33 | dY = aσ·dεp + (1−a)/2·τ·dγp | dY |
| `jiang_dD(sigma_mr, sigma_c, m, sigma, sigma_f, dY)` | 5.32 | dD = ⟨σmr − σc⟩^m·(1 + σ/σf)·dY | dD |
| `swt_energeticke_cykly(sigma_max, deps_1, sigma_f, b, eps_f, c, E)` | 5.34 | inverze energetického SWT | Nf |

---

## 9. Kapitola 6 – Lomová mechanika (`unava.lomova_mechanika`)

U funkcí s parametrem `mu=None` znamená `None` rovinnou napjatost a zadané µ rovinnou deformaci.

### 9.1 Součinitele intenzity napětí a napjatost u čela trhliny

| Funkce | Rov. | Popis | Vrací |
|---|---|---|---|
| `K_I(sigma, a, Y=1.0)` | 6.7 | σ·√(πa)·Y | KI |
| `K_II(tau, a, Y=1.0)` | 6.12 | τ·√(πa)·Y | KII |
| `K_III(tau, a, Y=1.0)` | 6.15 | τ·√(πa)·Y | KIII |
| `sigma_z_rovinna_deformace(sigma_x, sigma_y, mu)` | 6.1 | σz = µ(σx + σy) | σz |
| `kappa(mu, rovinna_deformace=True)` | 6.6 | 3 − 4µ, resp. (3 − µ)/(1 + µ) | κ |
| `pole_mod_I(KI, r, theta, G=None, mu=None, rovinna_deformace=True)` | 6.5 | napětí (a posuvy) u čela, mód I | `(σx, σy, τxy[, u, v])` |
| `pole_mod_II(KII, r, theta, G=None, mu=None, rovinna_deformace=True)` | 6.11 | totéž pro mód II | `(σx, σy, τxy[, u, v])` |
| `pole_mod_III(KIII, r, theta, G=None)` | 6.14 | totéž pro mód III | `(τxz, τyz[, w])` |
| `hlavni_napeti_rovinna(sigma_x, sigma_y, tau_xy)` | 6.16 | hlavní napětí | `(σ1, σ2)` |

Posuvy u, v (w) se vrací jen tehdy, když je zadán modul `G` (u módů I a II také `mu`).

### 9.2 Plastická zóna

| Funkce | Rov. | Popis | Vrací |
|---|---|---|---|
| `plasticka_zona_I(KI, sigma_Y, theta, mu=None)` | 6.20 | radiusvektor hranice, mód I | r |
| `plasticka_zona_II(KII, sigma_Y, theta, mu=None)` | 6.22 | radiusvektor hranice, mód II | r |
| `plasticka_zona_III(KIII, sigma_Y)` | 6.26 | kruhová zóna, mód III | r |
| `rp_I(KI, sigma_Y)` | 6.21 | (1/π)·(KI/σY)² | rp |
| `rp_II(KII, sigma_Y)` | 6.23 | (3/π)·(KII/σY)² | rp |
| `rp_III(KIII, sigma_Y)` | 6.27 | (1/π)·(KIII/τY)² | rp |

### 9.3 Lineární elastická lomová mechanika

| Funkce | Rov. | Popis | Vrací |
|---|---|---|---|
| `mala_plasticka_zona(sigma, sigma_Y)` | 6.28 | σ/σY < 0,3 | `bool` |
| `platnost_KIC(B, a, W, KIC, sigma_Y)` | 6.30 | B, a, (W − a) ≥ 2,5·(KIC/σY)² | `bool` |
| `kriterium_K(KI, KIC, k=2.0)` | 6.31 | KI < KIC/k | `bool` |
| `KIC_barson_rolfe(E, KV)` | 6.33 | 0,0149·√E·KV^(3/4) | KIC |
| `KIC_sailors_corten(E, KV)` | 6.34 | 0,0321·√(E·KV) | KIC |

### 9.4 Energetická kritéria

| Funkce | Rov. | Popis | Vrací |
|---|---|---|---|
| `energie_bez_trhliny(sigma, E, V)` | 6.35 | W0 = σ²/(2E)·V | W0 |
| `E_efektivni(E, mu=None)` | 6.40 | E*, resp. E/(1 − µ²) | E* |
| `griffith_WT(sigma, a, B, E_ef)` | 6.39 | WT = πBσ²a²/E* | WT |
| `plocha_trhliny(a, B, pruchozi=True)` | 6.42 | 4aB, resp. 2aB | Sp |
| `griffith_kriticke_napeti(E_ef, gamma, a)` | 6.43 | σc = √(2E*γ/(πa)) | σc |
| `griffith_kriticka_delka(E_ef, gamma, sigma)` | 6.44 | ac = 2E*γ/(πσ²) | ac |
| `hnaci_sila(sigma, a, E_ef)` | 6.45 | G = πσ²a/E* | G |
| `hnaci_sila_poddajnost(F, B, dC_da, konstantni_posuv=False)` | 6.48, 6.49 | G = ±F²/(2B)·dC/da | G |
| `G_I(KI, E, mu=None)` | 6.50 | KI²/E, resp. (1 − µ²)KI²/E | GI |
| `G_II(KII, E, mu)` | 6.51 | (1 − µ²)KII²/E | GII |
| `G_III(KIII, E, mu)` | 6.51 | (1 + µ)KIII²/E | GIII |
| `sih_Sc(KIC, mu, Gs)` | 6.53 | KIC²(1 − 2µ)/(4πGs) | Sc |
| `sih_Smin_II(tau, a, mu, Gs)` | 6.54 | τ²a/(12Gs)·(2(1 − µ) − µ²) | Smin |

Celková hnací síla je podle (6.52) součet `G_I + G_II + G_III`.

U empirických převodů KV → KIC (6.33, 6.34) skripta neuvádějí jednotky. Obvykle se dosazuje E v MPa a KV v J a výsledek se bere v MPa·√m; pokud je to pro zadání důležité, ověřte jednotky s vyučujícím. U `platnost_KIC` musí mít KIC/σY rozměr √délky, tedy pro KIC v MPa·√m se B, a a W zadávají v metrech.

### 9.5 Elasticko-plastická lomová mechanika

| Funkce | Rov. | Popis | Vrací |
|---|---|---|---|
| `ctod(V, W, a, r)` | 6.58 | δ = r(W − a)V/(r(W − a) + a) | δ |
| `ctod_elasticko_plasticky(KI, Ye, E_ef, Vp, W, a, rp)` | 6.59 | δe + δpl | δ |
| `J_z_K(KI, E_ef)` | 6.62 | J = KI²/E* | J |
| `J_z_ctod(sigma_Y, delta)` | 6.63 | J = σY·δ | J |
| `ramberg_osgood_hrr(sigma, sigma_0, E, alpha, n)` | 6.64 | ε = ε0·(σ/σ0 + α(σ/σ0)ⁿ) | ε |

---

## 10. Kapitola 7 – Svarové spoje (`unava.svary`)

| Funkce | Rov. | Popis | Vrací |
|---|---|---|---|
| `fat_C(FAT, m=3)` | 7.1 | C = FAT^m·2·10⁶ | C |
| `fat_cykly(delta_sigma, FAT, m=3, mez_cyklu=1e7, za_mezi="vodorovne", m2=22)` | 7.1 | N = C/Δσ^m | N |
| `iiw_omezeni(delta_sigma_n, sigma_n_max, sigma_Y)` | 7.2, 7.3 | Δσn ≤ 1,5σY a σn,max ≤ σY | `bool` |
| `delta_K_svar(a, d_mem, Ym, Mkm, d_b, Yb, Mkb)` | 7.4 | √(πa)·(Δσmem·Ym·Mk,m + Δσb·Yb·Mk,b) | ΔK |
| `paris_rychlost(C, delta_K, m)` | 7.5 | da/dN = C·ΔK^m | da/dN |
| `paris_cykly(C, m, Y, delta_sigma, a0, af)` | 7.6–7.8 | integrace Paris-Erdogana (m = 2 zvlášť) | N |
| `korekce_tloustky(t_eff, n, t_ref=25)` | 7.9 | f(t) = (tref/teff)^n | f(t) |
| `km_efektivni(km_vypoctene, km_zahrnute)` | 7.10 | km,eff = km,calc/km,covered | km,eff |
| `km_kombinovane(km_osove, km_uhlove)` | 7.11 | 1 + (km,axial − 1) + (km,angular − 1) | km |

`fat_cykly` za mezí `mez_cyklu`:
- `za_mezi="vodorovne"`: pod mezí únavy (CAFL) vrací `inf`.
- `za_mezi="sklon"`: křivka pokračuje se sklonem `m2` (výchozí 22).

Pro smyková napětí zadejte `m=5` a `mez_cyklu=1e8`. Efektivní vrubové napětí se posuzuje vůči FAT 225 (ocel) nebo FAT 71 (hliník).

---

## 11. Opravy oproti skriptům

Při implementaci se ukázalo několik překlepů ve skriptech. Knihovna používá opravené tvary:

| Rov. | Ve skriptech | V knihovně | Důvod |
|---|---|---|---|
| 3.36 | ½·√(Σ(σi − σj)²) | **√(½·Σ(σi − σj)²)** | pro jednoosý tah musí vyjít σa = σ1; shodné s (5.21) |
| 4.12 | C′ = C·σc^(−m) | **C′ = C·σc^m** | jinak není S-N křivka v bodě σc spojitá |
| 4.14 | třetí člen (σa3/σa2)^d | **(σa3/σa1)^d** | odpovídá rovnici 4.13 |
| 4.24 | σf·2(Nf)^b | **σ′f·(2Nf)^b** | shodný tvar s Basquinem (3.1) |
| 6.20 | sin³(θ/2) | **sin²(θ/2)** | správný tvar pro rovinnou napjatost |
| 6.26 | τY = √3·σY | **τY = σY/√3** | podle HMH |

Mez kluzu ve smyku podle opraveného vztahu používá i funkce `rp_III`.

---

## 12. Řešené příklady

### Proměnné zatížení – stékající déšť + Morrow + Palmgren-Miner

```python
import unava

napeti = [170, -270, 390, 130, 260, -160, 130, -420, 180, 60, 320, -90]
B, tab = unava.zivotnost_bloku(napeti, sigma_f=1020, b=-0.08432,
                               nazvy=list("ABCDEFGHIJKL"))
print(f"B = {B:.0f} bloků")
print(f"N = {B * len(tab):.0f} cyklů")
```

```
B = 32787 bloků
N = 196723 cyklů
```

### Parametry E-N křivky z naměřených dat

```python
import unava

E = 206000
eps_a = [0.02433, 0.01752, 0.01349, 0.00989, 0.00699, 0.00488, 0.00395, 0.00293, 0.00265]
sig_a = [626.9, 563.6, 500.5, 465.5, 423.9, 380.4, 361.8, 335.5, 322.5]
Nf    = [24, 56, 133, 594, 897, 5700, 8630, 124450, 536512]

sigma_f, b = unava.basquin_regrese(sig_a, Nf)
eps_pl = [e - s / E for e, s in zip(eps_a, sig_a)]
eps_f, c = unava.basquin_regrese(eps_pl, Nf)

print(sigma_f, b, eps_f, c)
print("n' =", unava.n_z_mcb(b, c), " K' =", unava.K_z_mcb(sigma_f, eps_f, b, c))
print("2NT =", unava.tranzitni_pulcykly(sigma_f, b, eps_f, c, E))
```

`basquin_regrese` funguje pro jakoukoli mocninnou závislost na 2Nf, takže stejnou funkcí získáte i ε′f a c.

### Vrub v nízkocyklové oblasti – Neuber + SWT

```python
import unava

E, K, n = 206000, 1370, 0.216
sigma_f, b, eps_f, c = 742.3, -0.0658, 0.0603, -0.3085

sigma, eps = unava.neuber(2.5, 200, E, K, n)
Nf = unava.swt_cykly(sigma, eps, sigma_f, b, eps_f, c, E)
print(f"σ = {sigma:.0f} MPa, ε = {eps:.5f}, Nf = {Nf:.0f}")
```

### Šíření trhliny ve svaru – Paris-Erdogan

```python
import unava

N = unava.paris_cykly(C=1e-12, m=3, Y=1.12, delta_sigma=100, a0=0.15e-3, af=10e-3)
print(f"N = {N:.3e}")
```

---

## 13. Testy

Repozitář obsahuje sadu testů ve složce `tests/`. Kontrolují mimo jiné:
- výsledek úlohy ze skript (B = 32 787),
- vzorový příklad stékajícího deště z normy ASTM E1049,
- inverze všech implicitních rovnic,
- spojitost Haibachovy křivky a Paris-Erdoganova vztahu v okolí m = 2,
- kompatibilitu parametrů MCB a Ramberg-Osgooda,
- chování na hraničních vstupech.

Spuštění:

```
pip install pytest
pytest
```

Implementace stékajícího deště byla navíc porovnána s knihovnou [`rainflow`](https://pypi.org/project/rainflow/) na 10 000 náhodných historiích, a to v obou režimech (jeden průchod i opakovaný blok). Výsledky se shodovaly ve všech případech.

