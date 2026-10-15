# unava

Vzorce z předmětu Únava materiálu podle skript *Životnost a spolehlivost* (Kořínek, Tošková, Fusek; VŠB-TUO).
Čísla rovnic jsou v docstringu každé funkce. Kompletní popis všech funkcí je v [DOKUMENTACE.md](DOKUMENTACE.md).

## Instalace

```
pip install unava
```

Případně přímo z GitHubu: `pip install git+https://github.com/Fradek/Unava.git`

Pro vývoj (změny se projeví bez přeinstalování):

```
pip install -e ".[test]"
```

## Moduly

| modul | kapitola |
|---|---|
| `unava.zaklady` | 2 – základní pojmy |
| `unava.sn` | 3 – vysokocyklová únava (S-N) |
| `unava.kumulace` | 4 – kumulace poškození, stékající déšť |
| `unava.en` | 5 – nízkocyklová únava (E-N) |
| `unava.lomova_mechanika` | 6 – lomová mechanika |
| `unava.svary` | 7 – svarové spoje |

Všechny funkce jdou volat i přímo z `unava`.

## Příklad

```python
import unava

B, _ = unava.zivotnost_bloku(
    [150, -230, 360, 110, 240, -140, 100, -380, 160, 50, 290, -70],
    sigma_f=900, b=-0.09)
print(f"B = {B:.0f}")   # B = 10658
```

## Testy

```
pytest
```
