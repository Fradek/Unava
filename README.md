# unava

Vzorce z předmětu Únava materiálu podle skript *Životnost a spolehlivost* (Kořínek, Tošková, Fusek; VŠB-TUO).
Čísla rovnic jsou v docstringu každé funkce. Kompletní popis všech funkcí je v [DOKUMENTACE.md](DOKUMENTACE.md).

## Instalace

Ve složce s `pyproject.toml`:

```
pip install .
```

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
    [170, -270, 390, 130, 260, -160, 130, -420, 180, 60, 320, -90],
    sigma_f=1020, b=-0.08432)
print(f"B = {B:.0f}")   # B = 32787
```

## Testy

```
pytest
```
