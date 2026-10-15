import math

_DOLNI_2N = 1e-3
_HORNI_2N = 1e30


def _koren(f, dolni, horni, iteraci=200):
    fd = f(dolni)
    for _ in range(iteraci):
        stred = math.sqrt(dolni * horni)
        fs = f(stred)
        if fd * fs <= 0:
            horni = stred
        else:
            dolni, fd = stred, fs
    return math.sqrt(dolni * horni)


def _pulcykly(f, cil):
    if cil <= 0:
        raise ValueError("cílová hodnota musí být kladná")
    if cil <= f(_HORNI_2N):
        return math.inf
    if cil > f(_DOLNI_2N):
        raise ValueError("zatížení je mimo rozsah křivky životnosti (2Nf < 1e-3)")
    return _koren(lambda x: math.log(f(x)) - math.log(cil), _DOLNI_2N, _HORNI_2N)
