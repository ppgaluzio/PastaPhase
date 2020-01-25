#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as pl
import numpy as np

import pathmagic                # noqa

import sigma as ss
import energy_and_pressure as ep

kk = np.linspace(0, 100, 50)
ms = 100
mw = 25
m = 2
gs = 150
w0 = 0.2

e = np.zeros_like(kk)
p = np.zeros_like(kk)

for i, k in enumerate(kk):
    sigma = ss.SolveSigma(gs, ms, m, k, tol=1.0e-5, n_seeds=100)
    print(i, 'k = ', k, ' sigma = ', sigma)
    e[i] = ep.energy(ms, sigma, mw, w0, k, gs, m)
    p[i] = ep.pressure(ms, sigma, mw, w0, k, gs, m)

pl.plot(p, e)
pl.grid()
pl.xlabel(r"$p$")
pl.ylabel(r"$\varepsilon$")
pl.show()
