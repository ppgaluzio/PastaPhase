#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as pl
import numpy as np

import pathmagic                # noqa

import sigma as ss
import energy_and_pressure as ep

hc = 197.33
rho_inf = 0.05
rho_sup = 1.00
Npts = 200

rho = np.linspace(rho_inf*hc**3,
                  rho_sup*hc**3,
                  Npts)

m = 939

g0 = 10.608                    #coupling constant omega-nucleon
gs = 9.5684
ms = 550
mw = 783

e = np.zeros_like(rho)
p = np.zeros_like(rho)

for i, r in enumerate(rho):
    k = (1.5 * (np.pi**2) * r)**(1/3)
    w0 = g0 * r / mw**2
    sigma = ss.SolveSigma(gs, ms, m, k, tol=1.0e-5, n_seeds=100)
    print(i, 'k = ', k, ' sigma = ', sigma)
    e[i] = ep.energy(ms, sigma, mw, w0, k, gs, m) / hc**3
    p[i] = ep.pressure(ms, sigma, mw, w0, k, gs, m) / hc**3

pl.plot(e, p)
pl.grid()
pl.ylabel(r"$p$")
pl.xlabel(r"$\varepsilon$")
pl.show()
