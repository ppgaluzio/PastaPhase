#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as pl
import numpy as np

import pathmagic                # noqa

import sigma as ss
import energy_and_pressure as ep

hc = 197.33                    # MeV.fm
nB_inf = 0.05                 # minimum value for the barionic density in fm^-3 
nB_sup = 0.70                 # maximum value for the barionic density in fm^-3 
Npts = 200                     # number of points to the uniform matter EoS table

#NLE set of parameters
m = 939.0                       #nucleon mass in MeV
gw = 12.8679990675              #coupling constant omega-nucleon
gs = 10.2170005383              #coupling constant sigma-nucleon
gRho = 8.92188320928            #coupling constant rho-nucleon
ms = 508.194                    # sigma mass in MeV    
mw = 782.501                    # omega mass in MeV   
mRho = 763.0                    # rho mass in MeV

#defines proton fraction: Yp = (number of protons)/(number of barions)
Yp = 0.3
#nota pro Paulo, mudei a nomenclatura, rho > nB, eh o mais costumeiro. e assim posso chamar o meson rho, de rho.
nB = np.linspace(nB_inf*hc**3,
                  nB_sup*hc**3,
                  Npts)

#nB_p = Yp * nB                        #density of protons
#nB_n = (1.0 - Yp ) * nB               #density of electrons


e = np.zeros_like(nB)                  # baryons energy density (MeV/fm^3)
p = np.zeros_like(nB)                  # barionic pressure  (MeV/fm^3)

for i, n in enumerate(nB):
    k = (1.5 * (np.pi**2) * n)**(1/3)    # fermi momentum
    w0 = gw * n / mw**2                  # meson omega
    nB_p = Yp * n                        # density of protons
    nB_n = (1.0 - Yp ) * n               # density of electrons
    rho = 0.5 * gRho * (nB_p - nB_n) / mw**2          # meson rho

    sigma = ss.SolveSigma(gs, ms, m, k, tol=1.0e-5, n_seeds=100)
    print(i, 'k = ', k, ' sigma = ', sigma)
    e[i] = ep.energy(ms, mRho, rho, sigma, mw, w0, k, gs, m) / hc**3
    p[i] = ep.pressure(ms, mRho, rho, sigma, mw, w0, k, gs, m) / hc**3

pl.plot(e, p)
pl.grid()
pl.ylabel(r"$p$")
pl.xlabel(r"$\varepsilon$")
pl.show()
