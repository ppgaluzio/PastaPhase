# -*- coding: utf-8 -*-

import numpy as np
from scipy.integrate import quad


def _int_e(k, m, gs, sigma):
    return np.sqrt(k**2 + (m - gs * sigma)**2) * k**2


def _int_p(k, m, gs, sigma):
    return k**4 / np.sqrt(k**2 + (m - gs * sigma)**2)


def pressure(ms, sigma, mw, w0, k, gs, m):
    """
    pressure
    --------

    calculate pressurey, eq. 4.160 from compact stars book

    parameters:

    ms: m_sigma
    sigma: scalar field
    mw: m_omega
    w0: omega_0
    k: fermi energy
    gs: g_sigma
    m: mass
    """

    p = - 0.5 * ms**2 * sigma**2 \
        + 0.5 * mw**2 * w0**2 \
        + (1/3) * (2/np.pi**2) \
        * quad(_int_p, 0, k, args=(m, gs, sigma))[0]

    return p


def energy(ms, sigma, mw, w0, k, gs, m):
    """
    energy
    ------

    calculate energy density, eq. 4.160 from compact stars book

    parameters:

    ms: m_sigma
    sigma: scalar field
    mw: m_omega
    w0: omega_0
    k: fermi energy
    gs: g_sigma
    m: mass
    """

    e = + 0.5 * ms**2 * sigma**2 \
        + 0.5 * mw**2 * w0**2 \
        + (2 / np.pi**2) \
        * quad(_int_e, 0, k, args=(m, gs, sigma))[0]

    return e
