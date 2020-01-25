# -*- coding: utf-8 -*-

"""
sigma.py
--------

Solves the self consistent equation 4.151 from the Glendenning book.
The goal is to find the value of the scalar field sigma in the equation
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize


class ConvergenceError(Exception):
    """
    Raised when the alogithm did not converge
    """
    pass


def _integrand(k, m, gs, s):
    return k**2 * (m - gs * s) / np.sqrt(k**2 + (m - gs * s)**2)


def _SquaredResidue(s, gs, ms, m, k):
    """
    Squared residue is an internal module function that returns the squared of
    the residue function.

    Fiven y = f(x), the residue function is such as g(x) = y - f(x) = 0,
    so solving for y = f(x) is the same is minimizing g(x)**2

    Parameters:
    s: sigma
    gs: g_sigma
    ms: m_sigma
    m: m (mass)
    k: fermi energy
    """

    y = gs * s
    Int = quad(_integrand, 0, k, args=(m, gs, s))[0]
    Residue = y - (gs/ms)**2 * (2/np.pi**2) * Int

    return Residue**2


def SolveSigma(gs, ms, m, k, n_seeds=10, verbose=False):
    """
    SolveSigma
    ----------

    Solve self consistent equation to find the values of the scalar field sigma

    Parameters:
    gs: g_sigma
    ms: m_sigma
    m: m (mass)
    k: fermi energy
    n_seeds: number of different seeds for the optimization
    """

    s_min = 0
    s_max = 10
    s_amp = s_max - s_min

    FoundSolution = False
    s = None
    min_res = np.inf

    s_seeds = np.random.rand(n_seeds) * s_amp + s_min

    for s0 in s_seeds:
        res = minimize(_SquaredResidue, s0, args=(gs, ms, m, k))

        if res.success:
            # if s is None or res.fun < min_res:
            s = res.x
            min_res = res.fun
            if verbose:
                print("Found solution with residue {}".format(min_res))
            FoundSolution = True
            break

    if FoundSolution:
        return s
    else:
        print("Minimal value found for residue = {}".format(min_res))
        raise ConvergenceError("Minimize did not converge")
