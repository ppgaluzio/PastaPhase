# -*- coding: utf-8 -*-

"""
sigma.py
--------

Solves the self consistent equation 4.151 from the Glendenning book.
The goal is to find the value of the scalar field sigma in the equation
"""

import numpy as np


def _SquaredResidue(s, gs, ms, m):
    """
    Squared residue is an internal module function that returns the squared of
    the residue function.

    Fiven y = f(x), the residue function is such as g(x) = y - f(x) = 0,
    so solving for y = f(x) is the same is minimizing g(x)**2

    Parameters:
    s: sigma
    gs: g_sigma
    ms: m_sigma
    m: m
    """


    return Residue**2
