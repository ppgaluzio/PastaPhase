# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
import warnings

from .sigma import SolveSigma
from .energy_and_pressure import energy
from .energy_and_pressure import pressure


# XXX: we need to make sure that the units are consistent throughout
# the code, and eventualy pass in the class which unit system the user
# will prefer

# XXX: We need to define what are the parameters that are common to
# the star as a whole and maybe include these into a star class as a
# parent class

class Core():
    """
    Core
    ----

    Solve model equations for Core of a neutron star

    input
    -----

    m :: mass (unit?)
    hc :: unity converting factor

    """
    def __init__(self, m, hc):
        super(Core, self).__init__()
        self.m = m
        self.hc = hc

    @property
    def m(self):
        return self.__m

    @m.setter
    def m(self, m):
        if m < 0:
            warnings.warn(f"Mass must be positive-valued, {m} was passed "
                          "setting new value to zero")
            self.__m = 0

    def define_density(self, rho_inf, rho_sup, npts=100):
        """
        define_density
        --------------

        Define a numpy array w/ the density of the star core

        input
        -----
        rho_inf :: inferior bound to the array
        rho_sup :: superior bound to the array
        npts (100) :: number of points to be calculated
        """

        if not isinstance(rho_inf, (int, float)):
            raise TypeError(f"rho_inf of type{type(rho_inf)}, "
                            "should be int or float}")

        if not isinstance(rho_sup, (int, float)):
            raise TypeError(f"rho_inf of type{type(rho_inf)}, "
                            "should be int or float}")

        self.rho = np.linspace(rho_inf * self.hc**3,
                               rho_sup * self.hc**3,
                               npts)
        self.npts = npts

        return None
