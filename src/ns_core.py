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

class Core(object):
    """Solve model equations for the core of a neutron star

    Attributes
    ==========

    m : positive float, nucleon mass (MeV)

    sf : pandas DataFrame, state function of the star core. Columns
    are barionic density, energy density and pressure. Defined after
    call of StateFunction method

    Methods
    =======

    define_density : Define a numpy array w/ the barionic density of
    the star core

    """

    def __init__(self, m, Yp=0.5, verbose=False):
        """Solve model equations for the core of a neutron star

        Parameters
        ==========

        m : positive float, nucleon mass (MeV)

        Yp : float (0 <= Yp <= 1, default 0.5), proton fraction

        verbose : boolean (default False), verbose output

        """
        super(Core, self).__init__()

        self.vprint = print if verbose else lambda *a, **k: None

        self.__m = m
        self.__Yp = Yp

        self.__hc = 197.33      # MeV.fm (?)

        self.__gw = 12.8679990675   # coupling constant omega-nucleon
        self.__gs = 10.2170005383   # coupling constant sigma-nucleon
        self.__gRho = 8.92188320928  # coupling constant rho-nucleon
        self.__ms = 508.194          # sigma mass in MeV
        self.__mw = 782.501          # omega mass in MeV
        self.__mRho = 763.0          # rho mass in MeV

        return None

    @property
    def m(self):
        return self.__m

    @m.setter
    def m(self, m):
        if not isinstance(m, (int, float)):
            raise TypeError(f"m must be float, {type(m)} passed")

        if m < 0:
            warnings.warn(f"Mass must be positive-valued, {m} was passed "
                          f"setting new value to zero")
            self.__m = 0

    @property
    def Yp(self):
        return self.__Yp

    @Yp.setter
    def Yp(self, Yp):
        if not isinstance(Yp, (int, float)):
            raise TypeError(f"Yp must be float, {type(Yp)} passed")
        if Yp < 0 or Yp > 1:
            warnings.warn(f"Proton fraction (Yp) must be between 0 and 1, "
                          f"Yp = {Yp} instead, changed to default 0.5")
            self.__Yp = 0.5

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
