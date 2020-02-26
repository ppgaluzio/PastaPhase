# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
import warnings

from sigma import SolveSigma
from energy_and_pressure import energy
from energy_and_pressure import pressure

# XXX: we need to make sure that the units are consistent throughout
# the code, and eventualy pass in the class which unit system the user
# will prefer

# XXX: We need to define what are the parameters that are common to
# the star as a whole and maybe include these into a star class as a
# parent class


class models(object):
    def __init__(self):
        self.update_f = None
        return

    def coupling_constants(self, gw, gs, grho):
        """define the coupling constants of the model

        Parameters
        ==========

        gw : float, coupling constant omega-nucleon
        gs : float, coupling constant sigma-nucleon
        gRho : float, coupling constant rho-nucleon

        """
        self.gw = gw
        self.gs = gs
        self.grho = grho

    def masses(self, mw, ms, mrho):
        """define the masses of the model

        Parameters
        ==========

        mw : float, omega mass in MeV
        ms : float, sigma mass in MeV
        mrho : float, rho mass in MeV

        """
        self.mw = mw
        self.ms = ms
        self.mrho = mrho

    def def_update(self, update_f):
        """define update of the model

        Parameter
        =========

        update_f : function
                   gw, gs, gRho, mw, ms, mRho = update_f(density)
        """
        self.update_f = update_f

    def __update(self, density):
        if self.update_f is not None:
            self.gw, self.gs, self.grho, self.mw, self.ms, self.mrho = \
                self.update_f(density)


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

    __MODELS = {}

    __MODELS['nl3'] = models()
    __MODELS['nl3'].coupling_constants(12.8679990675,
                                       10.2170005383,
                                       8.92188320928)
    __MODELS['nl3'].masses(782.501,
                           508.194,
                           763.0)

    def __init__(self, m, model, Yp=0.5, verbose=False):
        """Solve model equations for the core of a neutron star

        Parameters
        ==========

        m : positive float, nucleon mass (MeV)

        model: string, name of the model to be implemented. Valid values: 'nl3'

        Yp : float (0 <= Yp <= 1, default 0.5), proton fraction

        verbose : boolean (default False), verbose output

        """
        super(Core, self).__init__()

        self.vprint = print if verbose else lambda *a, **k: None

        self.__m = m
        self.__Yp = Yp

        if isinstance(model, str):
            if model in self.__MODELS.keys():
                self.model = self.__MODELS[model]
            else:
                raise ValueError(f"model must be one of {self.__MODELS.keys()}"
                                 f", {model} given")
        elif isinstance(model, models):
            self.model = model
        else:
            raise TypeError("model must be valid string or instance of models")

        self.__hc = 197.33      # MeV.fm (?)

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

    def __define_density(self, nB_inf, nB_sup, npts):
        """Define a numpy array w/ the barionic density of the star core

        Parameters
        ==========

        nB_inf : positive float, minimum value for the barionic
        density (fm^-3)

        nB_sup : positive float, maximum value for the barionic
        density (fm^-3). nB_sup must be larger than nB_inf.

        npts : int (default 100), number of density values calculated
        between nB_inf and nB_sup

        Raises
        ======

        ValueError : if nBsup <= nb_inf, or if both values are <= 0
        """

        if not isinstance(nB_inf, (int, float)):
            raise TypeError(f"nB_inf of type{type(nB_inf)}, "
                            f"should be int or float")

        if not isinstance(nB_sup, (int, float)):
            raise TypeError(f"rho_inf of type{type(nB_inf)}, "
                            f"should be int or float")

        if nB_inf <= 0:
            raise ValueError(f"nB_inf must be positive, "
                             f"{nB_inf} passed instead.")

        if nB_sup <= 0:
            raise ValueError(f"nB_sup must be positive, "
                             f"{nB_sup} passed instead.")

        if nB_sup <= nB_inf:
            raise ValueError(f"nB_sup must be > nB_inf, "
                             f"instead nB_sup = {nB_sup} "
                             f"and nB_inf = {nB_inf}")

        self.__nB = np.linspace(nB_inf * self.__hc**3,
                                nB_sup * self.__hc**3,
                                npts)

        return None

    def StateFunction(self, nB_inf, nB_sup, npts=100):
        """Calculate the state function of the star core

        Parameters
        ==========

        nB_inf : positive float, minimum value for the barionic
        density (fm^-3)

        nB_sup : positive float, maximum value for the barionic
        density (fm^-3). nB_sup must be larger than nB_inf.

        npts : int (default 100), number of density values calculated
        between nB_inf and nB_sup

        Attribute
        =========

        sf : pandas DataFrame, state function of the star core. Columns
        are barionic density, energy density and pressure

        """

        self.__define_density(nB_inf, nB_sup, npts=npts)

        e = np.zeros_like(self.__nB)   # baryons energy density (MeV/fm^3)
        p = np.zeros_like(self.__nB)   # barionic pressure  (MeV/fm^3)

        for i, n in enumerate(self.__nB):

            self.model.__update(n)      # update model

            k = (1.5 * (np.pi**2) * n)**(1/3)        # fermi momentum
            w0 = self.model.gw * n / self.model.mw**2  # meson omega
            nB_p = self.__Yp * n         # density of protons
            nB_n = (1.0 - self.__Yp) * n  # density of electrons

            # meson rho
            rho = 0.5 * self.model.grho * (nB_p - nB_n) / self.model.mw**2

            sigma = SolveSigma(self.model.gs, self.model.ms, self.m, k,
                               tol=1.0e-5, n_seeds=100)
            self.vprint(i, 'k = ', k, ' sigma = ', sigma)
            e[i] = energy(self.model.ms, self.model.mrho, rho, sigma,
                          self.model.mw, w0, k, self.model.gs, self.m)
            p[i] = pressure(self.model.ms, self.model.mrho, rho, sigma,
                            self.model.mw, w0, k, self.model.gs, self.m)

        e = e / self.__hc**3
        p = p / self.__hc**3

        self.sf = pd.DataFrame(np.array([self.__nB, e, p]).T,
                               columns=['nB', 'e', 'p'])

        return None

    def PlotaStateFunction(self, show=True, filename=None):
        """Plota the state function

        Parameters
        ==========

        show : boolean (default True), if true show the graphs at runtime

        filename : string (default None), name of the file to save the
        plot, if None only show the figure

        """

        if not hasattr(self, 'sf'):
            warnings.warn(f"StateFunction method must be called prior",
                          RuntimeWarning)
            return None

        fig, ax = pl.subplots()

        ax.plot(self.sf['e'], self.sf['p'])
        ax.grid()
        ax.set_ylabel(r'$p$')
        ax.set_xlabel(r'$\varepsilon$')

        if show:
            fig.show()

        if filename is not None:
            fig.savefig(filename)

        if show:
            pl.show()

        return None
