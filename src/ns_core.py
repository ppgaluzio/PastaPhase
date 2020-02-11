# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
import warnings

from .sigma import SolveSigma
from .energy_and_pressure import energy
from .energy_and_pressure import pressure


class Core():
    """Documentation for Core

    """
    def __init__(self, m):
        super(Core, self).__init__()
        self.m = m

    @property
    def m(self):
        return self.__m

    @m.setter
    def m(self, m):
        if m < 0:
            warnings.warn(f"Mass must be positive-valued, {m} was passed "
                          "setting new value to zero")
            self.__m = 0
