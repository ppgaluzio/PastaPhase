#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pathmagic                # noqa

from ns_core import Core
from ns_core import models


def update(density):
    gw = 12.8679990675
    gs = 10.2170005383
    grho = 8.92188320928
    mw = 782.501 * density
    ms = 508.194
    mrho = 763.0 * density
    return gw, gs, grho, mw, ms, mrho


model = models()
model.coupling_constants(12.8679990675,
                         10.2170005383,
                         8.92188320928)
model.masses(782.501,
             508.194,
             763.0)
model.def_update(update)

# crazy user model

ns = Core(m=939.0, model=model, verbose=True)
ns.StateFunction(nB_inf=0.05, nB_sup=0.70, npts=200)
ns.PlotaStateFunction()

# Standard nl3

nl3 = Core(m=939.0, model='nl3', verbose=True)
nl3.StateFunction(nB_inf=0.05, nB_sup=0.70, npts=200)
nl3.PlotaStateFunction()
