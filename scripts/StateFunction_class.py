#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pathmagic                # noqa

from ns_core import Core

ns = Core(m=939.0, verbose=True)
ns.StateFunction(nB_inf=0.05, nB_sup=0.70, npts=200)
ns.PlotaStateFunction()
