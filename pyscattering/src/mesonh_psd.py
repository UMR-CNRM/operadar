#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 12:07:14 2023

@author: augros
"""

import math
from math import gamma

def PSD(M, D, CCLOUD, P3, typeh, a, b, nuconst, alpha, c, x, mumax, Nmoments):
    if Nmoments == 2:
        No = P3
        if CCLOUD == 'LIMT' and typeh == 'rr':
            Dm = (M / (a * No)) ** (1 / b)
            mucalc = 38.66 * math.exp(-1.59 * (Dm * 10 ** 3))
            mu = max(min(mucalc, mumax), 0.0)
            nu = mu + 1
        else:
            nu = nuconst
        lamb = (a * No * gamma(nu + b / alpha) / (M * gamma(nu))) ** (1 / b)
    else:
        nu = nuconst
        lamb = (M * gamma(nu) / (a * c * gamma(nu + b / alpha))) ** (1 / (x - b))
        No = c * (lamb ** x)

    lam_exp = min(10.0**35, lamb**(alpha*nu))
    

    N = math.exp(-(lamb*D)**alpha)/math.gamma(nu)*(D**(alpha*nu-1)) * No*alpha*lam_exp



    return lamb, N