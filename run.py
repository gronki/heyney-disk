#!/usr/bin/env python

import numpy as np
from cgs import *
from pyminiconf import *
import matplotlib.pyplot as plt
from time import sleep
from out import compute_corrections
from col2python import col2python
from sys import argv

K = 2 * cgs_boltz / cgs_mhydr

dat,par = col2python(argv[1])

N = int(argv[2])

ztop = 16
z = np.linspace(0,par['zscale'] * ztop,N)

g = np.exp( -z**2 / (2 * par['zscale']**2) )
T = np.linspace(par['temp_0'], 0.841*par['temp_eff'],N)
F = np.linspace(0,par['flux_acc'],N)
P =  K * par['temp_0'] * par['rho_0'] * np.linspace(1,1e-12,N)


nslope = 8
niter = 16
for it in range(niter):

    # print ("Iteracja %d" % it)

    dP,dT,dF = compute_corrections(z,P,T,F,N,par)

    t = (1.0 + it) / nslope if it < nslope else 1.0

    fig = plt.figure(figsize=(18,8))
    fig.suptitle("Iteracja %d" % it)

    ax = plt.subplot(1,3,1)
    ax.set_xlim(0,ztop)
    ax.set_title("P gas")
    ax.plot(dat['h'],dat['pgas'], color='#8E8E8E')
    #ax.set_yscale('log')
    #ax.set_ylim(0,K * par['temp_0'] * par['rho_0'] * 1.5)
    ax.plot(z/par['zscale'],P+dP,'--', color='#409350')
    ax.plot(z/par['zscale'],P,'-',color='black', linewidth=2)
    P = P + dP * t
    P = np.where(P < 0, 0, P)

    ax = plt.subplot(1,3,2)
    ax.set_xlim(0,ztop)
    #ax.set_yscale('log')
    #ax.set_ylim(0,par['temp_0']*1.5)
    ax.set_title("Temperatura")
    ax.axhline(par['temp_eff']*0.841, linestyle='--', linewidth=0.5, color='#C2AA56')
    ax.plot(dat['h'],dat['tgas'], color='#8E8E8E')
    ax.plot(z/par['zscale'], T+dT, '--', color='#409350')
    ax.plot(z/par['zscale'],T,'-',color='black', linewidth=2)
    T = T + dT * t
    T = np.where(T < par['temp_eff']*0.2, par['temp_eff']*0.2, T)

    ax = plt.subplot(1,3,3)
    ax.set_xlim(0,ztop)
    ax.set_title("Flux")
    #ax.set_yscale('log')
    ax.axhline(par['flux_acc'], linestyle='--', linewidth=0.5, color='#C2AA56')
    ax.plot(dat['h'],dat['frad'], color='#8E8E8E')
    ax.set_ylim(0,par['flux_acc'] * 1.2)
    ax.plot(z/par['zscale'],F+dF,'--', color='#409350')
    ax.plot(z/par['zscale'],F,'-',color='black', linewidth=2)
    F = F + dF * t
    F = np.where(F < 0, 0, F)
    plt.savefig('plt/f-%d-%03d.png' % (np.log2(N),it))
    plt.close(fig)
