from sympy import *
# from sympy.tensor import Idx, IndexedBase
from sympy.printing import fcode

z = { i-4:Symbol('z__array__%d' % i) for i in range(9) }
P = { i-4:Symbol('P__array__%d' % i) for i in range(9) }
T = { i-4:Symbol('T__array__%d' % i) for i in range(9) }
F = { i-4:Symbol('F__array__%d' % i) for i in range(9) }

alpha = Symbol("alpha")
sigma = Symbol('cgs_stef')
c = Symbol('cgs_c')
Omega = Symbol("omega")
kappa = Symbol("kram_es")
temp_eff = Symbol("temp_eff")
flux_acc = Symbol("flux_acc")

init_printing()

k = Symbol('cgs_boltz')
mH = Symbol('cgs_mhydr')
miu = 0.5
# K = k / (mH * miu)
K = Symbol('K')

# kappa = Function('\\kappa')(rho,T)

def RP(z,P,T,F):
    rho = P / (K*T)
    return - rho * ( Omega**2 * z - kappa * F / c )
def RT(z,P,T,F):
    rho = P / (K*T)
    return - 3 * kappa * rho / ( 16 * sigma * T**3 ) * F
def RF(z,P,T,F):
    return alpha * Omega * ( P + 4 * sigma / (3*c) * T**4 )


def fixi(e,ind=None):
    if type(ind) == int:
        ix = Number(ind)
    ix = Idx('i')
    return e.subs((z[0]+z[1])/2,Symbol('zmid')) \
            .subs(z[1]-z[0],Symbol('dz')) \
            .subs(P[1]-P[0],Symbol('dP')) \
            .subs(T[1]-T[0],Symbol('dT')) \
            .subs(F[1]-F[0],Symbol('dF')) \
            .subs((P[0]+P[1])/2,Symbol('Pmid')) \
            .subs((T[0]+T[1])/2,Symbol('Tmid')) \
            .subs((F[0]+F[1])/2,Symbol('Fmid')) \
            .subs(z[-1],Indexed('z',ix-1)) \
            .subs(z[0],Indexed('z',ix+0)) \
            .subs(z[1],Indexed('z',ix+1)) \
            .subs(P[-1],Indexed('P',ix-1)) \
            .subs(P[0],Indexed('P',ix+0)) \
            .subs(P[1],Indexed('P',ix+1)) \
            .subs(T[-1],Indexed('T',ix-1)) \
            .subs(T[0],Indexed('T',ix+0)) \
            .subs(T[1],Indexed('T',ix+1)) \
            .subs(F[-1],Indexed('F',ix-1)) \
            .subs(F[0],Indexed('F',ix+0)) \
            .subs(F[1],Indexed('F',ix+1))


AP = (P[1] - P[0] ) / (z[1] - z[0]) \
    - RP((z[0]+z[1])/2,(P[0]+P[1])/2,(T[0]+T[1])/2,(F[0]+F[1])/2)

AT = (T[1] - T[0]) / (z[1] - z[0]) \
    - RT((z[0]+z[1])/2,(P[0]+P[1])/2,(T[0]+T[1])/2,(F[0]+F[1])/2)

AF = (F[1] - F[0] ) / (z[1] - z[0])  \
    - RF((z[0]+z[1])/2,(P[0]+P[1])/2,(T[0]+T[1])/2,(F[0]+F[1])/2)

E = [ P[0],P[1],T[0],T[1],F[0],F[1] ]

B1 = F[0]
B2 = F[1] - flux_acc
B3 = F[1] - 2 * sigma * T[1]**4

f = open('out.py','w')

f.write("""#!/usr/bin/env python

import numpy as np
from cgs import *

def compute_corrections(z,P,T,F,N,par):

    A = np.zeros(3*N)
    M = np.zeros([3*N,3*N])

    AP = A[0:(N-1)]
    AT = A[(N-1):2*(N-1)]
    AF = A[2*(N-1):3*(N-1)]

    MPP = M[0:(N-1),0:N]
    MPT = M[0:(N-1),N:2*N]
    MPF = M[0:(N-1),2*N:3*N]
    MTP = M[(N-1):2*(N-1),0:N]
    MTT = M[(N-1):2*(N-1),N:2*N]
    MTF = M[(N-1):2*(N-1),2*N:3*N]
    MFP = M[2*(N-1):3*(N-1),0:N]
    MFT = M[2*(N-1):3*(N-1),N:2*N]
    MFF = M[2*(N-1):3*(N-1),2*N:3*N]

    AB1 = A[-1]
    AB2 = A[-2]
    AB3 = A[-3]

    MB1P = M[-1,0:N]
    MB1T = M[-1,N:2*N]
    MB1F = M[-1,2*N:3*N]
    MB2P = M[-2,0:N]
    MB2T = M[-2,N:2*N]
    MB2F = M[-2,2*N:3*N]
    MB3P = M[-3,0:N]
    MB3T = M[-3,N:2*N]
    MB3F = M[-3,2*N:3*N]

    K = 2 * cgs_boltz / cgs_mhydr

    temp_eff = par['temp_eff']
    flux_acc = par['flux_acc']
    alpha = par['alpha']
    omega = par['omega']
    kram_es = par['kram_es']

    for i in range(N-1):
        zmid = (z[i] + z[i+1])/2
        Pmid = (P[i] + P[i+1])/2
        Tmid = (T[i] + T[i+1])/2
        Fmid = (F[i] + F[i+1])/2
        dz = z[i+1] - z[i]
        dP = P[i+1] - P[i]
        dT = T[i+1] - T[i]
        dF = F[i+1] - F[i]
        AP [i]     = %s
        MPP[i,i]   = %s
        MPP[i,i+1] = %s
        MPT[i,i]   = %s
        MPT[i,i+1] = %s
        MPF[i,i]   = %s
        MPF[i,i+1] = %s
        AT [i]     = %s
        MTP[i,i]   = %s
        MTP[i,i+1] = %s
        MTT[i,i]   = %s
        MTT[i,i+1] = %s
        MTF[i,i]   = %s
        MTF[i,i+1] = %s
        AF [i]     = %s
        MFP[i,i]   = %s
        MFP[i,i+1] = %s
        MFT[i,i]   = %s
        MFT[i,i+1] = %s
        MFF[i,i]   = %s
        MFF[i,i+1] = %s
        if i == 0:
            AB1       = %s
            MB1F[i]   = %s
        if i == N-2:
            AB2       = %s
            MB2F[i+1] = %s
            AB3       = %s
            MB3F[i+1] = %s
            MB3T[i+1] = %s

    #np.save('M',M)
    #np.save('A',A)

    delt = np.linalg.solve(M,-A)
    return (delt[0:N],delt[N:2*N],delt[2*N:3*N])

""" % tuple( fixi(e) for e in \
      [ AP ] + [ AP.diff(v) for v in E ] \
    + [ AT ] + [ AT.diff(v) for v in E ] \
    + [ AF ] + [ AF.diff(v) for v in E ] \
    + [ B1, B1.diff(F[0]) ] \
    + [ B2, B2.diff(F[1]) ] \
    + [ B3, B3.diff(F[1]), B3.diff(T[1]) ] \
    ))
    #   [ LBC(AP) ] + [ LBC(AP).diff(v) for v in E ] \
    # + [ LBC(AT) ] + [ LBC(AT).diff(v) for v in E ] \
    # + [ LBC(AF) ] + [ LBC(AF).diff(v) for v in E ] \
    # + [ RBC(AP) ] + [ RBC(AP).diff(v) for v in E ] \
    # + [ RBC(AT) ] + [ RBC(AT).diff(v) for v in E ] \
    # + [ RBC(AF) ] + [ RBC(AF).diff(v) for v in E ] \


f.close()
