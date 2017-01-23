#!/usr/bin/env python

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
        AP [i]     = dP/dz + Pmid*(-Fmid*kram_es/cgs_c + omega**2*zmid)/(K*Tmid)
        MPP[i,i]   = -1/dz + (-Fmid*kram_es/cgs_c + omega**2*zmid)/(2*K*Tmid)
        MPP[i,i+1] = 1/dz + (-Fmid*kram_es/cgs_c + omega**2*zmid)/(2*K*Tmid)
        MPT[i,i]   = -Pmid*(-Fmid*kram_es/cgs_c + omega**2*zmid)/(2*K*Tmid**2)
        MPT[i,i+1] = -Pmid*(-Fmid*kram_es/cgs_c + omega**2*zmid)/(2*K*Tmid**2)
        MPF[i,i]   = -Pmid*kram_es/(2*K*Tmid*cgs_c)
        MPF[i,i+1] = -Pmid*kram_es/(2*K*Tmid*cgs_c)
        AT [i]     = 3*Fmid*Pmid*kram_es/(16*K*Tmid**4*cgs_stef) + dT/dz
        MTP[i,i]   = 3*Fmid*kram_es/(32*K*Tmid**4*cgs_stef)
        MTP[i,i+1] = 3*Fmid*kram_es/(32*K*Tmid**4*cgs_stef)
        MTT[i,i]   = -3*Fmid*Pmid*kram_es/(8*K*Tmid**5*cgs_stef) - 1/dz
        MTT[i,i+1] = -3*Fmid*Pmid*kram_es/(8*K*Tmid**5*cgs_stef) + 1/dz
        MTF[i,i]   = 3*Pmid*kram_es/(32*K*Tmid**4*cgs_stef)
        MTF[i,i+1] = 3*Pmid*kram_es/(32*K*Tmid**4*cgs_stef)
        AF [i]     = -alpha*omega*(Pmid + 4*Tmid**4*cgs_stef/(3*cgs_c)) + dF/dz
        MFP[i,i]   = -alpha*omega/2
        MFP[i,i+1] = -alpha*omega/2
        MFT[i,i]   = -8*Tmid**3*alpha*cgs_stef*omega/(3*cgs_c)
        MFT[i,i+1] = -8*Tmid**3*alpha*cgs_stef*omega/(3*cgs_c)
        MFF[i,i]   = -1/dz
        MFF[i,i+1] = 1/dz
        if i == 0:
            AB1       = F[i]
            MB1F[i]   = 1
        if i == N-2:
            AB2       = -flux_acc + F[i + 1]
            MB2F[i+1] = 1
            AB3       = -2*cgs_stef*T[i + 1]**4 + F[i + 1]
            MB3F[i+1] = 1
            MB3T[i+1] = -8*cgs_stef*T[i + 1]**3

    np.save('M',M)
    np.save('A',A)
    delt = np.linalg.solve(M,-A)
    return (delt[0:N],delt[N:2*N],delt[2*N:3*N])

