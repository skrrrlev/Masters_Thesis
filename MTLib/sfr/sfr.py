from typing import Union
import numpy as np
import sympy
from latex2sympy2 import latex2sympy

def sargent(mstar: "Union[float, list[float], np.ndarray]", z: "float")-> "np.ndarray[float]":
    '''From Sargent et al., 2014'''
    mstar = np.array([mstar]).flatten() # solar masses

    N = 9.5E-11 * 10**(-0.21*np.log10(mstar/(5E10))) # yr^-1
    SFR = N * np.exp((2.05*z)/(1+0.16*z**1.54))*mstar # solar masses / yr
    
    return SFR

def sargent_uncertainty(mstar: "Union[float, list[float], np.ndarray]", z: "float", δmstar: "Union[float,list[float]]"):
    A, δA = 2.05, 0.265
    B, δB = 0.16, 0.11
    C, δC = 1.54, 0.32
    N, δN = 9.5E-11, 0.25E-11
    V, δV = -0.21, 0.04

    results = []
    mstar = np.array([mstar]).flatten()
    δmstar = np.array([δmstar]).flatten()
    for m,δm in zip(mstar,δmstar):
        K = 10**(V*np.log10(m*2E-11)) * np.exp((A*z)/(1+B*z**C))
        dfdM = lambda M,z: N * K * (V+1)
        dfdA = lambda M,z: N * K * M * z / (1+B*z**C)
        dfdB = lambda M,z: -N * K * M * A*z*z**C/ (1+B*z**C)**2
        dfdC = lambda M,z: -N * K * M * A*z*B*z**C*np.log(z)  / (1+B*z**C)**2
        dfdN = lambda M,z: K * M
        dfdV = lambda M,z:  N * K * M * np.log(M*2E-11)

        uncertainty = lambda M,z,δM: np.sqrt( (dfdM(M,z) * δM)**2 + (dfdA(M,z) * δA)**2 +(dfdB(M,z) * δB)**2 +(dfdC(M,z) * δC)**2 +(dfdN(M,z) * δN)**2 +(dfdV(M,z) * δV)**2 )
        results.append(uncertainty(m,z,δm))

    return results
    
def schreiber(mstar: "Union[float, list[float], np.ndarray]", z: "float")-> "np.ndarray[float]":
    '''From Schreiber et al., 2015'''
    mstar = np.array([mstar]).flatten() # solar masses

    m = np.log10(mstar/1E9)
    r = np.log10(1+z)
    s = np.array([max(0,mm-0.36-2.5*r) for mm in m])
    log_SFR = m-0.5+1.5*r-0.3*(s**2)

    SFR = 10**log_SFR
    
    return SFR

def schreiber_uncertainty(mstar: "Union[float,list[float]]", z: "float", δmstar: "Union[float,list[float]]"):
    # define constants
    A0, δA0 = 1.5, 0.15
    A1, δA1 = 0.3, 0.08
    A2, δA2 = 2.5, 0.6
    M0, δM0 = 0.5, 0.07
    M1, δM1 = 0.36, 0.3
    ln10 = np.log(10)
    results = []
    mstar = np.array([mstar]).flatten()
    δmstar = np.array([δmstar]).flatten()
    for m,δm in zip(mstar,δmstar):
        K = 10**( (np.log10(m*1E-9)) - M0 + A0*np.log10(z+1) - A1*(np.log10(m*1E-9) - M1 - A2*np.log10(z+1))**2 )
        dfdM = lambda M,z: -(K * ( -2*ln10*A1*M1-2*np.log(z+1)*A1*A2+2*np.log(M*1E-9)*A1-ln10)) / (ln10*M)
        dfda0 = lambda M,z:  K * np.log(z+1)
        dfda1 = lambda M,z: -K * ( np.log10(M*1E-9) - M1 - A2*np.log10(z+1) )**2 * ln10
        dfda2 = lambda M,z: 2 * K * A1 * ( np.log10(M*1E-9) - M1 - A2*np.log10(z+1)) * np.log(z+1)
        dfdm0 = lambda M,z: -K * ln10
        dfdm1 = lambda M,z: 2 * K * A1 * ( np.log10(M*1E-9) - M1 - A2*np.log10(z+1)) * ln10

        uncertainty = lambda M,z,δM: np.sqrt( (dfdM(M,z) * δM)**2 + (dfda0(M,z) * δA0)**2 + (dfda1(M,z) * δA1)**2 + (dfda2(M,z) * δA2)**2 + (dfdm0(M,z) * δM0)**2 + (dfdm1(M,z) * δM1)**2)
        results.append(uncertainty(m,z,δm))
    return results


if __name__ == '__main__':
    print('\nSchreiber: ')
    M, δM, z = 10**11, 10**10, 2
    print(schreiber(M,z))
    print(schreiber_uncertainty(M,z,δM))
    print('\nSargent: ')
    print(sargent(M,z))
    print(sargent_uncertainty(M,z,δM))