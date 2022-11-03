#!/usr/bin/env python
# Converted by Phillip J. Wolfram from a matlab script authored by Todd Ringler, 2014
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from numpy import linalg as LA

def baroclinicModes_function(z, N, latitude):
    """
    intent ins
    # assumes that z is negative and increasing (m)
    # N is the BruntVaisala Frequency in 1/s
    # latitude is in degrees (N is positive)
    # intent outs

    intent outs
    # ix sorted indicies of rrd eigenvalues
    # rrd Rossby radius of deformation (km)
    # V eigenvectors
    # i sorted indicies of eigenvectors
    """
    if np.isnan(np.sum(N)):
        return np.nan*N, None, None, None, None

    #specify number of vertical layers
    #++++ set to correct # of vertical levels
    nVertLevels=len(z)
    grav = 9.81
    omega = 2*np.pi/86400
    f = 2*omega*np.sin(latitude*np.pi/180.0)
    zBottom = z[-1]
    H0=zBottom
    #++++

    #specify the number of vertical layer edges
    nVertLevelsP1=nVertLevels+1
    nVertLevelsM1=nVertLevels-1

    dz = zBottom / nVertLevels

    # interpolate bvf
    zMid = np.arange(dz/2, zBottom, dz)
    interpSpline = UnivariateSpline(z[-1::-1], N[-1::-1], k=1, s=0)
    bvf = interpSpline(zMid)

    #allocate variables defined at layer centers
    #we will solve M x=\lambda x
    M = np.zeros((nVertLevelsM1,nVertLevelsM1))

    #++++need to generized for variable dz
    #++++or just interpolate bvf to uniform grid might be better to start
    for k in np.arange(1,nVertLevelsM1-1):
        M[k,k-1]=+1/(dz**2)
        M[k,k  ]=-2/(dz**2)
        M[k,k+1]=+1/(dz**2)

    M[0,0]=-2/(dz**2)
    M[0,1]=+1/(dz**2)

    M[-1,-2]=+1/(dz**2)
    M[-1,-1]=-2/(dz**2)

    for k in np.arange(nVertLevelsM1):
        M[k,:] = M[k,:] / ( (bvf[k]+bvf[k+1])/2 )**2

    #solve the eigenvalue problem
    try:
        [E,V]=LA.eig(M)
    except:
        print 'latititue=', latitude
        print 'z=', z
        print 'N=', N
        return

    b = np.sqrt(-1.0/E)
    i = np.arange(nVertLevels-1)

    ix = np.argsort(b)
    ix = ix[::-1]
    c = b[ix]
    rrd = c / f / 1000.

    V *= -1.0 # for consistency with matlab since only shape matters
    H = np.zeros((nVertLevelsM1,nVertLevelsM1))
    for iMode in np.arange(nVertLevelsM1):
        H[nVertLevelsM1-1,iMode]=0.0
        sum=0
        for k in np.arange(nVertLevelsM1-2,-1,-1):
            work = (bvf[k+1]*bvf[k+1])*grav*dz
            H[k,iMode] = H[k+1,iMode] + V[k+1,iMode]*work
            sum = sum + H[k,iMode]*dz
        H[:,iMode]=H[:,iMode]-sum/H0

    return rrd, ix, V, i, H


def test_baroclinic_modes():
    figType = '.jpeg'

    #specify number of vertical layers
    #++++ set to correct # of vertical levels
    nVertLevels=250
    #++++

    #specify the number of vertical layer edges
    nVertLevelsP1=nVertLevels+1
    nVertLevelsM1=nVertLevels-1

    #specify the bottom and find layer thickness
    #++++set bottom depth correctly
    zBottom = -4000.0
    #++++
    dz = -zBottom / nVertLevels
    H0=-zBottom

    #define reference values
    #----
    density_ref=1030.0
    density_delta=1.5
    density_linear=0.05
    #----
    grav = 9.81
    omega = 2*np.pi/86400
    #++++
    f = 2*omega*np.sin(55.0*np.pi/180.0)
    #++++

    #allocate variables defined at layer centers
    #we will solve M x=\lambda x
    M = np.zeros((nVertLevelsM1,nVertLevelsM1))
    bvf = np.zeros((nVertLevels,))
    zMid = np.zeros((nVertLevels,))
    #allocate variables defined at layer interfaces
    z = np.zeros((nVertLevelsP1,))
    density = np.zeros((nVertLevelsP1,))

    #compute the z-level and density at layer edges
    for k in np.arange(nVertLevelsP1):
        #++++ define correct z here
        z[k]=-dz*k
        density[k]=density_ref + density_delta * z[k] / zBottom
        density[k]=density_ref - \
                (1.0-density_linear)*density_delta*np.tanh(z[k]/300.0) - \
                          density_linear*density_delta*z[k] / H0

    #compute the buoyancy frequency
    for k in np.arange(nVertLevels):
        zMid[k] = (z[k]+z[k+1])/2.0
        bvf[k]= np.sqrt(-grav/density_ref*(density[k]-density[k+1])/dz)

    rrd, ix, V, i, H = baroclinicModes_function(0.5*(z[:-1]+z[1:]), bvf, 55.0)

    figName='RRD'
    plt.figure(2)
    plt.scatter(i[:11],rrd[:11], color='black', linewidth=4)
    plt.xlim((0,10))
    plt.title('Spectrum of Rossby Radii of Deformation')
    plt.xlabel('Eigenvalues Sorted in Descending Order')
    plt.ylabel('Rossby Radius (km)')
    plt.savefig(figName+figType)

    figName='verticalEigenvectors'
    plt.figure(3)
    plt.plot( V[:,ix[0]], z[1:-1], color='black', linewidth=4 )
    plt.hold(True)
    plt.plot( V[:,ix[1]], z[1:-1], color='blue', linewidth=3 )
    plt.plot( V[:,ix[2]], z[1:-1], color='red', linewidth=2 )
    plt.plot( V[:,ix[3]], z[1:-1], color='green', linewidth=1 )
    plt.plot( V[:,ix[4]], z[1:-1], color='yellow', linewidth=0.5 )
    plt.hold(False)
    plt.title('Structure of Vertical Eigenvectors')
    plt.xlabel('Amplitude')
    plt.ylabel('z coordinate (m)')
    plt.savefig(figName+figType)

    figName='densityProfile'
    plt.figure(4)
    plt.plot( density, z)
    plt.title('Density Profile')
    plt.xlabel('Density (kg/m3)')
    plt.ylabel('z coordinate (m)')
    plt.savefig(figName+figType)

    figName='bvf2'
    plt.figure(8)
    plt.plot( np.log10(bvf*bvf), zMid)
    plt.title('log10(N^2)')
    plt.xlabel('log10(N^2) (1/s^2)')
    plt.ylabel('z coordinate (m)')
    plt.savefig(figName+figType)

    figName='horizontalEigenvectors'
    plt.figure(1)
    plt.plot( H[:,ix[0]], z[1:-1], color='black', linewidth=4 )
    plt.hold(True)
    plt.plot( H[:,ix[1]], z[1:-1], color='blue',  linewidth=3 )
    plt.plot( H[:,ix[2]], z[1:-1], color='red',   linewidth=2 )
    plt.plot( H[:,ix[3]], z[1:-1], color='green', linewidth=1 )
    plt.plot( H[:,ix[4]], z[1:-1], color='yellow',linewidth=0.5 )
    plt.hold(False)
    plt.title('Structure of Horizontal Eigenvectors')
    plt.xlabel('Amplitude')
    plt.ylabel('z coordinate (m)')
    plt.savefig(figName+figType)

    plt.show()

if __name__=="__main__":
    test_baroclinic_modes()


