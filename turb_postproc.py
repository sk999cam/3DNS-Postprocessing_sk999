import os
import gc
import array
import struct
import numpy as np

#This module is used to write to paraview-readable .vtk file from python.
#If you do not want to write to vtk, you will not need this module
#In the past I have had some trouble installing it, so feel free to not use if it's too involved getting it to run.
# from evtk.hl import gridToVTK

import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator

from spec_nonun import *

def conditionRadSlc(namePath,gridPath,DIM,makePrdc,bndrylr):

    #Conditions radial slice for use in 4DNS.
    # You will probably not need to use this function.
    # Input arguments
    #namePath: path to + name of turbulence slice
    #gridPath: path to grid the slice was sampled on
    #DIM: dimension of slice in points (integers)
    #makePrdc: logical, indicator whether boundaries should be made periodic.
    #bndrylr: logoical, indicator whether to add a boundary layer at hub and casing.

    fC = namePath
    c = array.array('f')
    c.fromfile(open(fC, 'rb'), os.path.getsize(fC) // c.itemsize)

    C = np.frombuffer(c, dtype=np.float64)  
    C = np.reshape(C,(DIM[0],DIM[1],DIM[2],3))


    u = C[:,:,:,0]
    v = C[:,:,:,1]
    w = C[:,:,:,2]
    
    plt.pcolor(v[1,:,:])
    plt.show()
    ni = int(DIM[0])
    nj = int(DIM[1])
    nk = int(DIM[2])


    fG = gridPath
    g = array.array('d')
    g.fromfile(open(fG, 'rb'), os.path.getsize(fG) // g.itemsize)
    G = np.frombuffer(g, dtype=np.float64)   
    G = np.transpose(G)
    G = np.reshape(G,(3,ni,nj,nk),order ='F')
    x = G[0,:,:,:]
    y = G[1,:,:,:]
    z = G[2,:,:,:]

    
    t = np.arctan2(y, z)
        
    vr = np.cos(t)*w+np.sin(t)*v
    vt = np.cos(t)*v-np.sin(t)*w
    

    if bndrylr == 1:


        print('Dampening fluctuations towards wall')
        # add endwall boundary layers
        r = np.sqrt(z*z + y*y)
        r = np.squeeze(r[0,0,:])
        
        fspan = (r  - np.min(r))/(np.max(r) - np.min(r))
        d = 0.01

        ffspan = 1 - fspan
        ihub = np.argwhere(fspan/d < 1)
        icas = np.argwhere(ffspan/d < 1)
        
        
        fbl = np.ones([nk])
        fbl[ihub] = -(2/3)*(fspan[ihub]/d)**(2.5) + (5/3)*fspan[ihub]/d
        fbl[icas] = -(2/3)*(ffspan[icas]/d)**(2.5) + (5/3)*ffspan[icas]/d
        fbl = np.expand_dims(fbl,axis=0)
        fbl = np.expand_dims(fbl,axis=0)
        #print(np.shape(fbl))

        u = u*np.tile(fbl,(ni,nj,1))
        vr = vr*np.tile(fbl,(ni,nj,1))
        vt = vt*np.tile(fbl,(ni,nj,1))

    if (makePrdc==1):
    
        print('Making box bounadries periodic')
       
        #enforcing periodicity on the j-index only (theta index in polar)
        u[:,0,:] = 1/2*u[:,1,:]+1/2*u[:,nj-2,:]
        #u[:,:,0] = 1/2*u[:,:,1]+1/2*u[:,:,nk-2]

        #u[ni,:,:] = 1/2*u[1,:,:]+1/2*u[ni-1,:,:]
        u[:,nj-1,:] = 1/2*u[:,1,:]+1/2*u[:,nj-2,:]
        #u[:,:,nk-1] = 1/2*u[:,:,1]+1/2*u[:,:,nk-2]

        #v[0,:,:] = 1/2*v[1,:,:]+1/2*v[ni-1,:,:]
        vt[:,0,:] = 1/2*vt[:,1,:]+1/2*vt[:,nj-2,:]
        #vt[:,:,0] = 1/2*vt[:,:,1]+1/2*vt[:,:,nk-2]

        #v[ni,:,:] = 1/2*v[1,:,:]+1/2*v[ni-1,:,:]
        vt[:,nj-1,:] = 1/2*vt[:,1,:]+1/2*vt[:,nj-2,:]
        #vt[:,:,nk-1] = 1/2*vt[:,:,1]+1/2*vt[:,:,nk-2]

        #w[0,:,:] = 1/2*w[1,:,:]+1/2*w[ni-1,:,:]
        vr[:,0,:] = 1/2*vr[:,1,:]+1/2*vr[:,nj-2,:]
        #vr[:,:,0] = 1/2*vr[:,:,1]+1/2*vr[:,:,nk-2]

        #w[ni,:,:] = 1/2*w[1,:,:]+1/2*w[ni-1,:,:]
        vr[:,nj-1,:] = 1/2*vr[:,1,:]+1/2*vr[:,nj-2,:]
        #vr[:,:,nk-1] = 1/2*vr[:,:,1]+1/2*vr[:,:,nk-2]
    else:
        print('Box bounadries not periodic')

    v = np.cos(t)*vt+np.sin(t)*vr
    w = np.cos(t)*vr-np.sin(t)*vt

    
    return [u,v,w]

def conditionAxSlc(namePath,gridPath,DIM,makePrdc,calcSpc,SpcName):

    from spec_nonun import smallCubeSpec as sCS

    #Conditions axial slice / box for use in 3DNS.
    # You will probably want to use this function.
    # Input arguments
    #namePath: path to + name of turbulence slice
    #gridPath: path to grid the slice was sampled on
    #DIM: dimension of slice in points (integers)
    #makePrdc: logical, indicator whether boundaries should be made periodic.
    #bndrylr: logoical, indicator whether to add a boundary layer at hub and casing.
    #calcSpc: logoical, indicator whether to calculate spectra of largest possible cube within inlet box
    #SpcName: path+name of .mat file to save spectrum to. You will need to uncomment this lower down in the function if write-out is desired.
    
  
    fC = namePath
    c = array.array('f')
    c.fromfile(open(fC, 'rb'), os.path.getsize(fC) // c.itemsize)

    C = np.frombuffer(c, dtype=np.float64)  
    C = np.reshape(C,(DIM[0],DIM[1],DIM[2],3))


    u = C[:,:,:,0]
    v = C[:,:,:,1]
    w = C[:,:,:,2]
    
    #plt.pcolor(v[1,:,:])
    #plt.axis('equal')
    #plt.show()
    ni = int(DIM[0])
    nj = int(DIM[1])
    nk = int(DIM[2])


    fG = gridPath
    g = array.array('d')
    g.fromfile(open(fG, 'rb'), os.path.getsize(fG) // g.itemsize)
    G = np.frombuffer(g, dtype=np.float64)   
    G = np.transpose(G)
    G = np.reshape(G,(3,ni,nj,nk),order ='F')
    x = G[0,:,:,:]
    y = G[1,:,:,:]
    z = G[2,:,:,:]

    if (makePrdc==1):
    
        print('Making box bounadries periodic')
       
        #enforcing periodicity in all directions
        u[0,:,:] = 1/2*u[1,:,:]+1/2*u[ni-2,:,:]
        u[:,0,:] = 1/2*u[:,1,:]+1/2*u[:,nj-2,:]
        u[:,:,0] = 1/2*u[:,:,1]+1/2*u[:,:,nk-2]

        u[ni-1,:,:] = 1/2*u[1,:,:]+1/2*u[ni-2,:,:]
        u[:,nj-1,:] = 1/2*u[:,1,:]+1/2*u[:,nj-2,:]
        u[:,:,nk-1] = 1/2*u[:,:,1]+1/2*u[:,:,nk-2]

        v[0,:,:] = 1/2*v[1,:,:]+1/2*v[ni-1,:,:]
        v[:,0,:] = 1/2*v[:,1,:]+1/2*v[:,nj-2,:]
        v[:,:,0] = 1/2*v[:,:,1]+1/2*v[:,:,nk-2]

        v[ni-1,:,:] = 1/2*v[1,:,:]+1/2*v[ni-2,:,:]
        v[:,nj-1,:] = 1/2*v[:,1,:]+1/2*v[:,nj-2,:]
        v[:,:,nk-1] = 1/2*v[:,:,1]+1/2*v[:,:,nk-2]

        w[0,:,:] = 1/2*w[1,:,:]+1/2*w[ni-2,:,:]
        w[:,0,:] = 1/2*w[:,1,:]+1/2*w[:,nj-2,:]
        w[:,:,0] = 1/2*w[:,:,1]+1/2*w[:,:,nk-2]

        w[ni-1,:,:] = 1/2*w[1,:,:]+1/2*w[ni-2,:,:]
        w[:,nj-1,:] = 1/2*w[:,1,:]+1/2*w[:,nj-2,:]
        w[:,:,nk-1] = 1/2*w[:,:,1]+1/2*w[:,:,nk-2]
    else:
        print('Box bounadries not periodic')

    
    if (calcSpc==1):

        boxL = [0,0,0]
        boxL[0] = np.abs(x[-1,0,0] - x[0,0,0])
        boxL[1] = np.abs(y[0,-1,0] - y[0,0,0])
        boxL[2] = np.abs(x[0,0,-1] - z[0,0,0])
    
        print('Calculating 3D sub-cube spectra')
        [wn, spec] = sCS(u,v,w, boxL)
        
        #write to file if desired
        #scipy.io.savemat(SpcName,dict(E=spec,k=wn))
      
    return [u,v,w]

def checkTu(namePath,cubeShp,vref):
    
    #Can use this to check whether your sampled box has the inlet turbulence you want from it.

    from scipy.interpolate import LinearNDInterpolator
    
    [uc,vc,wc] = read3DNSorder(namePath,[cubeShp,cubeShp,cubeShp])
    
    Tu_u = np.sqrt(np.mean((uc**2).flatten()))/vref
    Tu_v = np.sqrt(np.mean((vc**2).flatten()))/vref
    Tu_w = np.sqrt(np.mean((wc**2).flatten()))/vref

    print(Tu_u, Tu_v, Tu_w)
    
    return

def write3DNSorder(namePath,u,v,w):

    #writes velocity out the same way Andy writes and reads in in 3DNS (I believe).
    #namePath: string, path+name of file velocity is to be written to.
    #u,v,w: 3D arrays of floats holding velocity components.

    u = np.swapaxes(u,1,2)
    v = np.swapaxes(v,1,2)
    w = np.swapaxes(w,1,2)

    u = np.expand_dims(u,3)
    v = np.expand_dims(v,3)
    w = np.expand_dims(w,3)

    vel = np.concatenate((u,v,w), axis = 3)

    fvel = vel.flatten()
    if np.any(np.isnan(fvel)):
        print('nans in flow')
    else:
        print('no nans in flow')
    output_file = open(namePath, 'wb')
    vel.tofile(output_file)
    output_file.close()
    
    return

def read3DNSorder(namePath,floShp):

    fC = namePath
    c = array.array('f')
    c.fromfile(open(fC, 'rb'), os.path.getsize(fC) // c.itemsize)

    C = np.frombuffer(c, dtype=np.float64)  
    C = np.reshape(C,(floShp[0],floShp[2],floShp[1],3))


    u = C[:,:,:,0]
    v = C[:,:,:,1]
    w = C[:,:,:,2]
    u = np.swapaxes(u,1,2)
    v = np.swapaxes(v,1,2)
    w = np.swapaxes(w,1,2)

    return u, v, w
    
def testRoundTrip_P1(namePath,cubeShp,testPath):
    #You will be unlikely to need to use this.
    
    #This is round-trip testing used together with turbTest.f
    #Use this routine as a sanity check that writing and reading are done in a way that 3DNS
    # gets delivered the inflow turbulence correctly (hopefully)
    
    fC = namePath

    c = array.array('f')
    c.fromfile(open(fC, 'rb'), os.namePath.getsize(fC) // c.itemsize)

    C = np.frombuffer(c, dtype=np.float64)  
    C = np.reshape(C,(cubeShp,cubeShp,cubeShp,3))

    uc = C[:,:,:,0]
    vc = C[:,:,:,1]
    wc = C[:,:,:,2]
    
    write3DNSorder(testPath,uc,vc,wc)
    
    #plt.figure
    #plt.pcolor(uc[:,77,:])
    #plt.gca().set_aspect('equal', adjustable='box')
    #plt.colorbar()
    #plt.show()
    
    print(u[0,0,0], u[100,200,5])
    print(v[0,0,0], v[100,200,5])
    print(w[0,0,0], w[100,200,5])
    
    return

def testRoundTrip_P2(cubeShp,outPath):

    #Again, you should be able to just ignore this, unless you want to test that I'm writing things out correctly.
    #Part 2 of round-trip testing
    fC = outPath

    c = array.array('f')
    c.fromfile(open(fC, 'rb'), os.path.getsize(fC) // c.itemsize)

    C = np.frombuffer(c, dtype=np.float64)  
    C = np.reshape(C,(cubeShp,cubeShp,cubeShp,3))

    uc = C[:,:,:,0]
    vc = C[:,:,:,1]
    wc = C[:,:,:,2]
    
    #plt.figure
    #plt.pcolor(uc[:,77,:])
    #plt.gca().set_aspect('equal', adjustable='box')
    #plt.colorbar()
    #plt.show()
    
    print(u[0,0,0], u[100,200,5])
    print(v[0,0,0], v[100,200,5])
    print(w[0,0,0], w[100,200,5])

    #To test C-order writing as used in turbTest2.f
    #u2 = np.expand_dims(u,3)
    #v2 = np.expand_dims(v,3)
    #w2 = np.expand_dims(w,3)

    #vel = np.concatenate((u2,v2,w2), axis = 3)

    #output_file = open('../../DNS/test-Tu5.dat', 'wb')
    #vel.tofile(output_file)
    #output_file.close()
    
    return
    
    

