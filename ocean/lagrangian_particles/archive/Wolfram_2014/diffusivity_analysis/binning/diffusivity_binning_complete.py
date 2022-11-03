#!/usr/bin/env python
"""

    Compute diffusivity, eddy speed, and decorrelation via binning procedure

    Phillip Wolfram
    LANL
    04/01/2015

"""
## plot on mustang
#import matplotlib as mpl
#mpl.use('Agg')
# import libraries / packages
import netCDF4
import xarray as xr
import numpy as np
from file_handling import check_folder, savenpz
from latlon_coordinate_transforms import signed_distances_numexpr as signed_distances, haversine_formula_numexpr as haversine_formula, fix_periodicity_numexpr as fix_periodicity
from compute_diffusivity import compute_dispersion_numexpr, compute_dispersion_planar_numexpr

def build_particle_file(particledata, samplepoints, radius=1.0e5, Nlen=30, layerrange=np.arange(0,11), degrees=False, planar=False, fraction=1.0,filteroutcropped=True, L=1.e6):  #{{{

    # load the file database (for a single realization that could span multiple files)
    ds = xr.open_mfdataset(particledata, concat_dim='Time')
    if degrees:
        lonname = 'lon'
        latname = 'lat'
        conv = np.radians
    else:
        lonname = 'lon'
        latname = 'lat'
        def creturn(c):
            return c
        conv = creturn

    rEarth = 6371220. # from netcdf file

    data  = netCDF4.Dataset(samplepoints,'r')
    if planar:
        x = data.variables['xCell'][:]
        y = data.variables['yCell'][:]
    else:
        x = data.variables['lonCell'][:] - 2*np.pi*(data.variables['lonCell'][:] > np.pi)
        y = data.variables['latCell'][:]
    Nc = x.shape[0]
    Nlen = len(ds.time)
    buoyancySurfaceVec = ds.Nb

    case = 'CaseBinningEnsembledComplete_r=%0.2f_Nlen=%d_Nc=%d_frac=%0.2f/'%(radius, Nlen, Nc, fraction)

    print 'Case is ' + case + ' for layerrange ', layerrange
    #check_folder(case)

    ## assumes small file size (should only operate on a single contiguous realization)
    lat = ds.lat.values
    lon = ds.lon.values
    dtdays = ds.dtdays.values
    dtdays = np.swapaxes(dtdays[:,:Nlen],0,1)[np.newaxis,:,:]
    notoutcropped = ds.notoutcropped.values

    # allocate memory #{{{
    Nt = Nlen-1
    Nb = len(layerrange)
    Nr = len(ds['Nr'])
    clat      = np.nan*np.zeros((Nc, Nlen, Nr, Nb))
    clon      = np.nan*np.zeros((Nc, Nlen, Nr, Nb))
    drdr      = np.nan*np.zeros((Nc, Nlen, Nr, Nb))
    dxdx      = np.nan*np.zeros((Nc, Nlen, Nr, Nb))
    dxdy      = np.nan*np.zeros((Nc, Nlen, Nr, Nb))
    dydy      = np.nan*np.zeros((Nc, Nlen, Nr, Nb))
    absdrdr   = np.nan*np.zeros((Nc, Nlen, Nr, Nb))
    absdxdx   = np.nan*np.zeros((Nc, Nlen, Nr, Nb))
    absdxdy   = np.nan*np.zeros((Nc, Nlen, Nr, Nb))
    absdydy   = np.nan*np.zeros((Nc, Nlen, Nr, Nb))
    K_rr      = np.nan*np.zeros((Nc, Nt, Nr, Nb))
    K_xx      = np.nan*np.zeros((Nc, Nt, Nr, Nb))
    K_xy      = np.nan*np.zeros((Nc, Nt, Nr, Nb))
    K_yy      = np.nan*np.zeros((Nc, Nt, Nr, Nb))
    absK_rr   = np.nan*np.zeros((Nc, Nt, Nr, Nb))
    absK_xx   = np.nan*np.zeros((Nc, Nt, Nr, Nb))
    absK_xy   = np.nan*np.zeros((Nc, Nt, Nr, Nb))
    absK_yy   = np.nan*np.zeros((Nc, Nt, Nr, Nb))
    npart     = np.nan*np.zeros((Nc))
    meanu     = np.nan*np.zeros((Nc, Nt, Nr, Nb))
    eddyu     = np.nan*np.zeros((Nc, Nt, Nr, Nb))
    meanv     = np.nan*np.zeros((Nc, Nt, Nr, Nb))
    eddyv     = np.nan*np.zeros((Nc, Nt, Nr, Nb))
    meanspeed = np.nan*np.zeros((Nc, Nt, Nr, Nb))
    eddyspeed = np.nan*np.zeros((Nc, Nt, Nr, Nb))
    up0up     = np.nan*np.zeros((Nc, Nt, Nr, Nb))
    up0vp     = np.nan*np.zeros((Nc, Nt, Nr, Nb))
    vp0vp     = np.nan*np.zeros((Nc, Nt, Nr, Nb))
    vp0up     = np.nan*np.zeros((Nc, Nt, Nr, Nb))
    up0up0    = np.nan*np.zeros((Nc, Nr, Nb))
    vp0vp0    = np.nan*np.zeros((Nc, Nr, Nb))
    #}}}

    # starting center points of cluster
    for nlayer in layerrange:
        nlayerid = nlayer - layerrange[0]

        for acluster in np.arange(Nc):
            if np.mod(acluster,np.round(Nc/100)) == 0:
                print 'on %d of %d'%(acluster, Nc)
            # compute particle clusters
            if planar:
              clusterpart = np.sqrt(\
                  (y[acluster] - lat[0,0,nlayer,:])**2.0 + \
                  (x[acluster] - fix_periodicity(lon[0,0,nlayer,:], x[acluster], L))**2.0\
                  ) < radius
            else:
              clusterpart = haversine_formula(y[acluster], lat[0,0,nlayer,:], x[acluster], lon[0,0,nlayer,:], rEarth) < radius

            # subsample particles
            clustids = np.where(clusterpart)[0]
            Ncorig = np.sum(clusterpart)
            # remove non-desired random fraction from cluster
            removeids = np.random.choice(clustids, int(np.floor(clustids.shape[0]*(1-fraction))), replace=False)
            clusterpart[removeids] = False
            # checks to make sure subsampling was correct
            Ncsample = np.sum(clusterpart)
            assert np.abs(Ncsample - np.round(fraction*Ncorig)) < 2, 'Error with subsampling'

            llon = np.swapaxes(lon[:,:Nlen,nlayer,clusterpart],0,2)
            llat = np.swapaxes(lat[:,:Nlen,nlayer,clusterpart],0,2)

            # mask particles that have outcropped with nans
            if filteroutcropped:
              mask = np.swapaxes(notoutcropped[:,:Nlen,nlayer,clusterpart],0,2)
              # make first time step count for realizations
              mask[:,:,:] = mask[:,0,:][:,np.newaxis,:]
              # this just removes points that are directly in the outcropping layers
              notvalid = np.logical_not(mask)
              llon[notvalid] = np.nan
              llat[notvalid] = np.nan

            # compute velocity
            if planar:
                velx = fix_periodicity(llon[:,1:,:], llon[:,:-1,:], L) - llon[:,:-1,:]
                vely = llat[:,1:,:] - llat[:,:-1,:]
            else:
                velx, vely = signed_distances(llat[:,1:,:], llat[:,:-1,:], \
                                              llon[:,1:,:], llon[:,:-1,:], rEarth)

            velx /= (24.*60.*60.*dtdays)
            vely /= (24.*60.*60.*dtdays)

            # compute the eddy speed and mean speed
            meanu[acluster,:,:,nlayerid] = np.mean(velx,axis=0)
            meanv[acluster,:,:,nlayerid] = np.mean(vely,axis=0)
            eddyu[acluster,:,:,nlayerid] = np.std(velx,axis=0)
            eddyv[acluster,:,:,nlayerid] = np.std(vely,axis=0)
            meanspeed[acluster,:,:,nlayerid] = np.sqrt(meanu[acluster,:,:,nlayerid]**2.0 + meanv[acluster,:,:,nlayerid]**2.0)
            eddyspeed[acluster,:,:,nlayerid] = np.sqrt(eddyu[acluster,:,:,nlayerid]**2.0 + eddyv[acluster,:,:,nlayerid]**2.0)

            # compute decorrelation coefficient rho (Nc, Nt).  Note (Np,Nt,Ne) is dim of velx
            up = velx - meanu[acluster,:,:,nlayerid][np.newaxis,:,:]
            vp = vely - meanv[acluster,:,:,nlayerid][np.newaxis,:,:]
            # could generalize to full tensor if desired (need cross terms for the complex autocorrelation)
            up0up[acluster,:,:,nlayerid]  = np.mean(up[:,0,:][:,np.newaxis,:]*up, axis=0)
            up0vp[acluster,:,:,nlayerid]  = np.mean(up[:,0,:][:,np.newaxis,:]*vp, axis=0)
            vp0vp[acluster,:,:,nlayerid]  = np.mean(vp[:,0,:][:,np.newaxis,:]*vp, axis=0)
            vp0up[acluster,:,:,nlayerid]  = np.mean(vp[:,0,:][:,np.newaxis,:]*up, axis=0)
            up0up0[acluster,:,nlayerid]   = np.mean(up[:,0,:]*up[:,0,:],axis=0)
            vp0vp0[acluster,:,nlayerid]   = np.mean(vp[:,0,:]*vp[:,0,:],axis=0)

            # compute mean dispersion over particle clusters (result are size Ncluster, Ntime, Nb)
            if planar:
                npart[acluster], clat[acluster,:,:,nlayerid], clon[acluster,:,:,nlayerid], \
                        drdr[acluster,:,:,nlayerid], dxdx[acluster,:,:,nlayerid], \
                        dxdy[acluster,:,:,nlayerid], dydy[acluster,:,:,nlayerid], \
                        absdrdr[acluster,:,:,nlayerid], absdxdx[acluster,:,:,nlayerid], \
                        absdxdy[acluster,:,:,nlayerid], absdydy[acluster,:,:,nlayerid] = \
                        compute_dispersion_planar_numexpr(llat, fix_periodicity(llon, x[acluster], L), ensmean=0)
            else:
                npart[acluster], clat[acluster,:,:,nlayerid], clon[acluster,:,:,nlayerid], \
                        drdr[acluster,:,:,nlayerid], dxdx[acluster,:,:,nlayerid], \
                        dxdy[acluster,:,:,nlayerid], dydy[acluster,:,:,nlayerid], \
                        absdrdr[acluster,:,:,nlayerid], absdxdx[acluster,:,:,nlayerid], \
                        absdxdy[acluster,:,:,nlayerid], absdydy[acluster,:,:,nlayerid] = \
                        compute_dispersion_numexpr(llat, llon, rEarth, ensmean=0)

        # compute diffusivity from dispersion 1/2 * d/dt (variance)
        K_rr[:,:,:,nlayerid] = 0.5*np.diff(drdr[:,:,:,nlayerid],axis=1)/(24.*60.*60.*dtdays)
        K_xx[:,:,:,nlayerid] = 0.5*np.diff(dxdx[:,:,:,nlayerid],axis=1)/(24.*60.*60.*dtdays)
        K_xy[:,:,:,nlayerid] = 0.5*np.diff(dxdy[:,:,:,nlayerid],axis=1)/(24.*60.*60.*dtdays)
        K_yy[:,:,:,nlayerid] = 0.5*np.diff(dydy[:,:,:,nlayerid],axis=1)/(24.*60.*60.*dtdays)

        # absolute diffusivity
        absK_rr[:,:,:,nlayerid] = 0.5*np.diff(absdrdr[:,:,:,nlayerid],axis=1)/(24.*60.*60.*dtdays)
        absK_xx[:,:,:,nlayerid] = 0.5*np.diff(absdxdx[:,:,:,nlayerid],axis=1)/(24.*60.*60.*dtdays)
        absK_xy[:,:,:,nlayerid] = 0.5*np.diff(absdxdy[:,:,:,nlayerid],axis=1)/(24.*60.*60.*dtdays)
        absK_yy[:,:,:,nlayerid] = 0.5*np.diff(absdydy[:,:,:,nlayerid],axis=1)/(24.*60.*60.*dtdays)

        print 'finished %d layer with %d clusters' % (nlayer, Nc)

    # save data and make sure to document how files were produced
    xr.Dataset({'clon'      : (['Nc', 'Nt', 'Nr', 'Nb'], clon),        \
                'clat'      : (['Nc', 'Nt', 'Nr', 'Nb'], clat),        \
                'drdr'      : (['Nc', 'Nt', 'Nr', 'Nb'], drdr),        \
                'dxdx'      : (['Nc', 'Nt', 'Nr', 'Nb'], dxdx),        \
                'dxdy'      : (['Nc', 'Nt', 'Nr', 'Nb'], dxdy),        \
                'dydy'      : (['Nc', 'Nt', 'Nr', 'Nb'], dydy),        \
                'absdrdr'   : (['Nc', 'Nt', 'Nr', 'Nb'], absdrdr),     \
                'absdxdx'   : (['Nc', 'Nt', 'Nr', 'Nb'], absdxdx),     \
                'absdxdy'   : (['Nc', 'Nt', 'Nr', 'Nb'], absdxdy),     \
                'absdydy'   : (['Nc', 'Nt', 'Nr', 'Nb'], absdydy),     \
                'K_rr'      : (['Nc', 'Nt-1', 'Nr', 'Nb'], K_rr),      \
                'K_xx'      : (['Nc', 'Nt-1', 'Nr', 'Nb'], K_xx),      \
                'K_xy'      : (['Nc', 'Nt-1', 'Nr', 'Nb'], K_xy),      \
                'K_yy'      : (['Nc', 'Nt-1', 'Nr', 'Nb'], K_yy),      \
                'absK_rr'   : (['Nc', 'Nt-1', 'Nr', 'Nb'], absK_rr),   \
                'absK_xx'   : (['Nc', 'Nt-1', 'Nr', 'Nb'], absK_xx),   \
                'absK_xy'   : (['Nc', 'Nt-1', 'Nr', 'Nb'], absK_xy),   \
                'absK_yy'   : (['Nc', 'Nt-1', 'Nr', 'Nb'], absK_yy),   \
                'meanu'     : (['Nc', 'Nt-1', 'Nr', 'Nb'], meanu),     \
                'meanv'     : (['Nc', 'Nt-1', 'Nr', 'Nb'], meanv),     \
                'eddyu'     : (['Nc', 'Nt-1', 'Nr', 'Nb'], eddyu),     \
                'eddyv'     : (['Nc', 'Nt-1', 'Nr', 'Nb'], eddyv),     \
                'meanspeed' : (['Nc', 'Nt-1', 'Nr', 'Nb'], meanspeed), \
                'eddyspeed' : (['Nc', 'Nt-1', 'Nr', 'Nb'], eddyspeed), \
                'up0up'     : (['Nc', 'Nt-1', 'Nr', 'Nb'], up0up),     \
                'up0vp'     : (['Nc', 'Nt-1', 'Nr', 'Nb'], up0vp),     \
                'vp0vp'     : (['Nc', 'Nt-1', 'Nr', 'Nb'], vp0vp),     \
                'vp0up'     : (['Nc', 'Nt-1', 'Nr', 'Nb'], vp0up),     \
                'up0up0'    : (['Nc', 'Nr', 'Nb'], up0up0),            \
                'vp0vp0'    : (['Nc', 'Nr', 'Nb'], vp0vp0),            \
                'dtdays'    : (['Nr','Nt-1'], ds.dtdays.values),       \
                'npart'     : (['Nc'], npart)                          \
                },\
                coords={'rlzn'       : (['Nr'], ds.rlzn.values),                    \
                        'Nb'         : buoyancySurfaceVec[layerrange],              \
                        'xcluster'   : (['Nc'], x),                                 \
                        'ycluster'   : (['Nc'], y),                                 \
                        'time'       : (['Nt'], ds.time.values[:Nlen]),             \
                        'yearoffset' : ds.yearoffset.values                         \
                        }                                                           \
                )\
                .to_netcdf('dispersion_calcs_rlzn%04d_layerrange_%04d-%04d.nc'%(ds.rlzn.values,np.min(layerrange),np.max(layerrange)))

    # close the files
    ds.close()
    data.close()

    return  case #}}}

if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser() #{{{
    parser.add_option("-f", "--file", dest="inputfilename",
                      help="file to open for appending \
                      particle data 'particle_' extension or for analysis",
                      metavar="FILE")
    parser.add_option("-s", "--samplepoints", dest="samplepoints", help="npz file with x and y sample points", metavar="FILE")
    parser.add_option("-r", "--radius", dest="radius",
            help="radius of the cluster in km", metavar="FLOAT")
    parser.add_option("-d", "--degrees", dest="degrees", help="Data in degrees? T or F", metavar="BOOL")
    parser.add_option("-l", "--layer", dest="layer", help="Layer number", metavar="INT")
    parser.add_option("-p", "--planar", dest="planar", help="if on a plane", metavar="INT")
    parser.add_option("-n", "--numlen", dest="nlen",
            help="number of dtau derivatives", metavar="FLOAT")
    parser.add_option("--frac", dest="fraction", help="Fraction of particles in radius to use", metavar="FLOAT")

    options, args = parser.parse_args()

    if not options.inputfilename:
        parser.error("Input filename ('-f') is a required input... e.g., -f all_rlzn_particle_data.npz")
    if not options.samplepoints:
        parser.error("Sample points input filename ('-s') is a required input... e.g., -s coarse_sampling.npz")
    if not options.layer:
      options.layer = np.arange(5)
    else:
      try:
        options.layer = eval(options.layer)
      except:
        options.layer = np.array([int(options.layer)])
    if not options.radius:
        parser.error("Need starting number for radius.")
    if not options.nlen:
        parser.error("Need number of points for dtau derivatives.")
    if not options.degrees:
        print 'Warning: assuming data is in radians!'
        options.degrees = False
    else:
        if options.degrees == 'T':
            options.degrees = True
        elif options.degrees == 'F':
            options.degrees = False
        else:
            parser.error('Specify degrees as "T" or "F"')
    if not options.planar:
        print 'Warning: assuming data is not planar!'
        options.planar = False
    else:
        if options.planar == 'T':
            options.planar = True
        elif options.planar == 'F':
            options.planar = False
        else:
            parser.error('Specify planar as "T" or "F"')
    if not options.fraction:
        options.fraction = 1.0
    else:
        options.fraction = float(options.fraction)
    #}}}
    print 'starting analysis'
    options.inputfilename = build_particle_file(options.inputfilename, options.samplepoints, float(options.radius), int(options.nlen),
            layerrange=options.layer, degrees=options.degrees, planar=options.planar, fraction=options.fraction) #, \
#                filteroutcropped=False)
# option to not filter outcropped particles 

