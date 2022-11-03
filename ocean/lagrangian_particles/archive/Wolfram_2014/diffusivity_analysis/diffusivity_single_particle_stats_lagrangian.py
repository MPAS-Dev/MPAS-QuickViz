#!/usr/bin/env python
"""

    Compute single particle statistics accross realizations.
    
    Phillip Wolfram
    LANL
    09/19/2014

"""
# import libraries / packages
import os
from subprocess import call
import warnings
import cPickle as pickle
import numpy as np
import netCDF4
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation, LinearTriInterpolator
from GeneralTriangulation import GeneralTriangulation
from scipy import spatial
from scipy.optimize import curve_fit
from scipy.integrate import trapz
import pylab
import glob

from diffusivity_realization_data import proj_lat_long

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':36})
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})


# function definitions
def get_realization_paths(rootdir, folder_prefix): #{{{
    """
    Get paths for realization folders
    Phillip Wolfram, LANL
    09/19/2014
    """
    fnames = []
    for root, dirs, files in os.walk(rootdir):
      if(root=='.'):
        dirs.sort()
        for adir in dirs:
          if(adir.startswith(folder_prefix)):
            fnames.append(rootdir + '/' + adir)

    return fnames #}}}

def get_num_ensembles(rootdir, folder_prefix): #{{{
    # get the number of realizations
    numRealizations = 0
    for root, dirs, files in os.walk(rootdir):
      if(root=='.'):
        dirs.sort()
        for adir in dirs:
          # for each folder starting with the prefix
          if(adir.startswith(folder_prefix)):
            numRealizations += 1
    return numRealizations #}}}

def open_netcdf_files(rlzn_path_list,name_prefix): #{{{
    """
    open each netcdf file for reading

    Phillip Wolfram
    LANL
    09/19/2014
    """

    fopen_list = []
    for path in rlzn_path_list:
        netcdf_name = glob.glob(path+'/'+name_prefix)
        fopen_list.append(netCDF4.Dataset(netcdf_name[0],'r'))

    return fopen_list #}}}

def close_netcdf_files(fopen_list): #{{{
    """
    close each netcdf file

    Phillip Wolfram
    LANL
    09/22/2014
    """
    for fopen in fopen_list:
        fopen.close()

    return #}}}

def compute_rlzn_ensemble(fopen_list, var_list, range_list): #{{{
    """
    Compute ensemble for each realization of the product of var_list 
    in the netCDF files opened in fopen_list.  range_list is used
    to constrain arrays to a particular value or range.

    This is basically an average over each realization.

    Phillip Wolfram 
    LANL
    09/19/2014
    """

    sumprod = 0.0

    for afile in fopen_list:
        prodvars = 1.0
        for avar,arange in zip(var_list,range_list):
            if isinstance(avar,str):
                # open value in netCDF4 database
                var_value = afile.variables[avar]
            else:
                var_value = avar
           
            if arange is None:
                if hasattr(var_value, "__len__"): 
                    try:
                        prodvars *= var_value[:]
                    except:
                        tmp = var_value[:]
                        prodvars *= tmp[np.newaxis,:]
                        del tmp
                else:
                    prodvars *= var_value
            else:
                prodvars *= var_value[arange]

        sumprod += prodvars

    # compute ensemble with unbiased estimator
    ensmbl = sumprod/(len(fopen_list)-1)

    return ensmbl #}}}

#@profile
def compute_autocorrelation_rlzn_ensemble(fopen_list, te):
    """
    Compute the autocorrelation from the direct formula
    over each of the realizations.

    Phillip Wolfram
    LANL
    09/23/2014
    """
    print 'Compute the autocorrelation'

    # initialize components of rho
    sumuu = 0.0
    sumvv = 0.0

    psiuu = 0.0
    psivv = 0.0

    sumup2 = 0.0
    sumvp2 = 0.0

    # get characteristics of mean velocity field
    fbs = netCDF4.Dataset('buoyancySurface.nc','r')
    lonCell = fbs.variables['lonCell']
    latCell = fbs.variables['latCell']
    lon = np.degrees(np.mod(lonCell[:]+np.pi,2*np.pi)-np.pi)
    lat = np.degrees(latCell[:])
    hull = spatial.ConvexHull(np.vstack((lon,lat)).T)                         
    triang = Triangulation(lon,lat)
    buoy_surf_zonal = fbs.variables['buoyancySurfaceVelocityZonal']
    buoy_surf_merid = fbs.variables['buoyancySurfaceVelocityMeridional']

 
    # build up layers for interpolation of particle layers
    interp_zonal = []
    interp_merid = []
    nlayers = len(fbs.dimensions['nBuoyancySurfaces'])
    for alayer in np.arange(nlayers):
        interp_zonal.append(LinearTriInterpolator(triang, buoy_surf_zonal[0,:,alayer]))
        interp_merid.append(LinearTriInterpolator(triang, buoy_surf_merid[0,:,alayer]))

    for num, afile in enumerate(fopen_list):
        print 'working on %d' % num
        # interpolate mean velocities onto points for the computation
        x = afile.variables['xParticle'][:te,:]
        y = afile.variables['yParticle'][:te,:]
        z = afile.variables['zParticle'][:te,:]
        latr, lonr = proj_lat_long(x,y,z)
        latr = np.degrees(latr)
        lonr = np.degrees(lonr)

        ubar = np.zeros(x.shape)
        vbar = np.zeros(x.shape)
        nparticle_layer = x.shape[1]/nlayers
        for alayer in np.arange(nlayers):
            ps = np.arange(alayer*nparticle_layer,(alayer+1)*nparticle_layer)
            ubar[:,ps] = interp_zonal[alayer](lonr[:,ps],latr[:,ps])
            vbar[:,ps] = interp_merid[alayer](lonr[:,ps],latr[:,ps])

        # compute portions of autocorrelation
        u = afile.variables['lonVel'][:te,:]
        up = u - ubar
        up0 = up[0,:]

        v = afile.variables['latVel'][:te,:]
        vp = v - vbar
        vp0 = vp[0,:]

        sumuu += up0*up
        sumvv += vp0*vp

        psiuu += up0*up0
        psivv += vp0*vp0
        
        sumup2 += np.nanmean(up**2.0, axis=0)
        sumvp2 += np.nanmean(vp**2.0, axis=0)
        

    fbs.close()

    # note division by psi removes need to divide the sums by the number of realizations
    sumuu /= psiuu 
    sumvv /= psivv

    sumup2 /= len(fopen_list)
    sumvp2 /= len(fopen_list)

    print 'done'

    return sumuu, sumvv, sumup2, sumvp2, lonr[0,:], latr[0,:], lon, lat, hull

def get_clusters(cluster_path): #{{{
    """
    Load the cluster information to link particles
    to clusters to compute dispersion.

    Phillip Wolfram, LANL
    09/19/2014
    """
    print 'loading cluster info'
    indicesToParticle = pickle.load(open(cluster_path+"/verticesToParticle.p","rb"))
    indicesOnCluster = pickle.load(open(cluster_path+"/verticesOnCell.p","rb"))
    maxIndices = pickle.load(open(cluster_path+"/maxVertices.p","rb"))
    print 'done'

    return indicesToParticle, indicesOnCluster, maxIndices #}}}

def compute_cluster_ensemble(var, indicesOnCluster, maxIndices, indicesToParticle): #{{{
    """
    Compute ensemble over particles in a cluster.

    Phillip Wolfram
    LANL
    09/19/2014
    """

    num_clusters = maxIndices.shape[0]
    if len(var.shape) == 1:
        meanvar = np.zeros((num_clusters,))
    elif len(var.shape) == 2:
        meanvar = np.zeros((var.shape[0],num_clusters))
    else:
        warnings.warn('did not have correct shape for ' + str(var) + ' with len(var.shape)='+ str(len(var.shape)))
        meanvar = None

    for aCluster, maxInd in enumerate(maxIndices):
        # get particles in cluster
        particles = indicesToParticle[indicesOnCluster[aCluster,0:maxInd]]

        # compute mean depending upon size of array
        if len(var.shape) == 1:
            meanvar[aCluster] = np.mean(var[particles])
        if len(var.shape) == 2:
            meanvar[:,aCluster] = np.mean(var[:,particles], axis=1)

    return meanvar #}}}

def compute_ensemble(fopen_list, var_list, range_list, indicesOnCluster, maxIndices, indicesToParticle): #{{{
    """
    Compute ensemble over each realization and averaging particles over their associated clusters.

    Phillip Wolfram
    LANL
    09/19/2014
    """

    rlzn_ensmbl = compute_rlzn_ensemble(fopen_list, var_list, range_list)

    # average is over 2nd dimension (if it exists), since first is time
    #ensmbl = compute_cluster_ensemble(rlzn_ensmbl, indicesOnCluster, maxIndices, indicesToParticle)

    ensmbl = rlzn_ensmbl

    return ensmbl #}}}

#@profile
def compute_autocorrelation_and_timescale(rootdir, folder_prefix, cluster_path, te): #{{{
    """
    Compute autocorrelation function, noting that we have to perform ensembles
    accross realizations.
    
    Phillip Wolfram
    LANL
    09/19/2014
    """

    print 'compute_autocorrelation_and_timescale'

    ####################################################################################################
    # set up paths and clusters
    ####################################################################################################

    rlzn_path_list = get_realization_paths(rootdir, folder_prefix)

    fopen_list = open_netcdf_files(rlzn_path_list,'output*nc')

    indicesToParticle, indicesOnCluster, maxIndices  = get_clusters(cluster_path)

    # just eddy part
    rhouu, rhovv, up2, vp2, lonp, latp, lon, lat, hull = compute_autocorrelation_rlzn_ensemble(fopen_list, te)
    
    np.save('rhouu'+str(te),rhouu)
    np.save('rhovv'+str(te),rhovv)
    np.save('up'+str(te),np.sqrt(up2))
    np.save('vp'+str(te),np.sqrt(vp2))
    np.save('lonp'+str(te),lonp)
    np.save('latp'+str(te),latp)
    np.save('lon'+str(te),lon)
    np.save('lat'+str(te),lat)
    np.save('hullsimplicies'+str(te),hull.simplices)
    
    rhouu = compute_cluster_ensemble(rhouu, indicesOnCluster, maxIndices, indicesToParticle)
    rhovv = compute_cluster_ensemble(rhovv, indicesOnCluster, maxIndices, indicesToParticle)
    up2   = compute_cluster_ensemble(up2,   indicesOnCluster, maxIndices, indicesToParticle)
    vp2   = compute_cluster_ensemble(vp2,   indicesOnCluster, maxIndices, indicesToParticle)
    lonp = compute_cluster_ensemble(lonp, indicesOnCluster, maxIndices, indicesToParticle)
    latp = compute_cluster_ensemble(latp, indicesOnCluster, maxIndices, indicesToParticle)

    np.save('rhouu_cluster'+str(te),rhouu)
    np.save('rhovv_cluster'+str(te),rhovv)
    np.save('up_cluster'+str(te),np.sqrt(up2))
    np.save('vp_cluster'+str(te),np.sqrt(vp2))
    np.save('lonp_cluster'+str(te),lonp)
    np.save('latp_cluster'+str(te),latp)

    close_netcdf_files(fopen_list)

    print 'compute_autocorrelation_and_timescale done'
    return rhouu, rhovv, np.sqrt(up2), np.sqrt(up2), lonp, latp, lon, lat, hull.simplices #}}}

def Rfit(x, Td, taue):
    """from Lumpkin2002Lagrangian, citing Garraffo et al (2001)"""
    return np.cos(np.pi*x/(2.*Td)) * np.exp(-(x/taue)*(x/taue))

def T_L(Td, taue):
    """from Lumpkin2002Lagrangian, citing Garraffo et al (2001)"""
    return np.sqrt(np.pi)/2.0 * taue * np.exp(-(np.pi*taue/(4*Td))*(np.pi*taue/(4*Td)))

def autocorrelation_curve_fit(xdata,ydata): #{{{

    # get good guess for Td
    test = np.where(ydata[:-1]*ydata[1:] < 0)
    i0 = test[0]
    if i0 != []:
        x0, x1 = xdata[i0[0]:i0[0]+2]
        y0, y1 = ydata[i0[0]:i0[0]+2]
        guess = x0 - y0/(y1-y0)*(x1-x0)
    else:
        guess = 90.
    
    # remove nan and bad values that shouldn't be used for the fit
    okvals = np.where(np.isfinite(xdata + ydata))
    x = xdata[okvals]
    y = ydata[okvals]

    # curve fit
    popt, pcov = curve_fit(Rfit, x, y, p0 = [guess,10.])

    # parameter return
    Td = popt[0]
    taue = popt[1]
    return Td, taue #}}}

def autocorr_curve_fit_guess(xdata,ydata): #{{{
    # look for first zero crossing and get value
    test = np.where(ydata[:-1]*ydata[1:] < 0)
    i0 = test[0]
    if i0 != []:
        x0, x1 = xdata[i0[0]:i0[0]+2]
        y0, y1 = ydata[i0[0]:i0[0]+2]
        Td = x0 - y0/(y1-y0)*(x1-x0)
    else:
        Td = 90.
    # least squares on remainder of function
    def Rfit_tau(x, taue):
        return Rfit(x,Td,taue)
    try:
        taue, pcov = curve_fit(Rfit_tau, xdata, ydata, p0=[7.])
    except:
        print 'trying simpler method, guess taue = 10.0 ...'
        taue = 10.0
    
    return Td, taue #}}}

def autocorr_curve_fit(xdata,ydata, p0=None): #{{{
    try:
        # remove nans
        okvals = np.where(np.isfinite(xdata + ydata))
        x = xdata[okvals]
        y = ydata[okvals]
        
        # get guess
        Td, taue = autocorr_curve_fit_guess(x,y)
        #print 'Td_guess = %f taue_guess = %f T_L_guess = %f' % (Td, taue, T_L(Td, taue))

        # get final solution
        popt, pcov = curve_fit(Rfit, x, y, p0=[Td, taue])

        Td = popt[0]
        taue = popt[1]
        
        #print 'Td= %f taue= %f T_L= %f' % (Td, taue, T_L(Td, taue))
        # make sure it is valid
        TL = np.abs(T_L(Td,taue))
    except:
        print 'failed with p0 = ', [Td, taue]
        Td = np.nan
        taue = np.nan
    return Td, taue  #}}}

def plot_save_fig(veleddy, autocorr, titlename, savename, lonp, latp, lon, lat, hulls, samplefreq=2., circ_radius=0, saveme=True): #{{{

    # factor of 2 for sampling frequency!
    x = samplefreq*np.arange(autocorr.shape[0])
  
    failed = False
    mask = False

    try:
        #Td, taue = autocorrelation_curve_fit(x, autocorr)
        Td, taue = autocorr_curve_fit(x, autocorr)

        #Td, taue = autocorrelation_curve_fit(x, autocorr)
        TL = np.abs(T_L(Td, taue))
        kappa = veleddy*veleddy*TL*86400
        
        fiterror = np.nanmean((Rfit(x,Td,taue) - autocorr)**2.0)

        if  (TL > 100. or fiterror > 0.03):
            mask = True
            #TL = np.nan
            #kappa = np.nan
            # silenced, too much output
            #print 'savename %s Td= %f taue= %f T_L= %f fiterror=%f' % (savename, Td, taue, TL, fiterror)
            #failed = True

        if np.isnan(Td) or np.isnan(taue): # or fiterror > 0.02:
          mask = True
          failed = True

    except:
        # the curve_fit failed, just use a nan to designate this
        finitepoints = np.sum(np.isfinite(autocorr))
        print 'warning: curve fit failed! nfinite points = ', finitepoints
        TL = np.nan
        kappa = np.nan
        Td = np.nan
        taue = np.nan
        mask = True
        if finitepoints > 0:
          failed = True
  
    # save data cases for failed runs that can be analyzed
    if failed:
        savename = savename.strip('.png') + '_failed.png'
        np.savez(savename.strip('.png'), x=x, y=autocorr)

    if saveme or failed:
        plt.figure()
        plt.plot(x, autocorr,'o')
        plt.plot(x, np.zeros(autocorr.shape),'k--')
        plt.ylabel(titlename)
        plt.xlabel(r"days")
        plt.ylim([-1.0, 1.0])
        plt.hold(True)
        xfine = np.linspace(0,samplefreq*autocorr.shape[0],100)
        plt.plot(xfine, Rfit(xfine,Td,taue),'k-')

        labelstr = r'$T_L=%.1f;~u=%.2f;~\kappa=%1.f$' % \
                (TL ,veleddy, kappa)
        plt.text(7, 0.75, labelstr, fontsize=25) 
        
        # inset for lat long / points
        #a = pylab.axes([0.72,0.7,0.15,0.15],axisbg='w')
        a = pylab.axes([0.72,0.30,0.15,0.15],axisbg='w')
        for simplex in hulls:
            plt.plot(lon[simplex],lat[simplex],'k-',lw=2,zorder=998)
        if circ_radius == 0:
            plt.plot(lonp,latp, 'rx',ms=12)
        else:
            circ = plt.Circle((lonp,latp),circ_radius, color='r',zorder=999)
            a.add_patch(circ)
        pylab.setp(a, xticks=[], yticks=[])
        #print 'saving %s' % savename
        plt.savefig(savename)
        plt.close('all')

    return TL, kappa, mask, Td, taue #}}}

def plot_overview(savename, titlename, lonp, latp, lon, lat, hulls, circ_radius=0): #{{{
    plt.figure()
    for simplex in hulls:
        p = plt.plot(lon[simplex],lat[simplex],'k-',lw=2,zorder=998)
    a = plt.gca()
    for x,y in zip(lonp,latp):
        circ = plt.Circle((x,y),circ_radius, color='r',zorder=999)
        a.add_patch(circ)
    plt.title(titlename)
    plt.xlabel('Lon')
    plt.ylabel('Lat')
    #print 'saving %s' % savename
    plt.savefig(savename)
    
    return #}}}

def plot_scalar(savename, titlename, lonp, latp, lon, lat, hulls, scalar, radius=500.0, clim=None): #{{{
    plt.figure()
    for simplex in hulls:
        p = plt.plot(lon[simplex],lat[simplex],'k-',lw=2,zorder=998)
    a = plt.gca()
    plt.scatter(lonp,latp,s=radius, c=scalar)
    plt.title(titlename)
    plt.xlabel('Lon')
    plt.ylabel('Lat')
    plt.colorbar()
    if clim is not None:
      plt.clim(clim)
    #print 'saving %s' % savename
    plt.savefig(savename)
    
    return #}}}

def plot_scalar_filtered(savename, title, triang, scalar, nfilt=0, clim=None, *args, **kwargs):
    plt.close()
    scalar = triang.smooth_laplacian(scalar,ntimes=nfilt, mean=np.nanmean)
    triang.plot_scalar(scalar, *args, **kwargs)
    if clim is not None:
      plt.clim(clim)
    plt.colorbar()
    plt.title(title)
    plt.savefig(savename)
    print 'saving %s' % (savename)
    plt.close()

def plot_points(savename, title, x, y, lon, lat, hulls):
    plt.close()
    for simplex in hulls:
        p = plt.plot(lon[simplex],lat[simplex],'k-',lw=2,zorder=998)
    plt.plot(x,y, 'r.')
    plt.xlim([-20., 20.])
    plt.ylim([20., 50.])
    plt.title(title)
    plt.savefig(savename)
    print 'saving %s' % (savename)
    plt.close()

def signed_log10(x): #{{{
    return np.sign(x)*np.log10(np.abs(x)) #}}}

if __name__ == "__main__":
    from optparse import OptionParser

    # Get command line parameters #{{{
    parser = OptionParser()
    parser.add_option("-r", "--root", dest="root",
                      help="folder root holding analysis folders",
                      metavar="FILE")
    parser.add_option("-p", "--prefix", dest="prefix",
                      help="folder prefix for analysis",
                      metavar="FILE")
    parser.add_option("-c", "--load_cluser_info_path", dest="clusterinfopath",
                      help="location of files with cluster info", 
                      metavar="DIR")
    parser.add_option("-d", "--dirname", dest="dirname",
                      help="dirname",
                      metavar="DIR")
    parser.add_option("-e", "--te", dest="te",
                      help="endingtime",
                      metavar="INT")

    options, args = parser.parse_args()

    if not options.root:
        #parser.error("Root directory is a required input.")
        options.root = '.'
    if not options.prefix:
        options.prefix = "analyze_output"
    if not options.clusterinfopath:
        options.clusterinfopath = None
    if options.te:
        options.te = int(options.te)
    if not options.dirname:
      # make host dir
      dirname = './clusters_lagrangian'
    else:
      dirname = options.dirname
    #}}}

# uncomment to build up *.npy files needed to make plots
    rhouu, rhovv, up, vp, lonp, latp, lon, lat, hulls = \
            compute_autocorrelation_and_timescale(options.root,options.prefix, options.clusterinfopath, options.te)
    
    rhouu = np.load('rhouu_cluster'+str(options.te)+'.npy')
    rhovv = np.load('rhovv_cluster'+str(options.te)+'.npy')
    up    = np.load('up_cluster'+str(options.te)+'.npy')
    vp    = np.load('vp_cluster'+str(options.te)+'.npy')
    lonp  = np.load('lonp_cluster'+str(options.te)+'.npy')
    latp  = np.load('latp_cluster'+str(options.te)+'.npy')
    lon   = np.load('lon'+str(options.te)+'.npy')
    lat   = np.load('lat'+str(options.te)+'.npy')
    hulls = np.load('hullsimplicies'+str(options.te)+'.npy')
    
    Nlayers = 5
    npart_layers = lonp.shape[0]/Nlayers
    print '%s clusters in %s layers' % (npart_layers, Nlayers)

    ### test #{{{
    ##triang = GeneralTriangulation(lonp[:npart_layers],latp[:npart_layers])

    ##def plot_scalar_filtered(triang, scalar, nfilt=0):
    ##    scalar = triang.smooth_laplacian(scalar,ntimes=nfilt)
    ##    plt.figure()
    ##    triang.plot_scalar(scalar)
    ##    plt.colorbar()

    ##plot_scalar_filtered(triang, up[:npart_layers], 0)
    ##plot_scalar_filtered(triang, vp[:npart_layers], 0)
    ##plt.show()
    ###}}}
  

    if not os.path.isdir(dirname):
        os.mkdir(dirname)
    for layer in np.arange(Nlayers):
        if not os.path.isdir(dirname+'/layer'+str(layer)):
            os.mkdir(dirname+'/layer'+str(layer))
    
    TL = np.zeros((npart_layers,Nlayers))
    Td = np.zeros((npart_layers,Nlayers))
    taue = np.zeros((npart_layers,Nlayers))
    kappa = np.zeros((npart_layers,Nlayers))
    uprime = np.zeros((npart_layers,Nlayers))
    mask = np.zeros((npart_layers,Nlayers))
    savefigs = False
    if True: 
        for layer in np.arange(Nlayers):
            print 'filter variables for layer %d' % (layer)
            triang = GeneralTriangulation(lonp[layer*npart_layers:(layer+1)*npart_layers], latp[layer*npart_layers:(layer+1)*npart_layers]) 
            
            #filter = triang.smooth_laplacian
            #nfilt = 10
            
            filter = triang.smooth_running
            nfilt = 1
            radius=1.0 # worked fairly well
            #radius=0.5# worked fairly well

            #filter = triang.smooth_none

            rhouu[:,layer*npart_layers:(layer+1)*npart_layers] = \
                    filter(rhouu[:,layer*npart_layers:(layer+1)*npart_layers].T, ntimes=nfilt, radius=radius, mean=np.nanmean).T
            rhovv[:,layer*npart_layers:(layer+1)*npart_layers] = \
                    filter(rhovv[:,layer*npart_layers:(layer+1)*npart_layers].T, ntimes=nfilt, radius=radius, mean=np.nanmean).T
            up[layer*npart_layers:(layer+1)*npart_layers] = \
                    filter(up[layer*npart_layers:(layer+1)*npart_layers], ntimes=nfilt, radius=radius, mean=np.nanmean)
            vp[layer*npart_layers:(layer+1)*npart_layers] = \
                    filter(vp[layer*npart_layers:(layer+1)*npart_layers], ntimes=nfilt, radius=radius, mean=np.nanmean)
            print 'done'

        np.savez('filtered_single_particle_stats'+str(options.te),filtername=filter.__name__, \
                radius=radius, nfilt=nfilt, rhouu=rhouu, rhovv=rhovv, up=up, vp=vp)

    if True:
        # load values
        data = np.load('filtered_single_particle_stats'+str(options.te)+'.npz')
        rhouu = data['rhouu']
        rhovv = data['rhovv']
        up = data['up']
        vp = data['vp']

        for layer in np.arange(Nlayers):
            for acluster in np.arange(npart_layers):
                # random sample for output
                savefigs = False
                if (np.random.random()*1000 < 2):
                    print 'plotting figure for layer %s/%s cluster %s/%s' % (layer+1, Nlayers, acluster+1, npart_layers)
                    savefigs = True

                # compute variables
                index = layer * npart_layers + acluster
                TLuu, kappauu, maskuu, Tduu, taueuu  = plot_save_fig(up[index], rhouu[:,index], \
                        r"$\rho_{uu}$", dirname+'/layer'+str(layer)+'/rhouumean'+("%.6d" % (acluster))+'.png',\
                        lonp[index], latp[index], lon, lat, hulls, saveme=savefigs)
                TLvv, kappavv, maskvv, Tdvv, tauevv = plot_save_fig(vp[index], rhovv[:,index], \
                        r"$\rho_{uu}$", dirname+'/layer'+str(layer)+'/rhovvmean'+("%.6d" % (acluster))+'.png',\
                        lonp[index], latp[index], lon, lat, hulls, saveme=savefigs)

                TL[acluster,layer] = 0.5*(TLuu+TLvv)
                Td[acluster,layer] = np.maximum(Tduu,Tdvv)
                taue[acluster,layer] = np.maximum(taueuu,tauevv)
                uprime[acluster, layer] = np.sqrt(np.mean(up[index]**2.0 +vp[index]**2.0))
                kappa[acluster,layer] = TL[acluster,layer] * 86400 * uprime[acluster, layer] * uprime[acluster, layer]
                mask[acluster,layer] = maskuu or maskvv
               
                if savefigs:
                    plt.close('all')
                if acluster % 1000 == 0:
                    print 'layer %s/%s cluster %s/%s' % (layer+1, Nlayers, acluster+1, npart_layers)

        np.save('kappa'+str(options.te),kappa)
        np.save('TL'+str(options.te),TL)
        np.save('Td'+str(options.te),Td)
        np.save('taue'+str(options.te),taue)
        np.save('uprime'+str(options.te),uprime)
        np.save('mask'+str(options.te),mask)

    kappa = np.load('kappa'+str(options.te)+'.npy')
    TL = np.load('TL'+str(options.te)+'.npy')
    Td = np.load('Td'+str(options.te)+'.npy')
    taue = np.load('taue'+str(options.te)+'.npy')
    uprime = np.load('uprime'+str(options.te)+'.npy')
    mask = np.load('mask'+str(options.te)+'.npy')

    plt.close('all')
    clim_kappa  = [3, 5]
    clim_TL     = [0, 30]
    clim_Td     = [0, 30]
    clim_taue   = [0, 30]
    clim_uprime = [0, 0.75]
    clim_L_L= [0, 500.0]
    nfilt = 10
    # figure out which points are not nans
    
    for layer in np.arange(Nlayers):
        kappai = kappa[:,layer]
        TLi = TL[:,layer]
        Tdi = Td[:,layer]
        tauei = taue[:,layer]
        uprimei = uprime[:,layer]
        maski = mask[:,layer]
        lonpi = lonp[layer*npart_layers:(layer+1)*npart_layers]
        latpi = latp[layer*npart_layers:(layer+1)*npart_layers]

        # remove points with nans so we "interpolate" holes
        okpoints = np.isfinite(kappai)
        badpoints = np.isnan(kappai)
        
        lonpibad = lonpi[badpoints]
        latpibad = latpi[badpoints]
        plot_points(dirname + '/bad_points'+str(layer)+'.png', 'NaNs', lonpibad, latpibad, lon, lat, hulls) 
        
        lonpi = lonpi[okpoints]
        latpi = latpi[okpoints]
        kappai = kappai[okpoints]
        maski = maski[okpoints]
        TLi = TLi[okpoints]
        Tdi = Tdi[okpoints]
        tauei = tauei[okpoints]
        uprimei = uprimei[okpoints]

        triang = GeneralTriangulation(lonpi, latpi)

        plot_scalar_filtered(dirname + '/kappa'+str(layer)+'_normed.png', '$\kappa$', triang, signed_log10(kappai),nfilt, clim=clim_kappa)
        plot_scalar_filtered(dirname + '/TL'+str(layer)+'_normed.png', '$T_L$', triang, TLi,nfilt, clim=clim_TL)
        plot_scalar_filtered(dirname + '/Td'+str(layer)+'_normed.png', '$T_d$', triang, Tdi,nfilt, clim=clim_Td)
        plot_scalar_filtered(dirname + '/taue'+str(layer)+'_normed.png', r"$\tau_e$", triang, tauei,nfilt, clim=clim_taue)
        plot_scalar_filtered(dirname + '/uprime'+str(layer)+'_normed.png', "$u'$", triang, uprimei, nfilt, clim=clim_uprime)
        plot_scalar_filtered(dirname + '/L_L'+str(layer)+'_normed.png', "$L_L~\\textrm{(km)}$", triang, TLi*uprimei*86400/1000., nfilt, clim=clim_L_L)

        plot_scalar_filtered(dirname + '/kappa'+str(layer)+'.png', '$\kappa$', triang, signed_log10(kappai), nfilt)
        plot_scalar_filtered(dirname + '/TL'+str(layer)+'.png', '$T_L$', triang, TLi, nfilt)
        plot_scalar_filtered(dirname + '/Td'+str(layer)+'.png', '$T_d$', triang, Tdi, nfilt)
        plot_scalar_filtered(dirname + '/taue'+str(layer)+'.png', r"$\tau_e$", triang, tauei, nfilt)
        plot_scalar_filtered(dirname + '/uprime'+str(layer)+'.png', "$u'$", triang, uprimei, nfilt)
        plot_scalar_filtered(dirname + '/L_L'+str(layer)+'.png', "$L_L~\\textrm{(km)}$", triang, TLi*uprimei*86400/1000., nfilt)

        plot_scalar_filtered(dirname + '/mask'+str(layer)+'.png', 'Mask', triang, maski, cmap=plt.get_cmap('Greys')) 
