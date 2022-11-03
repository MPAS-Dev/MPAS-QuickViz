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
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation, LinearTriInterpolator
from scipy import spatial
from scipy.optimize import curve_fit
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
def compute_autocorrelation_rlzn_ensemble(fopen_list):
    """
    Compute the autocorrelation from the direct formula
    over each of the realizations.

    Phillip Wolfram
    LANL
    09/23/2014
    """

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

    for afile in fopen_list:
        # interpolate mean velocities onto points for the computation
        x = afile.variables['xParticle'][:,:]
        y = afile.variables['yParticle'][:,:]
        z = afile.variables['zParticle'][:,:]
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
        u = afile.variables['lonVel'][:,:]
        up = u - ubar
        up0 = up[0,:]

        v = afile.variables['latVel'][:,:]
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
def compute_autocorrelation_and_timescale(rootdir, folder_prefix, cluster_path): #{{{
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

    fopen_list = open_netcdf_files(rlzn_path_list,'buoyancySurface*nc')

    indicesToParticle, indicesOnCluster, maxIndices  = get_clusters(cluster_path)

    # just eddy part
    rhouu, rhovv, up2, vp2, lonp, latp, lon, lat, hull = compute_autocorrelation_rlzn_ensemble(fopen_list)
    

    np.save('rhouu',rhouu)
    np.save('rhovv',rhovv)
    np.save('up',np.sqrt(up2))
    np.save('vp',np.sqrt(vp2))
    np.save('lonp',lonp)
    np.save('latp',latp)
    np.save('lon',lon)
    np.save('lat',lat)
    np.save('hullsimplicies',hull.simplices)
    
    rhouu = compute_cluster_ensemble(rhouu, indicesOnCluster, maxIndices, indicesToParticle)
    rhovv = compute_cluster_ensemble(rhovv, indicesOnCluster, maxIndices, indicesToParticle)
    up2   = compute_cluster_ensemble(up2,   indicesOnCluster, maxIndices, indicesToParticle)
    vp2   = compute_cluster_ensemble(vp2,   indicesOnCluster, maxIndices, indicesToParticle)
    lonp = compute_cluster_ensemble(lonp, indicesOnCluster, maxIndices, indicesToParticle)
    latp = compute_cluster_ensemble(latp, indicesOnCluster, maxIndices, indicesToParticle)

    np.save('rhouu_cluster',rhouu)
    np.save('rhovv_cluster',rhovv)
    np.save('up_cluster',np.sqrt(up2))
    np.save('vp_cluster',np.sqrt(vp2))
    np.save('lonp_cluster',lonp)
    np.save('latp_cluster',latp)

    close_netcdf_files(fopen_list)

    print 'compute_autocorrelation_and_timescale done'
    return rhouu, rhovv, np.sqrt(up2), np.sqrt(up2), lonp, latp, lon, lat, hull.simplices #}}}

def Rfit(x, Td, taue):
    """from Lumpkin2002Lagrangian, citing Garraffo et al (2001)"""
    return np.cos(np.pi*x/(2*Td)) * np.exp(-(x/taue)**2.0)

def T_L(Td, taue):
    """from Lumpkin2002Lagrangian, citing Garraffo et al (2001)"""
    return np.sqrt(np.pi)/2.0 * taue * np.exp(-(np.pi*taue/(4*Td))**2.0)

def plot_save_fig(veleddy, autocorr, titlename, savename, lonp, latp, lon, lat, hulls, circ_radius=0): #{{{
    plt.figure()
    x = 2*np.arange(autocorr.shape[0])
    xfine = np.linspace(0,2*autocorr.shape[0],100)
    autocorrfine = np.interp(xfine,x,autocorr)
    plt.plot(x, autocorr,'o')
    plt.plot(x, np.zeros(autocorr.shape),'k--')
    plt.ylabel(titlename)
    plt.xlabel(r"days")
    plt.ylim([-1.0, 1.0])
    plt.hold(True)
    popt, pcov = curve_fit(Rfit, xfine, autocorrfine, p0=(7.0,30.0))
    #labelstr = r'$\tau_e=%.1f~T_D=%.1f~T_L=%.1f$' % \
    #        (popt[0],popt[1],T_L(popt[0],popt[1]))
    #labelstr = r'$T_L=%.1f;~\tau_e=%.1f;~T_D=%.1f$' % \
    #        (T_L(popt[0],popt[1]),popt[1],popt[0])
    TL = T_L(popt[0],popt[1])
    kappa = veleddy*veleddy*TL*86400
    labelstr = r'$T_L=%.1f;~u=%.2f;~\kappa=%1.f$' % \
            (TL ,veleddy, kappa)
    plt.text(7, 0.75, labelstr, fontsize=25) 
            
    plt.plot(xfine, Rfit(xfine,popt[0],popt[1]),'k-')
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
    return TL, kappa #}}}

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

    options, args = parser.parse_args()

    if not options.root:
        #parser.error("Root directory is a required input.")
        options.root = '.'
    if not options.prefix:
        options.prefix = "analyze_output"
    if not options.clusterinfopath:
        options.clusterinfopath = None
    #}}}

    #rhouu, rhovv, up, vp, lonp, latp, lon, lat, hulls = \
    #        compute_autocorrelation_and_timescale(options.root,options.prefix, options.clusterinfopath)
    
    rhouu = np.load('rhouu_cluster.npy')
    rhovv = np.load('rhovv_cluster.npy')
    up    = np.load('up_cluster.npy')
    vp    = np.load('vp_cluster.npy')
    lonp  = np.load('lonp_cluster.npy')
    latp  = np.load('latp_cluster.npy')
    lon   = np.load('lon.npy')
    lat   = np.load('lat.npy')
    hulls = np.load('hullsimplicies.npy')

# boxes #{{{
#    # define the boxes
#    width_array = np.array([0.5, 0.5, 1.0])
#    cx_array = np.array([-5.0, -10.0, 0.0])
#    cy_array = np.array([32.5, 28.0, 35.0])
#
#    # make host dir
#    if not os.path.isdir('cluster_boxes'):
#        os.mkdir('cluster_boxes')
#    
#    # loop over each box
#    Nlayers = 5
#    npart_layers = lonp.shape[0]/Nlayers
#    lonp_layer = lonp[:npart_layers]
#    latp_layer = latp[:npart_layers]
#    for width,cx,cy in zip(width_array,cx_array,cy_array):
#        box = np.where(((cx-width/2) < lonp_layer) & ((cx+width/2) > lonp_layer) & ((cy-width/2) < latp_layer) & ((cy+width/2) > latp_layer))
#        print 'computing auto-correlation for %d clusters in box' % box[0].shape
#        dirname = 'cluster_boxes/' + 'X'+str(cx)+'_Y'+str(cy)+'_W'+str(width)
#        if not os.path.isdir(dirname):
#            os.mkdir(dirname)
#        for entry in box[0]:
#            for layer in np.arange(0,Nlayers):
#                if not os.path.isdir(dirname+'/layer'+str(layer)):
#                    os.mkdir(dirname+'/layer'+str(layer))
#                index = layer * npart_layers + entry
#                plot_save_fig(rhouu[1:,index],r"$\rho_{uu}$", dirname+'/layer'+str(layer)+'/rhouu'+"{0:06d}".format(entry)+'.png',\
#                        lonp[index], latp[index], lon, lat, hulls)
#                plot_save_fig(rhouv[1:,index],r"$\rho_{uv}$", dirname+'/layer'+str(layer)+'/rhouv'+"{0:06d}".format(entry)+'.png',\
#                        lonp[index], latp[index], lon, lat, hulls)
#                plot_save_fig(rhovu[1:,index],r"$\rho_{vu}$", dirname+'/layer'+str(layer)+'/rhovu'+"{0:06d}".format(entry)+'.png',\
#                        lonp[index], latp[index], lon, lat, hulls)
#                plot_save_fig(rhovv[1:,index],r"$\rho_{vv}$", dirname+'/layer'+str(layer)+'/rhovv'+"{0:06d}".format(entry)+'.png',\
#                        lonp[index], latp[index], lon, lat, hulls)
#                plt.close('all')
#
#        # compute mean over box
#        for layer in np.arange(0,Nlayers):
#            index = layer * npart_layers + box[0]
#            plot_save_fig(rhouu[1:,index].mean(axis=1),r"$\rho_{uu}$", dirname+'/layer'+str(layer)+'/rhouumean.png',\
#                    lonp[index], latp[index], lon, lat, hulls)
#            plot_save_fig(rhouv[1:,index].mean(axis=1),r"$\rho_{uv}$", dirname+'/layer'+str(layer)+'/rhouvmean.png',\
#                    lonp[index], latp[index], lon, lat, hulls)
#            plot_save_fig(rhovu[1:,index].mean(axis=1),r"$\rho_{vu}$", dirname+'/layer'+str(layer)+'/rhovumean.png',\
#                    lonp[index], latp[index], lon, lat, hulls)
#            plot_save_fig(rhovv[1:,index].mean(axis=1),r"$\rho_{vv}$", dirname+'/layer'+str(layer)+'/rhovvmean.png',\
#                    lonp[index], latp[index], lon, lat, hulls)
#            plt.close('all')
#
#        # make vertical layer comparisons via bash script with imagemagick
#        call('./combine_rho_curves.sh '+ dirname, shell=True)
#}}}

    # make host dir
    if not os.path.isdir('cluster_balls'):
        os.mkdir('cluster_balls')
    
    # loop over each ball
    Nlayers = 5
    npart_layers = lonp.shape[0]/Nlayers
    lonp_layer = lonp[:npart_layers]
    latp_layer = latp[:npart_layers]
    # define the circles
    Ndimx = 8
    Ndimy = 8
    x,y = np.meshgrid(np.linspace(-15,15,Ndimx),np.linspace(20,50,Ndimy))
    dist = np.sqrt((13./16.*x)**2 + (y-35)**2.0)
    #indomain = np.where(dist < 13.8)
    indomain = np.where(dist < 13.5)
    cx_array = (x[indomain]).ravel()
    cy_array = (y[indomain]).ravel()
    print 'computing for %s different clusters' % (len(cx_array))
    radius   = 2.0
    plot_overview('cluster_balls/cluster_locations.png', str(len(cx_array)) + ' Cluster locations', cx_array, cy_array, lon, lat, hulls, circ_radius=radius)

    # get nearest neighbors within radius on grided domain
    points = zip(lonp_layer.ravel(), latp_layer.ravel())
    tree = spatial.KDTree(points)
    
    ballpoints = zip(cx_array.ravel(), cy_array.ravel())
    balls = tree.query_ball_point(ballpoints,radius)
  
    TL = np.zeros((len(ballpoints),Nlayers))
    kappa = np.zeros((len(ballpoints),Nlayers))
    uprime = np.zeros((len(ballpoints),Nlayers))
    num = 0
    for center, ball in zip(ballpoints,balls):
        if not ball:
            continue
        print 'computing mean auto-correlation for %d clusters in ball' % len(ball)
        dirname = 'cluster_balls/X%.2f_Y%.2f_R%.2f' % (center[0],center[1],radius)
        if not os.path.isdir(dirname):
            os.mkdir(dirname)

        # compute mean over box
        for layer in np.arange(0,Nlayers):
            if not os.path.isdir(dirname+'/layer'+str(layer)):
                os.mkdir(dirname+'/layer'+str(layer))
            index = layer * npart_layers + ball
            uprime[num, layer] = np.nanmean(np.sqrt(up[index]**2.0 +vp[index]**2.0))
            TLuu, kappauu = plot_save_fig(np.nanmean(up[index]), np.nanmean(rhouu[1:,index],axis=1),r"$\rho_{uu}$", dirname+'/layer'+str(layer)+'/rhouumean.png',\
                    center[0], center[1], lon, lat, hulls, circ_radius=radius)
            TLvv, kappavv = plot_save_fig(np.nanmean(vp[index]), np.nanmean(rhovv[1:,index],axis=1),r"$\rho_{vv}$", dirname+'/layer'+str(layer)+'/rhovvmean.png',\
                    center[0], center[1], lon, lat, hulls, circ_radius=radius)

            Uprime2 = np.nanmean(vp[index]**2.0 + up[index]**2.0)
            TL[num,layer] = np.max([TLuu, TLvv])
            kappa[num,layer] = TL[num,layer]*86400*Uprime2
            
            plt.close('all')
        num += 1

    #for  center, ball in (zip(ballpoints,balls)): #{{{
    #    if not ball:
    #        continue
    #    print 'computing auto-correlation for %d clusters in ball' % len(ball)
    #    dirname = 'cluster_balls/X%.2f_Y%.2f_R%.2f' % (center[0],center[1],radius) # plot individual points
    #    for entry in ball:
    #        for layer in np.arange(0,Nlayers):
    #            index = layer * npart_layers + entry
    #            plot_save_fig(rhouu[1:,index],r"$\rho_{uu}$", dirname+'/layer'+str(layer)+'/rhouu'+"{0:06d}".format(entry)+'.png',\
    #                    lonp[index], latp[index], lon, lat, hulls)
    #            plot_save_fig(rhouv[1:,index],r"$\rho_{uv}$", dirname+'/layer'+str(layer)+'/rhouv'+"{0:06d}".format(entry)+'.png',\
    #                    lonp[index], latp[index], lon, lat, hulls)
    #            plot_save_fig(rhovu[1:,index],r"$\rho_{vu}$", dirname+'/layer'+str(layer)+'/rhovu'+"{0:06d}".format(entry)+'.png',\
    #                    lonp[index], latp[index], lon, lat, hulls)
    #            plot_save_fig(rhovv[1:,index],r"$\rho_{vv}$", dirname+'/layer'+str(layer)+'/rhovv'+"{0:06d}".format(entry)+'.png',\
    #                    lonp[index], latp[index], lon, lat, hulls)
    #            plt.close('all')
    # }}}

    np.save('kappa',kappa)
    np.save('T_L',T_L)

    clim_kappa  = [0, 20000]
    clim_TL     = [0, 20]
    clim_uprime = [0, 0.40]
    for layer in np.arange(0,Nlayers):
      plot_scalar('cluster_balls/kappa'+str(layer)+'_normed.png', '$\kappa$', cx_array, cy_array, lon, lat, hulls, kappa[:,layer],clim=clim_kappa)
      plot_scalar('cluster_balls/TL'+str(layer)+'_normed.png', '$T_L$', cx_array, cy_array, lon, lat, hulls, TL[:,layer],clim=clim_TL)
      plot_scalar('cluster_balls/uprime'+str(layer)+'_normed.png', "$u'$", cx_array, cy_array, lon, lat, hulls, uprime[:,layer], clim=clim_uprime)
      
      plot_scalar('cluster_balls/kappa'+str(layer)+'.png', '$\kappa$', cx_array, cy_array, lon, lat, hulls, kappa[:,layer])
      plot_scalar('cluster_balls/TL'+str(layer)+'.png', '$T_L$', cx_array, cy_array, lon, lat, hulls, TL[:,layer])
      plot_scalar('cluster_balls/uprime'+str(layer)+'.png', "$u'$", cx_array, cy_array, lon, lat, hulls, uprime[:,layer])

    for center, ball in (zip(ballpoints,balls)):
        dirname = 'cluster_balls/X%.2f_Y%.2f_R%.2f' % (center[0],center[1],radius) 
        # make vertical layer comparisons via bash script with imagemagick
        call('./combine_rho_curves.sh '+ dirname, shell=True)

    call('./combine_rho_curves_total.sh cluster_balls', shell=True)
