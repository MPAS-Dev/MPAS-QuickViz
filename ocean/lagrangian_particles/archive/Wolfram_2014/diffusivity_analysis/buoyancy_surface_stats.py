#!/usr/bin/env python

# import libraries / packages
import numpy as np
import netCDF4
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from GeneralTriangulation import GeneralTriangulation as Triangulation
from convert_ParaView_xml_to_matplotlib_colormap import convert_ParaView_xml_to_matplotlib_colormap


# function definitions

def compute_stats(fname_in):  #{{{

  print 'loading data'
  f_in = netCDF4.Dataset(fname_in, 'r')

  lonCell = np.degrees(np.mod(f_in.variables['lonCell'][:]+np.pi,2*np.pi)-np.pi)
  latCell = np.degrees(f_in.variables['latCell'][:])

  velzonal = f_in.variables['buoyancySurfaceVelocityZonal'][:,:,:]
  velmerid = f_in.variables['buoyancySurfaceVelocityMeridional'][:,:,:]
  depth = f_in.variables['buoyancySurfaceDepth'][:,:,:]
  print 'done'

  print 'compute zonal means'
  velzonal_mean = np.mean(velzonal, axis=0)
  print 'done'
  print 'compute meridional means'
  velmerid_mean = np.mean(velmerid, axis=0)
  print 'done'
  print 'compute depth means'
  depth_mean = np.mean(depth, axis=0)
  print 'done'
  
  print 'compute zonal std'
  velzonal_std = np.std(velzonal, axis=0)
  print 'done'
  print 'compute meridional std'
  velmerid_std = np.std(velmerid, axis=0)
  print 'done'
  print 'compute depth std'
  depth_std = np.std(depth, axis=0)
  print 'done'


  np.savez('buoyancySurfaceStats.npz',lonCell=lonCell, latCell=latCell, \
      velzonal_mean=velzonal_mean, velmerid_mean=velmerid_mean, depth_mean=depth_mean, \
      velzonal_std=velzonal_std, velmerid_std=velmerid_std, depth_std=depth_std)

  f_in.close()

  return  #}}}

def plot_data():
  data = np.load('buoyancySurfaceStats.npz')

  lonCell = data['lonCell']
  latCell = data['latCell']
  nlayers = data['velzonal_mean'].shape[1]
  npart_layers = data['velzonal_mean'].shape[0]/nlayers

  def plotme(tri,name, alayer, normalize=True, scalar=None):
    if scalar is None:
        scalar = data[name][:,alayer]
    if normalize:
        tri.plot_scalar(scalar, cmap=convert_ParaView_xml_to_matplotlib_colormap('redblue.xml'))
        cmax = np.max(np.abs(scalar))
        plt.clim([-cmax, cmax])
    else:
        tri.plot_scalar(scalar)
    plt.title(name)
    plt.colorbar()
    savename = 'buoyancy_surface_stats_layer' + "%.3d_" % (alayer)  + name + '.png'
    plt.savefig(savename)
    print 'saved %s' % (savename)
    plt.close()
  def quiverplot(tri, uname, vname, alayer):
      plt.figure()
      u = data[uname][:,alayer]
      v = data[vname][:,alayer]
      #uc = Triangulation(tri.x,tri.y,u).coarsen()
      #vc = Triangulation(tri.x,tri.y,v).coarsen()
      #Q = plt.quiver(uc.x, uc.y, uc.scalar, vc.scalar, units='inches', pivot='mid')
      Q = plt.quiver(tri.x, tri.y, u, v, units='width', pivot='mid', alpha=0.75)
      qk = plt.quiverkey(Q, 0.5, 0.95, 0.5, r'$0.5 \frac{m}{s}$', labelpos='E', coordinates='figure')
      savename = 'buoyancy_surface_stats_layer' + "%.3d_" % (alayer)  + 'quiver' + uname.lstrip('velzonal') + '.pdf'
      plt.savefig(savename)
      print 'saved %s' % (savename)
      plt.close()

  for alayer in np.arange(nlayers):
    plt.figure()
    tri = Triangulation(lonCell, latCell)
    plotme(tri,'velzonal_mean', alayer)
    plotme(tri,'velmerid_mean', alayer)
    quiverplot(tri,'velzonal_mean', 'velmerid_mean', alayer)
    quiverplot(tri,'velzonal_std', 'velmerid_std', alayer)
    plotme(tri,'depth_mean', alayer, normalize=False)
    plotme(tri,'velzonal_std', alayer, normalize=False)
    plotme(tri,'velmerid_std', alayer, normalize=False)
    plotme(tri,'depth_std', alayer, normalize=False)
    plotme(tri,'eddy_velocity', alayer, normalize=False, scalar = np.sqrt(data['velzonal_std'][:,alayer]**2.0 + data['velmerid_std'][:,alayer]**2.0))

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="inputfilename",
                      help="file to open for appending \
                      particle data 'particle_' extension",
                      metavar="FILE")
    options, args = parser.parse_args()

    if not options.inputfilename:
        parser.error("Input filename is a required input.")

    compute_stats(options.inputfilename)
    plot_data()
