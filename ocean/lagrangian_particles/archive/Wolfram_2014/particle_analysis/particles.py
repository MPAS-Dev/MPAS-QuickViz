#!/usr/bin/env python

import os
from subprocess import call
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import MultiLineString, LineString, shape, MultiPoint, Polygon
import shapely.geometry
from descartes import PolygonPatch
from scipy.stats import chi2
import time


def get_buffer(x, y, std, nstds):
    # build up multipoints for buffer
    abuffer = None
    for ax_s, ax_e, ay_s, ay_e, astd in zip(x[:-1], x[1:],y[:-1],y[1:], std):
        aline = LineString([(ax_s, ay_s),(ax_e, ay_e)])
        newbuffer = aline.buffer(astd*nstds)
        if abuffer is not None:
            abuffer = abuffer.union(newbuffer)
        else:
            abuffer = newbuffer
    return abuffer

def ellipse_points(xc, yc, cov, N=20, volume=0.118): #{{{
    """
    return collection of points for ellipse
    volume=0.118 for 1-sigma.
    Adapted from 
    http://www.nhsilbert.net/source/2014/06/bivariate-normal-ellipse-plotting-in-python/

    """
   
    # compute the width, height, and rotation for the ellipse
    def eigsorted(cov):
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        return vals[order], vecs[:,order]

    vals, vecs = eigsorted(cov)

    theta = np.arctan2(*vecs[:,0][::-1])
    width = 0
    height = 0
    width, height = 2 * np.sqrt(chi2.ppf(volume,2)) * np.sqrt(vals)

    #if (width < 0) or (height < 0):
    #    raise ValueError('width = %f and height %f must be positive!' % (width, height))

    # generate ellipse points
    t = np.linspace(0,2*np.pi,N,endpoint=False)
    # generate points
    xp = width/2.0  * np.cos(t)
    yp = height/2.0 * np.sin(t)
    ## rotate
    pr = np.vstack((xp,yp))
    Rot = np.array([[np.cos(theta), - np.sin(theta)],[np.sin(theta), np.cos(theta)]])
    pr = np.dot(Rot,pr)
    # translate
    xp = pr[0,:] + xc
    yp = pr[1,:] + yc

    #xp[-1] = xp[0]
    #yp[-1] = yp[0]

    return xp, yp #}}}

def test_ellipse_points(sigxx,sigyy, sigmaxy): #{{{
    x,y = ellipse_points(2.0,4.0,np.array([[sigxx,sigmaxy],[sigmaxy,sigyy]]))
    width = max(x)-min(x)
    height = max(y)-min(y)
    print width, height, width/height
    plt.plot(x,y)
    cov = np.cov(np.vstack((x,y)))
    print cov
    xt, yt = ellipse_points(x.mean(), y.mean(), cov)
    plt.plot(xt,yt,'o')
    # note, 'o''s dont have to be on line, but shape must be "parallel"

    plt.axis('equal')
    plt.show()
    return #}}}

def get_cov_shape(x, y, cov, nstds): #{{{
    abuffer = None
    for num, (ax, ay), in enumerate(zip(x,y)):
        xe, ye = ellipse_points(ax,ay,cov[:,:,num])
        points = [tuple(row) for row in np.vstack((xe, ye)).T]
        newbuffer = Polygon(points)

        # make sure that we aren't invalid due to machine precision
        if newbuffer.is_empty or newbuffer.area < 2*np.finfo(np.float).eps:
            continue

    # artificial indent for comments
    #try:
        if abuffer is not None:
            abuffer = abuffer.union(newbuffer)
        else:
            abuffer = newbuffer
    #except:
    #    ax = plt.gca()
    #    from shapely.validation import explain_validity
    #    for ab,cl in zip(abuffer,['r','g']):
    #        print explain_validity(ab)
    #        print ab
    #        import pdb; pdb.set_trace()
    #        x,y = ab.xy
    #        plot(x,y,'k-')
    #        ax.add_patch(PolygonPatch(ab,fc=cl,alpha=0.5))
    #    ax.add_patch(PolygonPatch(newbuffer,fc='b',alpha=0.5))
    #    plt.show()

    #    import pdb; pdb.set_trace()

    return abuffer #}}}


def get_shapely_hull(x,y):
    return MultiPoint([tuple(row) for row in np.vstack((x,y)).T]).convex_hull


class Particles:

    def __init__(self, file_name):
        data = np.load(file_name)
        self.lon = data['lon']
        self.lat = data['lat']
        self.Nt = self.lon.shape[0]
        self.Np = self.lon.shape[1]
        # hard coded 5 layrs for now!!!
        self.Npl = self.Np/5
        self.Nr = self.lon.shape[2]
        self.lonb = None
        self.latb = None
        self.hulls = None
        self.inittime = time.strftime("%c").replace(' ','_').replace(':','_')
        self.savefolder = './quick_save_figs/' 
        # make sure the save folder exists
        if not os.path.isdir(self.savefolder):
            os.mkdir(self.savefolder)

    def load_hull(self):
        self.lonb   = np.load('lon.npy')             
        self.latb   = np.load('lat.npy')
        self.hulls = np.load('hullsimplicies.npy')
        self.hulldef = {"type": "Polygon", "coordinates": [zip(self.lonb, self.latb)]}
        self.hullarea = shape(self.hulldef).area

    def load_hull_fine(self):
        self.lonb   = np.load('lon-1.npy')             
        self.latb   = np.load('lat-1.npy')
        self.hulls = np.load('hullsimplicies-1.npy')

    def get_particle_data(self,particle_num):
        """ Return particle data [Nt,Nr]
            
            Nt is the number of time
            snapshots and Nr is the 
            number of realziations.

            Phillip Wolfram
            10/02/2014
            LANL
        """
        return self.lon[:,particle_num,:], self.lat[:,particle_num,:]

    def mean_path(self,particle_num):
        #print (self.lon[:,particle_num,:].mean(axis=1)).shape
        return self.lon[:,particle_num,:].mean(axis=1), self.lat[:,particle_num,:].mean(axis=1)

    def plot_CI_asymmetric(self, particle_num, nstds=1): #{{{
        # get data
        lonm, latm = self.mean_path(particle_num)
        cov = self.cov_times(particle_num)
        abuffer = get_cov_shape(lonm, latm, cov, nstds)
       
        # plot CI for 1-std
        ax = plt.gca()
        #import pdb; pdb.set_trace()
        #try:
        if type(abuffer) is shapely.geometry.polygon.Polygon:
            patch = PolygonPatch(abuffer, fc='#a0a0a0', ec='#a0a0a0', alpha=0.5, zorder=997)
            ax.add_patch(patch)
        #except:
        else:
            for ab in abuffer:
                patch = PolygonPatch(ab, fc='#a0a0a0', ec='#a0a0a0', alpha=0.5, zorder=997)
                ax.add_patch(patch)

        # plot mean
        plt.plot(lonm,latm, c='k', ls='-', lw=1, zorder=998)#, alpha=0.5)
        plt.arrow(lonm[-2],latm[-2],np.diff(lonm[-2:])[0],np.diff(latm[-2:])[0],zorder=998, width=0.010, color='k')#, alpha=0.5 )
        #plt.plot(lonm[0],latm[0], 'bx', zorder=997)
        return #}}}

    def plot_CI_symmetric(self, particle_num, nstds=1): #{{{
        # get data
        lonm, latm = self.mean_path(particle_num)
        std = np.std(np.sqrt((self.lon[:,particle_num,:]-lonm[:,np.newaxis])**2.0 + (self.lat[:,particle_num,:]-latm[:,np.newaxis])**2.0), axis=1)
        abuffer = get_buffer(lonm, latm, std, nstds)
       
        # plot CI for 1-std
        patch = PolygonPatch(abuffer, fc='w', ec='k', alpha=0.5, zorder=997)
        ax = plt.gca()
        ax.add_patch(patch)
        # plot mean
        plt.plot(lonm,latm, c='k', ls='-', lw=1, zorder=997)
        plt.plot(lonm[0],latm[0], 'bx', zorder=997)
        return #}}}


    def plot_particle_paths(self,particle_num,labelme=False): #{{{
        for rlzn in np.arange(self.Nr):
            plt.plot(self.lon[:,particle_num,:], self.lat[:,particle_num,:],label='rlzn = '+str(rlzn))
        
        self.plot_CI(particle_num)

        plt.plot(self.lon[0,particle_num,:], self.lat[0,particle_num,:], 'x', c='#f9f9f9', ms=10) 
        if labelme:
            t = plt.text(self.lon[0,particle_num,0], self.lat[0,particle_num,0], s=str(particle_num), color='k',zorder=999)
            t.set_bbox(dict(color='k', facecolor='white', alpha=0.75,zorder=998))
        return

    def plot_boundary(self): #{{{
        if self.lonb is None:
            self.load_hull()
            self.load_hull_fine()
        for simplex in self.hulls:
            plt.plot(self.lonb[simplex], self.latb[simplex],'k-',lw=2)
        plt.hold(True)
        return  #}}}

    def cov_times(self, aparticle):
        cov =  np.zeros((2,2,self.Nt))
        for at in np.arange(self.Nt):
            cov[:,:,at] = np.cov(np.vstack((self.lon[at,aparticle,:],self.lat[at,aparticle,:])))
            #plt.figure()
            #plt.plot(self.lon[at,aparticle,:].ravel(), self.lat[at,aparticle,:].ravel(),'.')
            #print cov[:,:,at]
            #plt.show()
        return cov


    def compute_cov_area(self, particlelist, nstds=1.0): #{{{

        arealist = np.zeros((len(particlelist),))
        for num, aparticle in enumerate(particlelist):
            lonm, latm = self.mean_path(aparticle)
            cov = self.cov_times(aparticle)
            abuffer = get_cov_shape(lonm, latm, cov, nstds)
            arealist[num] = abuffer.area

        # save particle list
        np.save(self.savefolder+'particle_list'+self.inittime, particlelist)
        np.save(self.savefolder+'area_list'+self.inittime, arealist)

        return #}}}

    def plot_filled_space_CIs_asymmetric(self, particlelayer, nstds, frachullarealist, printprogress=True): #{{{
        particlelist = []
        arealist = []
        area = 0.0
        afrachullarea = frachullarealist.pop()

        # get first buffer
        aparticle = particlelayer[0]
        lonm, latm = self.mean_path(aparticle)
        cov = self.cov_times(aparticle)
        filledspace = get_cov_shape(lonm, latm, cov, nstds)
        particlelist.append(aparticle)
        arealist.append(filledspace.area)
      
        # loop over others
        for aparticle in particlelayer[1:]:
            lonm, latm = self.mean_path(aparticle)
            cov = self.cov_times(aparticle)
            abuffer = get_cov_shape(lonm, latm, cov, nstds)
            
            if filledspace.intersection(abuffer).is_empty:
                # otherwise add it to the list
                filledspace = filledspace.union(abuffer)
                
                # and compute total area
                area += abuffer.area
                #if printprogress:
                #    print '%% area = %.1f at particle %d' % (area/self.hullarea*100, aparticle - min(particlelayer))
                print '%% area = %.1f at particle %d' % (area/self.hullarea*100, aparticle - min(particlelayer))

                # note we are using it 
                particlelist.append(aparticle)
                arealist.append(abuffer.area)

                # plot it
                #self.plot_CI_asymmetric(aparticle, nstds)

            # break 
            if area > afrachullarea*self.hullarea:
                #plt.title(str(afrachullarea))
                #self.show()
                # get next checkpoint
                if len(frachullarealist):
                    afrachullarea = frachullarealist.pop()
                else:
                    # no more checkpoints
                    break

        # save particle list
        particlelistarray = np.asarray(particlelist)
        arealistarray = np.asarray(arealist)
        np.save(self.savefolder+'particle_list'+self.inittime, particlelistarray)
        np.save(self.savefolder+'area_list'+self.inittime, arealistarray)

        for aparticle in particlelist:
            self.plot_CI_asymmetric(aparticle, nstds)
        self.show()
        #print 'Done finding CIs'
        ## make plots
        #for aparticle in particlelist:
        #    self.plot_CI_asymmetric(aparticle, nstds)

        return #}}}

    def plot_filled_space_CIs_symmetric(self, particlelayer, nstds, frachullarealist): #{{{
        particlelifst = []
        area = 0.0
        afrachullarea = frachullarealist.pop()

        # get first buffer
        aparticle = particlelayer[0]
        lonm, latm = self.mean_path(aparticle)
        std = np.std(np.sqrt((self.lon[:,aparticle,:]-lonm[:,np.newaxis])**2.0 + \
                (self.lat[:,aparticle,:]-latm[:,np.newaxis])**2.0), axis=1)
        filledspace = get_buffer(lonm, latm, std, nstds)
        particlelist.append(aparticle)
      
        # loop over others
        for aparticle in particlelayer[1:]:
            lonm, latm = self.mean_path(aparticle)
            std = np.std(np.sqrt((self.lon[:,aparticle,:]-lonm[:,np.newaxis])**2.0 +\
                    (self.lat[:,aparticle,:]-latm[:,np.newaxis])**2.0), axis=1)
            abuffer = get_buffer(lonm, latm, std, nstds)
            
            if filledspace.intersection(abuffer).is_empty:
                # otherwise add it to the list
                filledspace = filledspace.union(abuffer)
                
                # and compute total area
                area += abuffer.area
                print '%% area = %.1f at particle %d' % (area/self.hullarea*100, aparticle - min(particlelayer))

                # note we are using it 
                particlelist.append(aparticle)

                # plot it
                self.plot_CI_symmetric(aparticle, nstds)

            # break 
            if area > afrachullarea*self.hullarea:
                plt.title(str(afrachullarea))
                self.show()
                # get next checkpoint
                if len(frachullarealist):
                    afrachullarea = frachullarealist.pop()
                else:
                    # no more checkpoints
                    break

        return #}}}

    def plot_hulls(self, particle_num):
        hull = get_shapely_hull(self.lon[:,particle_num,:].ravel(), \
                self.lat[:,particle_num,:].ravel())
        patch = PolygonPatch(hull, fc='w', ec='k', alpha=0.5, zorder=997)
        ax = plt.gca()
        ax.add_patch(patch)
        return 
    
    def plot_filled_space_hulls(self, particlelayer, frachullarea): #{{{
        particlelist = []
        area = 0.0

        # get first buffer
        aparticle = particlelayer[0]
        filledspace = get_shapely_hull(self.lon[:,aparticle,:].ravel(),self.lat[:,aparticle,:].ravel())
        particlelist.append(aparticle)
      
        # loop over others
        for aparticle in particlelayer[1:]:
            abuffer = get_shapely_hull(self.lon[:,aparticle,:].ravel(), \
                    self.lat[:,aparticle,:].ravel())
            
            if filledspace.intersection(abuffer).is_empty:
                # otherwise add it to the list
                filledspace = filledspace.union(abuffer)
                
                # and compute total area
                area += abuffer.area
                print '%% area = %.1f at particle %d' % (area/self.hullarea*100, aparticle - particlelayer[0])

                # note we are using it 
                particlelist.append(aparticle)

            # break 
            if area > frachullarea*self.hullarea:
                break

        print 'Done finding hulls'
        # make plots
        for aparticle in particlelist:
            self.plot_particle_paths(aparticle)
            self.plot_hulls(aparticle)
            self.plot_CI(aparticle,1.0)
        print 'Done printing hulls'

        return #}}}

    def show(self, showme=False):
        name = self.savefolder + time.strftime("%c").replace(' ','_').replace(':','_') + '.png'
        plt.savefig(name)
        if showme:
            plt.show()

def git_commit():
    commitname = "'Automatic commit from particles.py at " + time.strftime("%c").replace(' ','_').replace(':','_') + " from " + os.getcwd() +"'"
    cmd = 'cd ~/Documents/MPAS-pwolfram_fork/; git commit -m ' + commitname + ' ~/Documents/MPAS-pwolfram_fork/src/core_ocean/analysis_members/particle_analysis/particles.py'
    print cmd
    call(cmd,shell=True)  


if __name__ == "__main__":

    git_commit()

    particles = Particles('all_rlzn_particle_data.npz')
    #particles.plot_boundary()
    #particles.plot_particle_paths(1)
    #particles.plot_CI_symmetric(1)
    #ignore = np.array([13, 15, 5, 16, 18, 20, 22, 23, 24, 10, 7, 28, 19, 8, 17, 26, 30])
    #for i in np.arange(35):
    #    if i not in ignore:
    #        particles.plot_CI(i)
    layernum = 1 
    particle_range = np.arange(layernum*particles.Npl, (layernum+1)*particles.Npl)
    np.random.shuffle(particle_range)
    print particle_range
    #particles.plot_filled_space_CIs_symmetric(particle_range, 2.0, 0.10)
    #particles.plot_filled_space_hulls(particle_range, 0.12)
    perclist = list(np.arange(0.08,0.12,0.01))
    perclist.reverse()
    np.save(particles.savefolder+'particle_range'+particles.inittime, particle_range)
    plt.close('all')
    particles.plot_boundary()
    # note, bug here, nstds is not used but volume = 0.118 is used instead
    particles.plot_filled_space_CIs_asymmetric(particle_range, 2.0, perclist)
    

