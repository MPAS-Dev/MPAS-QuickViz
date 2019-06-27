"""
Created Friday Dec 14 2018

@author: Phillip J. Wolfram

version: 1.0

"""

# libraries
import numpy as np
from scipy.ndimage import label
import numpy.ma as ma
import matplotlib.pyplot as plt


def ow_eddy_labeled_mask(ow, owmin=-0.2, mineddycells=0):
    """
    Computes labels corresponding to unique areas of owmin
    isosurfaces.

    Phillip J. Wolfram
    12/14/2018
    """


    # get discrete regions of eddy activity
    mask = (ow < owmin)
    mask, nfeatures = label(mask)

    # filter regions to ensure each is larger than
    # mineddycells
    if mineddycells > 0:
        for ii in 1 + np.arange(nfeatures):
            if np.sum(mask == ii) < mineddycells:
                # remove the feature
                mask[mask==ii] = 0

        mask, nfeatures = label(mask > 0)

    return mask, nfeatures


def get_min_id(data, featuremask):
    """ 
    Returns the id of the minimum point within the featuremask of the data

    Phillip J. Wolfram
    12/14/2018
    """
    # returns the id of minvalue in data, only considering 
    # where featuremask = True
    return np.unravel_index(np.array(ma.array(data.data, mask=~featuremask).argmin(), 
                                     dtype='i'), data.shape)


def neighlist(me, nx, ny):
    """
    Returns an indexable list of neighbors for use with region masks.

    Phillip J. Wolfram
    12/14/2018
    """
        
    # ensure correct input type
    me = np.array(me)
    
    # check the neighboors for next min
    up = me + np.array([0,1])
    down = me - np.array([0,1])
    right = me + np.array([1,0])
    left = me - np.array([1,0])
    
    # combine directions
    neighs = np.vstack((me, up, down, right, left))
    
    # handle boundaries
    neighs[:,0] = np.clip(neighs[:,0], 0, nx)
    neighs[:,1] = np.clip(neighs[:,1], 0, ny)
    
    # remove repeated cell ids
    neighs = np.unique(neighs, axis=0)
    
    # return in index format so these 
    # values can be used for assignments
    return (neighs[:,0], neighs[:,1])


def r2check(ow, da, minr2points=50, plot=False):
    """
    Verifies that area contained within isosurface meets r2 criteria.

    Phillip J. Wolfram
    12/14/2018
    """
    # sort the 1D-coordinate
    ids = np.argsort(ow)

    # remove masked values
    ids = ids[~ow[ids].mask]

    # compute 1D coordinate values to compute r2
    a = np.cumsum(da[ids])
    ow = ow[ids]

    if plot:
        plt.figure()
        plt.plot(ow, a, '.')
        print('r^2= {:f}'.format(np.corrcoef(ow, a)[0,1]**2.0))

    # always "succeed" if there aren't eoungh points to test
    if len(ow) < minr2points:
        return 1.0
    else:
        # return r2 value
        return np.corrcoef(ow, a)[0,1]**2.0


def find_eddy(feature, ow, da, minr2points=30, mineddycells=100, r2cond=0.9):
    """
    Returns a found eddyfeature provided it meets r2 criteria 
    based on a feature of a ow field with paired area da for a
    minimum r2 of r2cond requiring minr2points to make the comparison.

    Phillip J. Wolfram
    12/14/2018
    
    """

    # initialize empty eddy feature
    eddyfeature = feature*False
    foundeddy = False

    # find staring minimum seed
    owmin = get_min_id(ow, feature)
    eddyfeature[owmin] = True

    # we will check the found min point and build list from there
    # note that at some time in the code, here aren't enough points
    # to assess whether there is an eddy
    r2 = 1.0
    while r2 > r2cond:
        
        # get neighboors to found point
        neighs = neighlist(np.array(np.where(eddyfeature)).T, 
                eddyfeature.shape[0]-1, eddyfeature.shape[1]-1)
        
        # clean up mask for testing, flagging neighbors 
        # that aren't currently found min or out of feature
        eddytest = feature*False
        eddytest[neighs] = True
        eddytest[eddyfeature] = False
        eddytest[~feature] = False

        # find minimum of neighs and add to eddy feature
        owmin = get_min_id(ow, eddytest)
        eddyfeature[owmin] = True
        
        r2 = r2check(ow[eddyfeature], da[eddyfeature], minr2points=minr2points)
        
        # thus, this consistently adds the minimum of the neighbor to the feature
        # to build out increasing isosurfaces w and paired areas
        
        # if there aren't any additional eddies to test, return result
        if np.sum(eddytest) == 0:
            # print 'No neighbooring eddies that can be tested.'
            eddyfound = (r2 > r2cond and np.sum(eddyfeature) >= mineddycells) 
            return eddyfound, eddyfeature
    
    # roll back including point that broke accuracy of estimate
    eddyfeature[owmin] = False

    # if eddy is found record it
    r2 = r2check(ow[eddyfeature], da[eddyfeature], minr2points=-1)
    eddyfound = (r2 > r2cond and np.sum(eddyfeature) >= mineddycells) 
    return eddyfound, eddyfeature


def find_all_eddies(ow, da, owmin=-0.2, mineddycells=100, minr2points=30, r2cond=0.9):
    """
    Finds all the eddies in the dataset and places them in a mask where an integer
    labels the eddy number.

    inputs
    ------
    ow:             Okubo Weiss array
    da:             cell area array
    owmin:          minimum Okubo Weiss threshold
    mineddycells:   minimum number of cells to be considered an eddy
    minr2points:    minimum number of points to be used to assess r2 criterion
    r2cond:         minimum acceptable r2 to be considered an eddy

    outputs
    -------
    alleddies:      an integer array numbering eddy label corresponding to each
                    eddy found.  0 indicates no eddy was found.

    Phillip J. Wolfram
    12/14/2018
    """

    # get the overall eddy mask
    mask, nfeatures = ow_eddy_labeled_mask(ow, owmin, mineddycells)


    # loop over each of the large features to find the eddies
    neddies = 0
    alleddies = np.zeros_like(mask)
    for afeature in 1 + np.arange(nfeatures):
        # print('On feature {:d} of {:d}'.format(afeature, nfeatures))

        # iterate on feature until all possible eddies 
        # (taking minpoint cells at a time are found)
        feature = (mask == afeature)

        # make sure there are values to compute
        if ow.mask[feature].min():
            # print 'Cannot analyze feature for eddies because it is masked.'
            continue
        
        # find the eddies if there are values to compute
        while np.sum(feature) > 0:
            foundeddy, eddyfeature = find_eddy(feature, ow, da, 
                                               minr2points, mineddycells, r2cond)

            if foundeddy:
                # store eddy information
                neddies += 1
                # print('Found eddy {:d}'.format(neddies))
                alleddies[eddyfeature] = neddies

            # remove eddy information from feature
            feature[eddyfeature] = False
            # print('{:d} points left in feature'.format(feature.sum()))
            if np.sum(feature) < minr2points:
                # terminate check because there aren't enough points to assess
                feature = False

    return alleddies


def eddy_centers(alleddies, ow):
    """
    Find the eddy centers for an alleddy index array.

    inputs
    ------
    alleddies:      an integer array numbering eddy label corresponding to each
                    eddy found.  0 indicates no eddy was found.
    ow:             Okubo Weiss array

    outputs
    -------
    eddycentermask: mask showing location of eddy center where integer indicates number
                    corresponding to alleddies index label.

    Phillip J. Wolfram
    12/14/2018
    """

    eddycenter = np.zeros_like(alleddies)

    # for each eddy
    for aeddy in 1 + np.arange(alleddies.max()):
        # find minimum index for eddy
        owmin = get_min_id(ow, (aeddy == alleddies))
        # mark eddy center
        eddycenter[owmin] = aeddy

    return eddycenter

