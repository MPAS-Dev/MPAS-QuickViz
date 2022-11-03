#!/usr/bin/env python
# Phillip J. Wolfram
# 01/28/2016

import numpy as np
import matplotlib.pyplot as plt

def lines(plot=True):
    x = np.linspace(0,10,20)
    l1 = np.log(x)
    l2 = np.ones_like(x)*1.5

    if plot:
        plt.figure()
        plt.hold(True)
        plt.plot(x,l1,'.-', label='l1')
        plt.plot(x,l2,'.-', label='l2')
        #plt.show()

    return x, l1, l2

def first_line_intersection(x, l1, l2):
    """
    Finds the first intersection of two segmental lines with linear
    interpolation in segments
    Phillip J. Wolfram
    01/28/2016
    """

    y = l2 - l1

    # need to get starting index of first crossing
    shiftpoint = np.where(np.abs(np.diff(y > 0)))[0]

    if len(shiftpoint) > 1:
      print 'Intersection code not capable of finding multiple intersections, assuming first is the true intersection'
      shiftpoint = shiftpoint[0]

    x0 = x[shiftpoint]
    x1 = x[shiftpoint+1]
    y0 = y[shiftpoint]
    y1 = y[shiftpoint+1]

    # find x point of intersection
    xi = x0 - y0*(x1-x0)/(y1-y0)

    # find y point of intersection
    l10 = l1[shiftpoint]
    l11 = l1[shiftpoint+1]
    yi = (xi - x0)*(l11-l10)/(x1-x0) + l10

    # could be a test to make sure result is the same
    #l20 = l1[shiftpoint]
    #l21 = l1[shiftpoint+1]
    #yi2 = x0 + l20*(x1-x0)/(l21-l20)

    # handle nan cases
    if xi.size == 0:
      xi = np.array([np.nan])
    if yi.size == 0:
      yi = np.array([np.nan])

    return xi, yi

def test(): #{{{

    x, l1, l2 = lines(plot=True)
    xi, yi = first_line_intersection(x, l1, l2)

    plt.plot(xi, yi,'r.',label='intersection')
    plt.legend(loc='best')
    plt.show()
    return #}}}

if __name__ == "__main__":
    test()

