#!/usr/bin/env python

from aggregate_diffusivity import write_paraview

if __name__ == "__main__":
    from optparse import OptionParser

    # Get command line parameters #{{{
    parser = OptionParser()
    parser.add_option("-e", "--ensemble", dest="ensemble",
                      help="number of ensemble for output",
                      metavar="INTEGER")

    options, args = parser.parse_args()

    if not options.ensemble:
        parser.error("Input ensemble number is a required input.")
    #}}}

    # write the files needed to compute diffusivity
    write_paraview('./ensembles/', options.ensemble,'ensemble_of_'+str(options.ensemble),15,30,2)

