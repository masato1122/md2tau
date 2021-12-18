# -*- coding: utf-8 -*-
import numpy as np
from optparse import OptionParser
from ase.build import (nanotube, molecule)

def split_chirality(chirality):
    import re
    nm = re.split('[:,/]', chirality)
    n = int(nm[0])
    m = int(nm[1])
    return n, m

def main(options):
    nc, mc = split_chirality(options.chirality)
    cnt1 = nanotube(nc, mc, length=options.Mswnt)
    fullname = "C%d"%(options.full)
    full1 = molecule(fullname)
    
    print(cnt1.cell)
    print(full1.cell)

    
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("--chirality", dest="chirality", type="string",
            default="10,10", help="chirality (n,m)")
    parser.add_option("--fullerene", dest="full", type="int",
            default="60", help="type of fullerene")
    
    parser.add_option("--Mswnt", dest="Mswnt", type="int",
            default="4", help="# of units of swnt")
    parser.add_option("--Mfull", dest="Mfull", type="int",
            default="1", help="# of units of peapod")
    parser.add_option("--Nbuffer", dest="Nbuffer", type="int",
            default=0, help="# of units of SWNT for an adiabatic layer")
    parser.add_option("--Lxy", dest="Lxy", type="float",
            default=None, help="Lxy [A]")
    (options, args) = parser.parse_args()

    split_chirality(options.chirality)
    main(options)

