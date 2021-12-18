# -*- coding: utf-8 -*-
import os, sys
import os.path
import math
import numpy as np
from optparse import OptionParser

def file_check(file):
    if os.path.exists(file) == False:
        print("Cannot find %s"%file)
        sys.exit()
    return 0

class crystal:
    def __init__(self):
        self.Na = ""
        self.cube = np.zeros(6)
        self.el = ""
        self.coord = ""
        self.type = ""
        self.cent = np.zeros(3)

    def output(self):
        print("%s "%(self.el))
        for j in range(3):
            print("%f "%(self.coord[j]), end="")
        print("")
        return 0
    
    def make_xyz(self, OFILE):
        ofs = open(OFILE, "w")
        ofs.write("%d\n"%(self.Na))
        for i in range(6):
            ofs.write("%15.8f "%(self.cube[i]))
        ofs.write("\n")
        for ia in range(self.Na):
            ofs.write("%s "%(self.el[ia]))
            for j in range(3):
                ofs.write("%15.8e "%(self.coord[ia,j]))
            ofs.write("%d "%(self.type[ia]))
            ofs.write("\n")
        ofs.close()
        print("Output:", OFILE)
        return 0

    def readxyz(self, file):
        if os.path.exists(file) == False:
            print("Cannot find %s"%file)
            sys.exit()

        ifs = open(file)
        line = ifs.readline().split()
        self.Na = int(line[0])
        self.cube = np.zeros(6)
        line = ifs.readline().split()
        for i in range(6):
            self.cube[i] = float(line[i])

        self.el = []
        self.type = []
        self.coord = np.zeros(self.Na*3).reshape(self.Na, 3)
        for i in range(self.Na):
            line = ifs.readline().split()
            Nlen = len(line)
            self.el.append( line[0] )
            for j in range(3):
                self.coord[i,j] = float(line[1+j])
            if len(line) >= 5:
                self.type.append( int(line[4]) )
            else:
                self.type.append( 1 )
        return 0
    
    def cal_center(self):
        self.cent = np.zeros(3)
        for ia in range(self.Na):
            self.cent[:] += self.coord[ia,:] / float(self.Na)
        return 0
    
    def offset_xy(self, Lxy):
        self.cal_center()
        trans = np.zeros(3)
        trans[0] = Lxy * 0.5 - self.cent[0]
        trans[1] = Lxy * 0.5 - self.cent[1]
        trans[2] = self.cube[4]
        self.cube[0] = 0.; self.cube[1] = Lxy
        self.cube[2] = 0.; self.cube[3] = Lxy
        self.cube[4] += trans[2]
        self.cube[5] += trans[2]
        for ia in range(self.Na):
            self.coord[ia,:] += trans[:]
        return 0
        
    def make_long_swnt(self, sw0, Mswnt, Nbuffer):
        size = np.zeros(3)
        for i in range(3):
            size[i] = sw0.cube[2*i+1] - sw0.cube[2*i]
        
        Mall = Mswnt + 2*Nbuffer
        self.Na   = sw0.Na * Mall
        self.cube = np.zeros(6)
        self.cube[:] = sw0.cube[:]
        self.cube[5] = sw0.cube[4] + (sw0.cube[5] - sw0.cube[4]) * float(Mall)
        self.el = []; self.type = []
        self.coord = np.zeros(self.Na*3).reshape(self.Na,3)
        for iu in range(Mall):
            if iu < Nbuffer or iu > Mall-1-Nbuffer:
                TYPE = 2
                el = "Sn"
            else:
                TYPE = 1
                el = "C"
            trans = np.zeros(3)
            trans[2] = size[2] * float(iu)
            for ia in range(sw0.Na):
                self.el.append( el )
                self.type.append( TYPE )
                self.coord[ sw0.Na*iu + ia, : ] = sw0.coord[ia,:] + trans[:]

    def make_peapod(self, sw0, Mswnt, Nbuffer, full, Mfull):
        size = np.zeros(3)
        for i in range(3):
            size[i] = sw0.cube[2*i+1] - sw0.cube[2*i]
        
        Mall = Mswnt + 2*Nbuffer
        self.Na = sw0.Na * Mall + full.Na * Mfull
        self.cube = np.zeros(6)
        self.cube = sw0.cube
        self.cube[5] = sw0.cube[4] + size[2] * float(Mall)
        self.el = []; self.type = []
        self.coord = np.zeros(self.Na*3).reshape(self.Na,3)
        #--- 1. for outer swnt
        count = 0
        for iu in range(Mall):
            #'''
            #--- ver.1
            if iu < Nbuffer or iu > Mall-1-Nbuffer:
                TYPE = 4
                el = "Sn"
            else:
                TYPE = 1
                el = "C"
            #'''
            #--- ver.2
            #TYPE = 1
            #el = "C"
            
            trans = np.zeros(3)
            trans[2] = size[2] * float(iu)
            for ia in range(sw0.Na):
                self.el.append( el )
                self.type.append( TYPE )
                self.coord[ count, : ] = sw0.coord[ia,:] + trans[:]
                count += 1
        #--- 2. for peapod
        z0 = self.cube[4] + size[2] * float(Nbuffer)
        z1 = self.cube[4] + size[2] * float(Nbuffer + Mswnt)
        Lff = (z1 - z0) / float(Mfull)
        zfull = np.linspace(z0 + Lff*0.5, z1 - Lff*0.5, Mfull)
        for ifull in range(Mfull):
            trans = np.zeros(3)
            trans[2] = zfull[ifull] - full.cent[2]
            TYPE = 2 + ifull%2
            '''
            #--- ver.1
            if ifull%2 == 0:
                el = "O"
            else:
                el = "Si"
            '''
            #--- ver.2
            el = "C"
            for ia in range(full.Na):
                self.el.append( el )
                self.type.append( TYPE )
                self.coord[count, :] = full.coord[ia, :] + trans[:]
                count += 1
        return 0

def main(options):
    
    file_check(options.FILE1)
    file_check(options.FILE2)
        
    sw0  = crystal()
    sw0.readxyz(options.FILE1)
    full = crystal()
    full.readxyz(options.FILE2)
    
    sw0.offset_xy(options.Lxy)
    full.offset_xy(options.Lxy)
    
    #--- make SWNT
    sw_long = crystal()
    sw_long.make_long_swnt(sw0, options.Mswnt, options.Nbuffer)
    OFILE = "swnt.xyz"
    sw_long.make_xyz(OFILE)
    
    #--- make peapod
    pea = crystal()
    pea.make_peapod(sw0, options.Mswnt, options.Nbuffer, full, options.Mfull)
    OFILE = "pea.xyz"
    pea.make_xyz(OFILE)
     

#===== STRAT MAIN ========#
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("--Fswnt", dest="FILE1", type="string",\
        help="SWNT file")
    parser.add_option("--Ffull", dest="FILE2", type="string",\
        help="fullerene file")
    parser.add_option("--Mswnt", dest="Mswnt", type="int",\
        help="# of units of swnt")
    parser.add_option("--Mfull", dest="Mfull", type="int",\
        help="# of units of peapod")
    parser.add_option("--Nbuffer", dest="Nbuffer", type="int",\
        default=0, help="# of units in an adiabatic layer")
    parser.add_option("--Lxy", dest="Lxy", type="float",\
        help="Lxy [A]")
    (options, args) = parser.parse_args()
    if options.FILE1 == None or options.FILE2 == None:
        print("Input file name:")
        sys.exit()
    
    main(options)
  
