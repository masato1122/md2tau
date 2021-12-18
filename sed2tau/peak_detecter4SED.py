#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os, sys
import os.path
import gc
import numpy as np
from optparse import OptionParser
import multiprocessing as mp
import detect_peaks
import tips
import read_file

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

#--- for multiprocess
def argwrapper(args):
    return args[0](*args[1:])

#///// Definition ////////
#SED_DIR = "SED_DATA"
#/////////////////////////

#-- directory for the log
#directory = SED_DIR    
#if os.path.exists(directory):
#    command = "rm -r ./%s; mkdir -p %s"%(directory, directory)
#else:

def make_sed_dir(directory):
    if os.path.exists(directory):
        command = ""
    else:
        command = "mkdir -p ./%s"%(directory)
        os.system(command)
    return 0

#--- check file
def file_check(file):
    if os.path.exists(file) == False:
        print "Cannot find %s"%file
        sys.exit()

#--- output peak file
'''
ALLDATA[nblock,ndat,nline]
PEAKXY[block][3][npeak(ib)]
'''
def output_peaks(OFILE, I_BLOCK, PEAK_IXY, IFREQ, ISED):
    ofs = open(OFILE, "w")
    nblock = len(I_BLOCK)
    ofs.write("#Label  ID   Frequency[THz]     SED[-]\n")
    
    for ib_peak in range(nblock):

        ib_all  = I_BLOCK[ib_peak]

        #--- 1.initial line for each block
        ofs.write("# Iblock: %d\n"%(ib_all))

        #--- 2.output peaks
        npeak = len(PEAK_IXY[ib_peak][0])
        for ip in range(npeak):
            num =  int(PEAK_IXY[ib_peak][0][ip])
            xp = float(PEAK_IXY[ib_peak][1][ip])
            yp = float(PEAK_IXY[ib_peak][2][ip])
            ofs.write("%3d  %5d  %15.7f  %15.7e "%(ip+1, num, xp, yp))
            ofs.write("\n")
        ofs.write("\n")

    print "Output:", OFILE
    return 0

def read_SED_block_file(SFILE, IFREQ, ISED):
    ifs = open(SFILE, "r")
    xx = []; yy = []; count = 0
    while ifs:
        line = ifs.readline()
        data = line.split()
        if len(data) == 0: break
        if line[0] == "#": continue
        xx.append(float(data[IFREQ]))
        yy.append(float(data[ISED]))
        count += 1
    ifs.close()
    xdat = np.zeros(count)
    ydat = np.zeros(count)
    for i in range(count):
        xdat[i] = xx[i]
        ydat[i] = yy[i]
    return xdat, ydat

def Detect_Peaks(iblock, LFLAG, pwidth, pmin, PMIN_FLAG, IFREQ, ISED, WUNIT, SED_DIR):

    OFILE = "./%s/block%d.txt"%(SED_DIR, iblock)
    xdat, ydat = read_SED_block_file(OFILE, IFREQ, ISED) 
    if WUNIT == "Hz":
        xdat *= 1e-12
        
    if PMIN_FLAG == 2:
        ytemp = np.zeros(len(ydat))
        ytemp[:] = ydat[:]
        ytemp.sort()
        iget = int(float(len(ydat)) * (1. - PEAK_MIN))
        pmin = ytemp[iget]
        del ytemp 
    
    print " iblock: %3d   peak min.SED: %6.3f"%(iblock, pmin)
    Ipeak = detect_peaks.ymin(xdat, ydat, pwidth, pmin)
    xpeak = xdat[Ipeak]
    ypeak = ydat[Ipeak]
    del xdat; gc.collect()
    del ydat; gc.collect()   
    return Ipeak, xpeak, ypeak

#-----
def divide_SED_data(sfile, ISED, SED_DIR):
    ifs = open(sfile, "r")
    nline = sum(1 for line in open(sfile))
    iblock = 0; FLAG = 0; LABEL = "#"
    sed_list = []
    for il in range(nline):           
        line = ifs.readline()
        data = line.split()
        if len(data) == 0:
            ofs.close()
            FLAG = 0
            continue
        if line[0] == "#":
            LABEL = line
            continue
        if FLAG == 0:
            ofile = "./%s/block%d.txt"%(SED_DIR, iblock)
            ofs = open(ofile, "w")
            ofs.write("%s"%(LABEL))
            iblock += 1
            nfreq = 0
        ofs.write("%s"%(line))
        sed_list.append(float(data[ISED]))
        FLAG = 1; nfreq += 1
    nblock = iblock
    return nblock, nfreq, sed_list


#--------- START MAIN ------------
if __name__ == "__main__":

    print "----------------------------------"
    print ""
    print "     peak_detecter ver.1.0"
    print ""
    print "----------------------------------"

    parser = OptionParser()

    #--- Input file name
    parser.add_option("-f", "--ifile",  dest="ifile", type="string", \
            help="Input file name: ")

    parser.add_option("-s", "--sedfile",  dest="sfile", type="string", \
            help="Input SED file name: ")
    parser.add_option("--seddir",  dest="sedir", type="string", \
            help="Directory of SED")
    parser.add_option("-n", "--nprocs",  dest="nprocs", type="int", \
            help="Number of parallelization: ")

    (options, args) = parser.parse_args()
    if options.ifile == None:
        print "Input input file name (-f, --ifile)."
        tips.mkinputfile("peak_example.in")
        sys.exit()
    IFILE = options.ifile
    SED_DIR = options.sedir
    file_check(IFILE)
   
    make_sed_dir(SED_DIR)

    #---- read input parameters
    LFLAG = 1
    #-- for SED file
    if options.sfile == None:
        SFILE = tips.read_input_file(IFILE, "SEDFILE", 1)
    else:
        SFILE = options.sfile
    data  = tips.read_input_file(IFILE, "IFREQ", 1);    IFREQ  = int(data) - 1
    data  = tips.read_input_file(IFILE, "ISED", 1);     ISED   = int(data) - 1
    WUNIT = tips.read_input_file(IFILE, "FREQUNIT", 1)
    #-- to make the calculation fast
    if options.nprocs == None:
        data  = tips.read_input_file(IFILE, "NPROCS", 0);   nprocs = int(data)
    else:
        nprocs = options.nprocs
    data  = tips.read_input_file(IFILE, "NKSKIP", 1);   NKSKIP = int(data)
    #-- peak info.
    data = tips.read_input_file(IFILE, "PEAKMIN_FLAG", 1); PMIN_FLAG  = int(data)
    data = tips.read_input_file(IFILE, "PEAKWIDTH", 1);    PEAK_WIDTH = float(data)
    data = tips.read_input_file(IFILE, "PEAKMIN", 1);      PEAK_MIN   = float(data)

    file_check(SFILE)
    
    #-- wavevector
    #Lunit = options.Lunit
    #kmax = 2. * np.pi / Lunit
    
    #--- 1. divide blocks
    nblock, nfreq, sed_list = divide_SED_data(SFILE, ISED, SED_DIR)
    nbused = int(nblock/NKSKIP)
    
    #--- print parameters
    print "Nblock(all)     : %5d"%(nblock)
    print "Nblock(used)    : %5d"%(nbused)
    print "Nkskip          : %5d"%(NKSKIP)
    print "Nfreq           : %5d"%(nfreq)
    print "Peak width [THz]: %7.4f"%(PEAK_WIDTH)

    if PMIN_FLAG == 0: 
        print "min.SED for peak: %7.4f"%(PEAK_MIN)
    
    elif PMIN_FLAG == 1:
        print "min.SED will be determined using top x% data of all data"
        print "                : %4.1f%%"%(100.*PEAK_MIN)
        
        #-- ver.2: cal min.SED for peak
        nall = len(sed_list); sed = np.zeros(nall)
        sed[:] = sed_list[:]
        sed.sort()
        iget = int(float(nall) * (1. - PEAK_MIN))
        PEAK_MIN = sed[iget]
    
    elif PMIN_FLAG == 2:
        print "min.SED will be determined using top x% data at each k"
        print "                : %4.1f%%"%(100.*PEAK_MIN)
        if PEAK_MIN < 0. or 1. < PEAK_MIN:
            print "Error"
            sys.exit()
    else:
        print "Error, PMIN_FLAG == %d"%(PMIN_FLAG)
        sys.exit()
 
    #--- 2. Fitting for each block: PEAK_IXY_ALL[nblock][2][npeaks(ib)]
    #- 2-1. for parallel
    if nprocs >= 2:
        print "\n===  Parallel calculation (%d)  ===\n"%(nprocs)
        pool = mp.Pool(nprocs)
        func_args = []
        for iblock in xrange(nbused):
            func_args.append((Detect_Peaks, iblock*NKSKIP, LFLAG, PEAK_WIDTH, PEAK_MIN, PMIN_FLAG, \
                    IFREQ, ISED, WUNIT, SED_DIR))
        PEAK_IXY_ALL = pool.map(argwrapper, func_args)
    
    #- 2-2. for serial
    else:
        print "\n===  Serial calculation  ===\n"
        PEAK_IXY_ALL = []
        for iblock in xrange(nbused):
            PEAK_IXY_ALL.append(Detect_Peaks(iblock*NKSKIP, LFLAG, PEAK_WIDTH, PEAK_MIN, PMIN_FLAG, \
                    IFREQ, ISED, WUNIT, SED_DIR))
    
    #--- 3-a. prepare
    I_block = []
    for ib in range(nbused):
        I_block.append(ib * NKSKIP)
    
    #--- 3. Output detected peaks
    line = SFILE.split("."); line = line[len(line)-2].split("/")
    label = line[len(line)-1]
    OFILE = "%s_PEAKS.txt"%(label)
    output_peaks(OFILE, I_block, PEAK_IXY_ALL, IFREQ, ISED)

    print "All done!\n"
    
