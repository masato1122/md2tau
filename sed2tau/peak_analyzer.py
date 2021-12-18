#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os, sys
import os.path
import numpy as np
from optparse import OptionParser

#--- for parallel calculation
import multiprocessing as mp
def argwrapper(args):
    return args[0](*args[1:])

#--- import my files
import tips
import read_file
import multi_Lorentzian_fitting
import plot_data

#-- directory for the log
#SED_DIR = "SED_DATA"
LOG_DIR = "./LOG_PEAK"
if os.path.exists(LOG_DIR):
    pass
else:
    command = "mkdir -p %s"%(LOG_DIR)
    os.system(command)

#--- check file
def file_check(file):
    if os.path.exists(file) == False:
        print "Cannot find %s"%file
        sys.exit()

#--- function for parallel calculation (a little messy)
'''
- ib_det: 0,1,2,..., Iblock[nb_det]: 0, m, 2m, 3m, ...
- params = [FILE_FLAG, IFREQ, ISED, JTYPE, WUNIT, PWIDTH, ILABEL4PEAK, IFREQ4PEAK]
  FILE_FLAG: 0. _PEAKS.txt file, 1. ./LOG_PEAK/block*.txt file
- sed_data[nblock_all, npeak, nlines+2]: If FILE_FLAG==1, dummpy (== -1).
/// for peak file
- peak_data[nb_detected][npeak][nlines] {0.label, 1.ID, 2-. data}
- Iblock[nb_detected]
'''
def get_peak_nearestID(xpeak, xdat):
    npeak = len(xpeak)
    ndat = len(xdat)
    I_PEAK = [int(npeak) for i in range(npeak)]
    iinit = 0; pcount = 0
    for ip in range(npeak):
        for idat in range(iinit, ndat-1, 1):
            if xdat[idat] <= xpeak[ip] and xpeak[ip] <= xdat[idat+1]:
                if abs(xpeak[ip] - xdat[idat]) <= abs(xdat[idat+1] - xpeak[ip]):
                    I_PEAK[ip] = idat
                else:
                    I_PEAK[ip] = idat + 1
                iinit = idat
                break
    return I_PEAK

def parallel_calculation(IBLOCK_DET, params, sed_data, Iblock, peak_data_iblock):
    
    FILE_FLAG  = int(params[0])
    IFREQ      = int(params[1])
    ISED       = int(params[2])
    JOB_TYPE   = params[3]
    WUNIT      = params[4]
    PWIDTH     = params[5]
    PERROR_LORENTZ = params[6]
    ILABEL4PEAK = int(params[7])
    IFREQ4PEAK  = int(params[8])
    SED_DIR     = params[9]
    
    #print Iblock
    #print IBLOCK_DET
    #sys.exit()
    
    #-- 1-1. get peak data: {L, W}_PEAK[nb_det][npeak] (L: label, W: frequency)
    L_PEAK = []; W_PEAK= []
    npeak = len(peak_data_iblock)
    for ip in range(npeak):
        #print IBLOCK_DET, ip, peak_data_iblock[ip][IFREQ4PEAK]
        L_PEAK.append(str(peak_data_iblock[ip][ILABEL4PEAK]))
        W_PEAK.append(float(peak_data_iblock[ip][IFREQ4PEAK])) 
    
    #-- 2. get raw data
    if FILE_FLAG == 0:
        ndat = len(sed_data[Iblock[IBLOCK_DET], :, IFREQ])
        xdat = np.zeros(ndat); ydat = np.zeros(ndat)
        xdat = sed_data[Iblock[IBLOCK_DET], :, IFREQ]
        ydat = sed_data[Iblock[IBLOCK_DET], :, ISED]
    elif FILE_FLAG == 1:
        BFILE = "%s/block%d.txt"%(SED_DIR, Iblock[IBLOCK_DET])
        IX = IFREQ
        IY = ISED
        xdat, ydat = read_file.read_block_file(BFILE, IX, IY)
        if WUNIT == "Hz":
            xdat *= 1e-12
    
    #--- 3. get Ipeak (peak ID in raw data)
    I_PEAK = get_peak_nearestID(W_PEAK, xdat)
    
    #--- band ID
    iband = 0
    if IBLOCK_DET == 0:
        iband = 0
    else:
        for i in range(IBLOCK_DET):
            if Iblock[i+1] <= Iblock[i]:
                iband += 1
    #print Iblock
    #sys.exit()
    #------------------------- Start each job --------------------------------#
    #--- Case1. Plot
    if JOB_TYPE == "plot":

        fname = "%s/block_%d-%d.png"%(LOG_DIR, iband, Iblock[IBLOCK_DET])
        '''
        if IBLOCK_DET != len(Iblock)-1 and Iblock[IBLOCK_DET] == Iblock[IBLOCK_DET+1]:
            fname = "%s/block%d-1.png"%(LOG_DIR, Iblock[IBLOCK_DET])
        if IBLOCK_DET != 0 and Iblock[IBLOCK_DET] == Iblock[IBLOCK_DET-1]:
            fname = "%s/block%d-2.png"%(LOG_DIR, Iblock[IBLOCK_DET])
        '''
        plot_data.Plot_Data_Peaks(fname, xdat, ydat, L_PEAK, W_PEAK, I_PEAK)
        return 0

    #--- Case2. Lorentzian Fitting
    elif JOB_TYPE == "fitting":
        #-- param_all[npeak][3] (0.center, 1.sigma, 2.amplitude)
        param_all = multi_Lorentzian_fitting.fitting_long(IBLOCK_DET, xdat, ydat, W_PEAK, I_PEAK, PWIDTH, PERROR_LORENTZ)
        return param_all
    
    #--- Case3. Plot results of Lorentzian fitting
    elif JOB_TYPE == "plot_lorentz":
        #-------------------
        ISIGMA4LOR = 4 - 1
        IAMP4LOR = 5 - 1
        #-------------------

        #--- get parameters of Lorentzian fitting from peak data
        #--- {S,A}_LOR[npeak]. center is W_PEAK
        S_LOR = []; A_LOR = []
        for ip in range(npeak):
            S_LOR.append(float(peak_data_iblock[ip][ISIGMA4LOR]))
            A_LOR.append(float(peak_data_iblock[ip][IAMP4LOR]))
        
        #--- plot dat with Lorentzian function
        fname = "%s/LOR_block_%d-%d.png"%(LOG_DIR, iband, Iblock[IBLOCK_DET])
        '''
        if IBLOCK_DET != len(Iblock)-1 and Iblock[IBLOCK_DET] == Iblock[IBLOCK_DET+1]:
            fname = "%s/LOR_block%d-1.png"%(LOG_DIR, Iblock[IBLOCK_DET])
        if IBLOCK_DET != 0 and Iblock[IBLOCK_DET] == Iblock[IBLOCK_DET-1]:
            fname = "%s/LOR_block%d-2.png"%(LOG_DIR, Iblock[IBLOCK_DET])
        '''
        plot_data.Plot_Data_Lorentzian(fname, xdat, ydat, L_PEAK, W_PEAK, I_PEAK, S_LOR, A_LOR)
        
        
    else:
        print "Error, job type, %s, is not defined."%(JOB_TYPE)
        sys.exit()

    return 0

#------ Output fitting results: Iblock[nblock], result[nblock][npeak][3]
def output_fitting_results(RFILE, Iblock, result, L2_peak):
    ofs = open(RFILE, "w")
    ofs.write("#ID  dummy  center[THz]  sigma   amplitude   tau[ps]\n")
    nblock = len(Iblock)
    for ib in range(nblock):
        ofs.write("# Iblock: %d\n"%(Iblock[ib]))
        np_loc = len(result[ib])
        for ip in range(np_loc):
            tau = 0.5 / result[ib][ip][1]   #-- tau[ps] = 1/2/sigma 
            ofs.write("%d  %s  "%(ip+1, L2_peak[ib][ip]))

            #--- (center, sigma, amp, tau)
            for j in range(3):
                ofs.write("%15.7e  "%(result[ib][ip][j]))
            ofs.write("%15.7e  "%(tau))
            ofs.write("\n")
        ofs.write("\n")

    print "Output:", RFILE
    ofs.close()
    return 0

#--------- START MAIN ------------
if __name__ == "__main__":

    print "------------------------------------------"
    print ""
    print " peak_analyzer (plot or fitting) ver.1.0"
    print ""
    print "------------------------------------------"

    #----- get input file name
    parser = OptionParser()
    parser.add_option("-f", "--ifile", dest="ifile", type="string", \
            help="Input input file name.")
    
    #- optional: SED file, peak file, and job type
    parser.add_option("-s", "--sedfile", dest="sfile", type="string", \
            help="Input SED file name.")
    parser.add_option("--seddir", dest="sdirectory", type="string", \
            help="Directory where SED block files are located")
    parser.add_option("-p", "--peakfile", dest="pfile", type="string", \
            help="Input peak file name.")
    parser.add_option("-j", "--jobtype", dest="jtype", type="string", \
            help="plot, fitting, or plot_lorentz")
    parser.add_option("-n", "--nprocs", dest="nprocs", type="int", \
            help="Number of parallelization")
    (options, args) = parser.parse_args()
    if options.ifile == None:
        tips.mkinputfile("peak_example.in")
        print "Input input file name (-f, --ifile)."
        sys.exit()
    
    IFILE = options.ifile
    file_check(IFILE)

    #----- read input parameters
    #- for SED file
    if options.sfile == None:
        SFILE = tips.read_input_file(IFILE, "SEDFILE", 1)
    else:
        SFILE = options.sfile
    if options.sdirectory == None:
        SED_DIR = tips.read_input_FILE(IFILE, "SED_DIR", 1)
    else:
        SED_DIR = options.sdirectory
    data  = tips.read_input_file(IFILE, "IFREQ", 1);    IFREQ  = int(data) - 1
    data  = tips.read_input_file(IFILE, "ISED", 1);     ISED   = int(data) - 1
    WUNIT = tips.read_input_file(IFILE, "FREQUNIT", 1)
    #- to make the calculation fast
    if options.nprocs == None:
        data  = tips.read_input_file(IFILE, "NPROCS", 0);   nprocs = int(data)
    else:
        nprocs = options.nprocs
    data  = tips.read_input_file(IFILE, "NKSKIP", 1);   NKSKIP = int(data)
    
    #----- job type
    if options.jtype == None:
        JTYPE = tips.read_input_file(IFILE, "JOBTYPE", 1)
    else:
        JTYPE = options.jtype
    
    #--- peak file
    if options.pfile == None:
        PFILE = tips.read_input_file(IFILE, "PEAKFILE", 1)
    else:
        PFILE = options.pfile
    
    file_check(SFILE)
    file_check(PFILE)

    #--- output
    if nprocs <= 1:
        print ">> Serical calculation."
    else:
        print ">> Parallel calculation (%d)"%(nprocs)
    print ""
    print "=====  JOB TYPE: %s  ====="%(JTYPE)
    print ""

    if JTYPE == "fitting":
        data = tips.read_input_file(IFILE, "FITTINGWIDTH", 1)
        PWIDTH = float(data)
        data = tips.read_input_file(IFILE, "LORENTZIANBUFFER", 1)
        PERROR = float(data)
    else:
        PWIDTH = 0.0    #--- This will NOT be used
        PERROR = float(data)
    
    #--- 1. read PEAK file
    '''
    Iblock[nblock_detected]  <-- "# Iblock: *" (at _PEAKS.txt file)
    peak_data[nb_detected][npeak][nlines+2]
    {0.label, 1.data id, 2- data  }
    >>>
    L_PEAK[nblock_detected][npeak]  [THz]
    W_PEAK[nblock_detected][npeak]
    '''
    print ">> Read peak file, %s"%(PFILE)
    Iblock, peak_data = read_file.read_SED_peaks(PFILE)
    nb_det = len(Iblock)

    #-------------------#
    ILABEL4PEAK = 1 - 1
    IFREQ4PEAK  = 3 - 1
    #-------------------#
    '''
    #- 1-2. take data
    L_PEAK = []; W_PEAK= []  #-- mat[nb_det][npeak]
    for ib in range(nb_det):
        L_PEAK.append([]); W_PEAK.append([])
        npeak = len(peak_data[ib])
        for ip in range(npeak):
            W_PEAK[ib].append(float(peak_data[ib][ip][IFREQ4PEAK]))
            L_PEAK[ib].append(str(peak_data[ib][ip][ILABEL4PEAK]))
    '''
    #--- 2-0. check ./LOG_PEAK/block#.txt
    FILE_FLAG = 1
    for ib in range(nb_det):
        fcheck = "%s/block%d.txt"%(SED_DIR, Iblock[ib])
        if os.path.exists(fcheck) == False:
            print "cannot find %s"%(fcheck)
            FILE_FLAG = 0
            break
    
    #--- 2-1. ./SED_PEAK/block*.txt files will be used.
    if FILE_FLAG == 1:
        sed_data = -1    #-- make dummy
        print "%s/block*.txt will be used."%(SED_DIR)
        print ""
        
    #--- 2-2. read SED file. if ./${SED_DIR}/block#.txt does NOT exit.
    elif FILE_FLAG == 0:
        '''
        sed_data[nb_all,ndata,nlines]
        '''
        print ""
        print ">> Read %s"%(SFILE)
        print "========== CAUTION!! =========="
        print "  Cannot find %s/block#.txt"%(SED_DIR)
        print "  This may use large memory."
        print "  All data at %s will be read."%(SFILE)
        print "==============================="
        sed_data = read_file.read_SED_file(SFILE)
        nb_all = len(sed_data)
        print "Nblock(all): ", nb_all

    #////// OBTAINED DATA ///////#
    '''
    1. Iblock[nblock][npeak][nlines+2], {0.label of peak, 1.data iD, 2- data}

    2-1. If ./LOG_PEAK/block#.txt does NOT exit (FILE_FLAG == 0),
    - sed_data[nb_all,ndata,nlines]

    2-2. otherwise, ./LOG_PEAK/block#.txt will be read.
    '''
    
    #----- 3. PLOT or FITTING
    params = [FILE_FLAG, IFREQ, ISED, JTYPE, WUNIT, PWIDTH, PERROR, ILABEL4PEAK, IFREQ4PEAK, SED_DIR]
    
    #-- 3-1. parallel
    if nprocs > 1:
        pool = mp.Pool(nprocs)
        func_args = []
        for ik in xrange(nb_det):
            func_args.append((parallel_calculation, ik, params, sed_data, Iblock, peak_data[ik]))
        result = pool.map(argwrapper, func_args)
    
    #-- 3-2. serial
    else:
        result = []
        for ik in xrange(nb_det):
            print "IK: ", ik
            result.append(parallel_calculation(ik, params, sed_data, Iblock, peak_data[ik]))
    
    #------ 4. For fitting: result[nb_det][npeak(ib)][3]
    if JTYPE == "fitting":
        
        ILABEL2 = 2 - 1

        #--- label2: L2_peak[nb_det][npeak]
        L2_peak = []
        for ik in range(nb_det):
            L2_peak.append([])
            npeak = len(peak_data[ik])
            for ip in range(npeak):
                label = peak_data[ik][ip][ILABEL2]
                L2_peak[ik].append(label)
        
        #--- output
        line = PFILE.split("."); 
        if len(line) < 2:
            label = "DUMMY"
        else:
            label = line[len(line)-2]
        RFILE = "%s_LORENTZ.txt"%(label)
        output_fitting_results(RFILE, Iblock, result, L2_peak)

    print "All done!\n"
     
