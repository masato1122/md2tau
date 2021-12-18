import os, sys
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def Make_Block_File(OFILE, xdat, ydat):
    
    NN = 1
    yave = np.zeros(NN)
    for i in range(NN):
        yave[i] = average_extopk(ydat, i)
    
    ofs = open(OFILE, "w")
    ofs.write("#Frequency[THz]  SED[-]  SED_average\n")
    for i in range(len(xdat)):
        ofs.write("%13.8e %13.8e "%(xdat[i], ydat[i]))
        for j in range(NN):
            ofs.write("%10.7e "%(yave[j]))
        ofs.write("\n")
    print "Output:", OFILE
    return 0

def Plot_Result(FNAME, xdat, ydat, xpeak, ypeak, param):
    npeak = len(xpeak)
    fig = plt.figure(figsize=matplotlib.figure.figaspect(0.3))

    plt.plot(xdat, ydat, marker=".")

    for ip in range(npeak):

        #--- ver.1
        #xx = np.linspace(xpeak[ip]-0.5*PEAK_WIDTH, xpeak[ip]+0.5*PEAK_WIDTH, 20)
        #yy = param[ip](xx)
        #plt.plot(xx, yy, '-')

        #--- ver.2
        y0 = 1e-2; y1 = 1e2
        plt.plot([xpeak, xpeak], [y0, y1], 'r')
    plt.yscale("log")
    plt.savefig(FNAME, format = 'png', dpi=300)
    print "Output:", FNAME

    return 0

def bubble_sort(array, SORT_TYPE):
    n = len(array)
    asort = np.zeros(n)
    for i in range(n):
        asort[i] = array[i]
    for i in range(n-1):
        for j in range(n-1, i, -1):
            FLAG = 0
            if SORT_TYPE == "UP":
                if asort[j] < asort[j-1]:
                    FLAG = 1
            else:
                if asort[j] > asort[j-1]:
                    FLAG = 1
            if FLAG == 1:
                tmp = asort[j]
                asort[j] = asort[j-1]
                asort[j-1] = tmp
    return asort

def bubble_sort_down_ID(array):
    n = len(array)
    ID = np.arange(n)
    asort = np.zeros(n)
    for i in range(n):
        asort[i] = array[i]
    for i in range(n-1):
        for j in range(n-1, i, -1):
            if asort[ID[j]] > asort[ID[j-1]]:
                itemp = ID[j]
                ID[j] = ID[j-1]
                ID[j-1] = itemp
    return ID

def average_extopk(data, ntop):
    ndat = len(data)
    if ntop > ndat-1:
        print "Error during taking average excluding top-k data"
        print "Change ntop value (%d to %d)"%(ntop, ndat)
        ntop = ndat
    dtemp = np.zeros(ndat)
    dtemp[:] = data[:]
    dtemp = bubble_sort(dtemp, "up")
    average = np.average(dtemp[0:ndat-1-ntop])
    return average

#----------- file handling
def mkinputfile(FNAME):
    ofs = open(FNAME, "w")
    ofs.write("#------------------------------------------------------\n")
    ofs.write("# 1.Shared\n")
    ofs.write("# NPROCS:   # of parallelization\n")
    ofs.write("# SEDFILE:  SED file name\n")
    ofs.write("# IFREQ, ISED: Line number for frequency and SED (>=1)\n")
    ofs.write("# FREQUNIT: Hz or THz\n")
    ofs.write("#\n")
    ofs.write("# 2. for peak_detecter\n")
    ofs.write("# NKSKIP: each NKSKI kpoints will be analyzed\n")
    ofs.write("# PEAKMIN_FLAG: 0. use PEAKMIN, 1. use top-{PEAKMIN} data\n")
    ofs.write("# PEAKMIN     : will be used to determin min SED for peaks\n")
    ofs.write("# PEAKWIDTH   : min. distance btw/ peaks\n")
    ofs.write("#\n")
    ofs.write("# 3. for peak_analyzer\n")
    ofs.write("# JOBTYPE:  plot, fitting, or plot_lorentz\n")
    ofs.write("# PEAKFILE: peak file name\n")
    ofs.write("#-----------------------------------------------------\n")
    ofs.write("\n")
    ofs.write("#----- For peak_analyzer (1.plot, 2.fitting, and 3.plot_lorentz)\n")
    ofs.write("#-- 1. Plot raw data and initial peak positions\n")
    ofs.write("#JOBTYPE   plot\n")
    ofs.write("#PEAKFILE  SED_rtz1_PEAKS.txt\n")
    ofs.write("\n") 
    ofs.write("#-- 2. Do fitting!\n")
    ofs.write("#JOBTYPE   fitting\n")
    ofs.write("#PEAKFILE  SED_rtz1_PEAKS.txt\n")
    ofs.write("\n") 
    ofs.write("#-- 3. Plot result of Lorentzian fitting\n")
    ofs.write("#JOBTYPE   plot_lorentz\n")
    ofs.write("#PEAKFILE  SED_rtz1_PEAKS_LORENTZ.txt\n")
    ofs.write("\n")
    ofs.write("#----- for parallel calculation\n")
    ofs.write("#NPROCS    2\n")
    ofs.write("NKSKIP    12\n")
    ofs.write("\n")
    ofs.write("#--- for SED file\n")
    ofs.write("#SEDFILE   SED_rtz1.txt\n")
    ofs.write("IFREQ     3\n")
    ofs.write("ISED      4\n")
    ofs.write("FREQUNIT  THz\n")
    ofs.write("\n")
    ofs.write("#--------- For peak_detecter\n")
    ofs.write("#- If PEAKMIN_FLAG == 0, min.SED for peaks is directly PEAKMIN.\n")
    ofs.write("#- If PEAKMIN_FLAG == 1, min.SED of top **% of all data will be used.\n")
    ofs.write("#- If PEAKMIN_FLAG == 2, min.SED of top **% of data at each k-point will be used.\n")
    ofs.write("#PEAKMIN_FLAG  0\n")
    ofs.write("#PEAKMIN       10.0\n")
    ofs.write("\n")
    ofs.write("PEAKMIN_FLAG  2\n")
    ofs.write("PEAKMIN       0.02\n")
    ofs.write("\n")
    ofs.write("PEAKWIDTH     0.1\n")
    ofs.write("FITTINGWIDTH  0.2\n")
    ofs.write("LORENTZIANBUFFER  0.05\n")
    ofs.write("\n")
    ofs.close()
    print "See %s"%(FNAME)
    return 0

def read_input_file(IFILE, SLINE, I_priority):
    ifs = open(IFILE, "r")
    lines = ifs.readlines()
    ifs.close()
    for line in lines:
        data = line.split()
        if len(data) <= 1: continue
        if line[0] == "#": continue
        if data[0] == SLINE:
            return data[1]
    if I_priority == 1:
        print "cannot find %s at %s"%(SLINE, IFILE)
        sys.exit()
    return -1
    
