#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os, sys
import numpy as np

#--- for matplotlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
from matplotlib import cm
#import seaborn as sns

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

def bubble_sort_down(xdat, ydat):
    ndat = len(xdat)
    newdat = np.zeros((ndat, 2))
    for i in range(ndat):
        for j in range(2):
            newdat[i][0] = xdat[i]
            newdat[i][1] = ydat[i]
    
    newdat = sorted(newdat, key=lambda x: x[1], reverse=True)
    xnew = np.zeros(ndat)
    ynew = np.zeros(ndat)
    for i in range(ndat):
        xnew[i] = newdat[i][0]
        ynew[i] = newdat[i][1]
    return xnew, ynew

#-------------------------
NUM_MULTIPLE = 5
NUM_SECTION = 500
#-------------------------

def cal_ID_col(ID_sort, Ncol):
    ndat = len(ID_sort)
    Neach = int(ndat / NUM_SECTION)
    ID_col = [int(0) for i in range(ndat)]
    for ii in range(ndat):
        num = ID_sort[ii]
        ID_col[num] = Ncol - int(num / Neach)
        if ID_col[num] < 0.:
            ID_col[num] = 0.
    return ID_col

def cal_ID_col2(ndat, Ncol):
    Neach = int(ndat / NUM_SECTION)
    ID_col = np.arange(ndat)
    for ii in range(ndat):
        ID_col[ii] = Ncol - int(ii / Neach)
        if ID_col[ii] < 0.:
            ID_col[ii] = 0.
    return ID_col

def cal_ID_plot(ndat, nskip):
    II = []; count = 0
    for i in range(ndat):
        if i <= ndat * 0.1:
            II.append(i)
            count += 1
        else:
            if i%nskip == 0:
                II.append(i)
                count +=1
    IDplot = np.zeros(count, dtype=int)
    for i in range(count):
        IDplot[i] = II[i]
    return IDplot


#--- plot data
def Plot_Data_Peaks(FNAME, xdat, ydat, L_PEAK, W_PEAK, I_PEAK):

    ndat = len(xdat)
    npeak = len(L_PEAK)
    ymin = np.amin(ydat) * 0.1
    ymax = np.amax(ydat) * 2.0

    #---------
    Ndiv = NUM_MULTIPLE
    #---------
    xmax = 55.

    #--- size & aspect ratio
    ratio = 0.2 * float(Ndiv)
    xwidth = 40/Ndiv + 2
    ywidth = xwidth * 0.2 * Ndiv
    fig = plt.figure(figsize=(xwidth, ywidth))

    #-- sort ydata for color
    Ncol = 11
    PER_SEC = 100./float(NUM_SECTION) #--- (100/Nsec)% data is in a section.
    Neach = int(ndat / NUM_SECTION)

    #--- ver.1
    #ID_sort = bubble_sort_down_ID(ydat)
    #ID_col = cal_ID_col(ID_sort, Ncol)
    
    #--- ver.2
    xsort, ysort = bubble_sort_down(xdat, ydat)
    ID_col = cal_ID_col2(len(xsort), Ncol)

    
    #-- skip data
    nplot = 2000; nskip = int(len(xsort) / nplot)
    ID_plot = cal_ID_plot(len(xsort), nskip) 
    '''
    if nskip == 0:
        nskip = 1
    if nplot > len(xsort):
        nplot = len(xsort)
    ID_plot = np.zeros(nplot, dtype=int)
    for i in range(nplot):
        ID_plot[i] = i * nskip
    '''
    
    plt.rcParams["font.size"] = 10
    
    #--- multiple plot
    for isec in range(Ndiv):
        plt.subplot(Ndiv,1,isec+1)
        x0 = xmax * float(isec)/ float(Ndiv) - 0.2
        x1 = xmax * float(isec+1)/ float(Ndiv) + 0.2
        if isec == Ndiv - 1:
            x1 = 52.

        ax = fig.add_subplot(Ndiv,1,isec+1)
        ax.grid(color='grey', axis='both', linestyle='solid', linewidth=0.4, which='major')
        ax.grid(color='grey', axis='x', linestyle='dashed', linewidth=0.2, which='minor')
        plt.gca().xaxis.set_major_locator(tick.MultipleLocator(0.5))
        plt.gca().xaxis.set_minor_locator(tick.MultipleLocator(0.1))
        plt.tick_params(axis='both', which='both', direction='in')
        plt.xlim([x0, x1])

        #-- ylabel
        plt.ylim([ymin, ymax])
        plt.yscale("log")
        plt.ylabel('SED (-)')

        #-- 2. raw dat
        cm = plt.cm.get_cmap('rainbow', Ncol)
        #sc = plt.scatter(xdat[ID_sort2], ydat[ID_sort2], c=ID_col2, s=5, cmap=cm)
        sc = plt.scatter(xsort[ID_plot], ysort[ID_plot], c=ID_col[ID_plot], s=5, cmap=cm)
        if isec == Ndiv-1:
            cbar = plt.colorbar(sc)
            #cbar.set_label(size=14)
            cbar.ax.set_ylabel('%4.2f%% data \nin a section'%(PER_SEC), fontsize=10)
            plt.xlabel('Frequency (THz)')

        #--- 3. peaks
        pcount = 0
        for ip in range(npeak):
            if W_PEAK[ip] <= x0 or x1 <= W_PEAK[ip]:
                continue
            c = plt.cm.Vega10(float(ip) / float(npeak))
            plt.axvline(W_PEAK[ip], linewidth=1.0, color=c, label=I_PEAK[ip])
            plt.text(W_PEAK[ip], ydat[I_PEAK[ip]], L_PEAK[ip], ha = 'left', va = 'bottom', size=6)
            pcount += 1
        if pcount != 0:
            plt.legend(ncol=18, loc='lower left', handlelength=0.5, fontsize=10)

    plt.savefig(FNAME, format = 'png', dpi=200)
    while os.path.exists(FNAME) == False:
        pass
    print "Output:", FNAME
    return 0

#--- plot data with Lorentzian function
def cal_lorentzian(x, center, sigma, amp):
    p1 = amp / sigma / np.pi
    p2 = 1. + np.power((x - center) / sigma, 2)
    return p1 / p2

def Plot_Data_Lorentzian(FNAME, xdat, ydat, L_PEAK, W_PEAK, I_PEAK, S_LOR, A_LOR):

    ndat = len(xdat)
    npeak = len(L_PEAK)
    ymin = np.amin(ydat) * 0.1
    ymax = np.amax(ydat) * 2.0

    #---------
    Ndiv = NUM_MULTIPLE
    #---------

    xmax = 55.
    
    #--- cal. tau [ps]
    TAU_LABEL = []
    for i in range(npeak):
        tau = "%.0fps"%(1./2./S_LOR[i])
        TAU_LABEL.append(tau)

    #--- size & aspect ratio
    ratio = 0.2 * float(Ndiv)
    xwidth = 40/Ndiv + 2
    ywidth = xwidth * 0.2 * Ndiv
    fig = plt.figure(figsize=(xwidth, ywidth))

    #-- sort ydata for color
    Ncol = 11
    PER_SEC = 100./float(NUM_SECTION) #--- (100/Nsec)% data is in a section.
    Neach = int(ndat / NUM_SECTION)


    #--- ver.1
    #ID_sort = bubble_sort_down_ID(ydat)
    #ID_col = cal_ID_col(ID_sort, Ncol)
    
    #--- ver.2
    xsort, ysort = bubble_sort_down(xdat, ydat)
    ID_col = cal_ID_col2(len(xsort), Ncol)
     
    #--- skip data
    nplot = 2000; nskip = int(len(xsort) / nplot)
    ID_plot = cal_ID_plot(len(xsort), nskip)

    '''
    if nskip == 0: nskip = 1
    if nplot > len(xsort):
        nplot = len(xsort)
    ID_plot = [int(0) for i in range(nplot)]
    for i in range(nplot):
        ID_plot[i] = i * nskip
    '''

    plt.rcParams["font.size"] = 10
    
    #--- multiple plot
    for isec in range(Ndiv):
        plt.subplot(Ndiv,1,isec+1)
        x0 = xmax * float(isec)/ float(Ndiv) - 0.2
        x1 = xmax * float(isec+1)/ float(Ndiv) + 0.2
        if isec == Ndiv - 1:
            x1 = 52

        ax = fig.add_subplot(Ndiv,1,isec+1)
        ax.grid(color='grey', axis='both', linestyle='dashed', linewidth=0.2)
        plt.gca().xaxis.set_major_locator(tick.MultipleLocator(1))
        plt.gca().xaxis.set_minor_locator(tick.MultipleLocator(0.1))
        plt.tick_params(axis='both', which='both', direction='in')
        plt.xlim([x0, x1])

        #-- ylabel
        #plt.ylim([ymin, ymax])
        plt.yscale("log")
        plt.ylabel('SED (-)')

        #-- 2. raw dat
        cm = plt.cm.get_cmap('rainbow', Ncol)
        #sc = plt.scatter(xdat[ID_sort2], ydat[ID_sort2], c=ID_col2, s=5, cmap=cm)
        sc = plt.scatter(xsort[ID_plot], ysort[ID_plot], c=ID_col[ID_plot], s=2, cmap=cm)
        if isec == Ndiv-1:
            cbar = plt.colorbar(sc)
            cbar.ax.set_ylabel('%4.2f%% data \nin a section'%(PER_SEC), fontsize=10)
            plt.xlabel('Frequency (THz)')

        #--- 3. peaks
        pcount = 0
        for ip in range(npeak):
            if W_PEAK[ip] <= x0 or x1 <= W_PEAK[ip]:
                continue
            c = plt.cm.Vega10(float(ip) / float(npeak))
            
            #--- local Lorentzian curve
            xlor_min = W_PEAK[ip] - 0.3; xlor_max = W_PEAK[ip] + 0.3
            if xlor_min < x0:
                xlor_min = x0
            if xlor_max > x1:
                xlor_max = x1
            xlor = np.linspace(xlor_min, xlor_max, 100)
            ylor = cal_lorentzian(xlor, W_PEAK[ip], S_LOR[ip], A_LOR[ip])
           
            xlabel = xlor[int(len(xlor) * 0.6)]
            ylabel = ylor[int(len(xlor) * 0.6)]
            
            #--- plot
            plt.plot(xlor, ylor, color=c, label=L_PEAK[ip])
            plt.text(xlabel, ylabel, TAU_LABEL[ip], ha='left', va='bottom', size=5)
            pcount += 1

        #if pcount != 0:
        #    plt.legend(ncol=18, loc='lower left', handlelength=0.5, fontsize=10)
        
    plt.savefig(FNAME, format = 'png', dpi=200)
    while os.path.exists(FNAME) == False:
        pass
    print "Output:", FNAME
    return 0



