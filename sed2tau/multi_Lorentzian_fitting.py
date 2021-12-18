#!/usr/bin/env python
import numpy as np
from lmfit.models import LorentzianModel, ExponentialModel
import sys

LPARA = ["center", "sigma", "amplitude"]

def Set_Multi_Lorentzian(center, sigma, amp):
    Npeak = len(center)

    #--- prepare label matrix
    LABEL = [];
    for i in range(3):
        LABEL.append("")

    #--- start fitting
    lorentz = []
    for i in range(Npeak):
        PRE = "L%d_"%(i)
        for j in range(3):
            LABEL[j] = "%s%s"%(PRE, LPARA[j])
        lorentz.append(LorentzianModel(prefix=PRE))

        if i == 0:
            pars = lorentz[0].make_params()
            pars[LABEL[0]].set(center[i][0], min=center[i][1], max=center[i][2])
            pars[LABEL[1]].set(sigma[i][0], min=sigma[i][1], max=sigma[i][2])
            pars[LABEL[2]].set(amp[i][0], min=amp[i][1], max=amp[i][2])
            mod = lorentz[0]
        else:
            pars.update(lorentz[i].make_params())
            pars[LABEL[0]].set(center[i][0], min=center[i][1], max=center[i][2])
            pars[LABEL[1]].set(sigma[i][0], min=sigma[i][1], max=sigma[i][2])
            pars[LABEL[2]].set(amp[i][0], min=amp[i][1], max=amp[i][2])
            mod += lorentz[i]

    return mod, pars


def Multi_Lorentzian_Fitting(x, y, center, sigma, amp):
    Npeak = len(center)
    mod, pars = Set_Multi_Lorentzian(center, sigma, amp)
    init = mod.eval(pars, x=x)
    out  = mod.fit(y, pars, x=x)
    comps = out.eval_components(x=x)
    return out


def Extract_Fitting_Result(npeak, out):
    param = np.zeros((npeak, 3))
    LABEL = []
    for i in range(3):
        LABEL.append("")
    for ip in range(npeak):
        PRE = "L%d"%(ip)
        for j in range(3):
            LABEL[j] = "%s_%s"%(PRE, LPARA[j])
            param[ip][j] = float(out.params[LABEL[j]].value)
    return param


'''
xdat, ydat: raw data which mahve multiple peaks.
x_peak: peak position
x_peak[ip] is near xdat[i_near[ip]]. In general, they are NOT exactly the same.
'''
def fitting_long(IBLOCK, xdat, ydat, x_peak, i_peak, PWIDTH, PERROR):

    #--- 0. prepare
    ndat = len(xdat)
    np_tot = len(x_peak)
    
    #'''
    Nwidth = int(PWIDTH / (xdat[1] - xdat[0]))
    if Nwidth < 5:
        PWIDTH = (xdat[1] - xdat[0]) * 5.
        Nwidth = int(PWIDTH / (xdat[1] - xdat[0]))

    #--- 1. divide into blocks: Ip_bloci[nblock][2] (Iinit, Iend) (0 to len(x_peak)-1)
    Ip_block = []; bcount = 0
    Ip_block.append([-1,-1]); Ip_block[0][0] = 0
    for ip in range(1,np_tot,1):
        if x_peak[ip] > x_peak[ip-1] + PWIDTH:
            Ip_block[bcount][1] = ip - 1
            Ip_block.append([ip, -1])
            bcount += 1 
    Ip_block[bcount][1] = np_tot - 1
    nblock = bcount + 1
    
    #--- perform fitting for each block: {0.center, 1.sigma, 2.amplitude}
    param_all = np.zeros((np_tot, 3))
    pcount = 0
    
    #--- 2. Lorentzian fitting for each block
    for iblock in range(nblock):
        I0_raw = i_peak[Ip_block[iblock][0]] - Nwidth
        I1_raw = i_peak[Ip_block[iblock][1]] + Nwidth
        
        if I1_raw < I0_raw:
            print "Error: w0 > w1 (%f > %f) (%d > %d in peak file)"%\
                    (xdat[I0_raw], xdat[I1_raw], I0_raw, I1_raw)
            sys.exit()
        
        #if I0_raw < 0:
        #    I0_raw = 0
        if I1_raw > ndat:
            I1_raw = ndat - 1
                
        #-- fitting data
        nfit = I1_raw - I0_raw + 1
        xfit = np.zeros(nfit); yfit = np.zeros(nfit)
        #FLAG = 0
        for i in range(nfit):
            idata = I0_raw + i
            if idata < 0:
                idata = -idata
                xfit[i] = - xdat[idata]
                yfit[i] = ydat[idata]
                FLAG = 1
            else:
                xfit[i] = xdat[idata]
                yfit[i] = ydat[idata]
        
        '''
        if FLAG == 0:
            for i in range(len(xfit)):
                print "%7.4f %10.4"%(xfit[i], yfit[i])
            sys.exit()
        '''

        #-- initial setting of Lorentzian parameters
        np_loc = Ip_block[iblock][1] - Ip_block[iblock][0] + 1
        center = np.zeros((np_loc, 3))
        sigma = np.zeros((np_loc, 3))
        amp  = np.zeros((np_loc, 3))
        for ip in range(np_loc):
            idat = i_peak[Ip_block[iblock][0] + ip]
            #-- A. center
            center[ip][0] = xdat[idat]
            center[ip][1] = xdat[idat] - PERROR
            center[ip][2] = xdat[idat] + PERROR
            if center[ip][1] < 0.:
                center[ip][1] = 0.
            #-- B. sigma
            sigma[ip][0] = 0.5 / 20.
            sigma[ip][1] = 0.
            sigma[ip][2] = np.inf
            #-- C. amplitude
            amp[ip][0] = ydat[idat] * np.pi * sigma[ip][0]
            amp[ip][1] = 0.
            amp[ip][2] = np.inf

        #--- Lorentzian fitting!!!
        out = Multi_Lorentzian_Fitting(xfit, yfit, center, sigma, amp)

        #--- extract parameters: param_block[np_loc][3]: {0.center, 1.sigma, 2.amplitude}
        param_block = Extract_Fitting_Result(np_loc, out)
        
        for i in range(np_loc):
            param_all[pcount] = param_block[i]
            pcount += 1

    if pcount != np_tot:
        print "Error: # of peaks are inconsistent (%d != %d)"%(pcount, np_tot)
        sys.exit()

    return param_all

