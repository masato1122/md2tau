import sys
import numpy as np
import tips

#--- 1. quadratic fitting
'''
- using width of peak, quadratic fitting
'''
def quadratic(xdat, ydat, xwidth):

    ndat = len(xdat)
    nwidth = int(xwidth / (xdat[1] - xdat[0]))
    if nwidth%2 == 0:
        nwidth += 1
    nhalf = (nwidth - 1) / 2

    print "Nwidth:", nwidth

    #--- prepare
    xlocal = np.zeros(nwidth)
    ylocal = np.zeros(nwidth)
    xpeak = []; ypeak = []; pcount = 0
    ALLPARAM = []

    #--- detect!
    for idat in range(nhalf, ndat-nhalf-1, 1):
        for i in range(nwidth):
            xlocal[i] = xdat[idat-nhalf+i]
            ylocal[i] = ydat[idat-nhalf+i]

        param = np.poly1d( np.polyfit(xlocal, ylocal, 2) )
        xcent = - param[1]/2./param[2]

        if param[2] < 0 and xlocal[0] < xcent and xcent < xlocal[nwidth-1]:
            if pcount != 0:
                if abs(xcent-xpeak[pcount-1] < xwidth):
                    continue
            xpeak.append(xcent)
            ycent = param(xcent)
            ypeak.append(ycent)
            ALLPARAM.append(param)
            pcount += 1

    return xpeak

#----- 2. use threshold
'''
- using ymin for peak
'''
'''
def get_peak_min(ydat, topratio): 
    ID_sort = tips.bubble_sort_down_ID(ydat)
    iget = int(len(ydat) * topratio) 
    return ydat[ID_sort[iget]]
'''
#---- ver.1: first version
def ymin(xdat, ydat, peak_width, minSED4PEAK):
    '''
    OFILE = "check.txt"
    ofs = open(OFILE, "w")
    for i in range(len(xdat)):
        ofs.write("%f %f %f\n"%(xdat[i], ydat[i], minSED4PEAK))
    '''
    
    #--- 0-2. peak width
    ndat = len(xdat)
    Nbuffer = int(peak_width/(xdat[1] - xdat[0]))
    if Nbuffer == 0:
        Nbuffer = 1

    #--- 1. get ID of data exceeding peak_min
    Iexceed = []; icount = 0; bcount = 0
    for idat in range(ndat):
        if ydat[idat] > minSED4PEAK:
            Iexceed.append(idat)
            if icount > 0:
                if Iexceed[icount] > Iexceed[icount-1] + Nbuffer:
                    bcount += 1
            icount += 1
    nex = icount
    nblock = bcount + 1
    if nblock == 0:
        return []

    #--- 2. get peak position: Iexceed[] => Ipeak[]
    Ipeak = []; pcount = 0;
   
    #-- ver.1
    ymax = -100.; FLAG = 0
    Itemp = -1
    for idat in range(nex):
        
        #- new peak position
        if ydat[Iexceed[idat]] > ymax:
            Itemp = Iexceed[idat]
            ymax = ydat[Iexceed[idat]]

        #- determin peak posiiton
        if idat < nex-1:
            if Iexceed[idat+1] - Iexceed[idat] > Nbuffer:
                Ipeak.append(Itemp)
                ymax = -100.; Itemp = -1
                pcount += 1
        else:
            Ipeak.append(Itemp)
            ymax = -100.; Itemp = -1
            pcount += 1

    #-- confirm
    if pcount != nblock:
        print "Error: %d != %d"%(pcount, nblock)
        sys.exit()

    return Ipeak


#----
#def newalgorism(xdat, ydat, peak_width, minSED4PEAK):
#    return 0

#----- 3. remove singurality points
def remove_singurality(xdat, ydat, Iinput):
    return 0 
