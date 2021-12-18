#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os, sys
import numpy as np

#----- read SED file and output dall[nblock,ndata,nlines]
def read_SED_file(SFILE):
    Nline = sum(1 for line in open(SFILE))
    ifs = open(SFILE, "r")
    dall = []; bcount = 0; dcount = 0
    dall.append([])

    #--- read file
    for il in range(Nline):
        line = ifs.readline()
        data = line.split()
        if len(data) == 0:
            dall.append([])
            bcount += 1
            if il > Nline - 5: break
            if dcount == 0: break
            dcount = 0
            continue
        if line[0] == "#":
            continue

        #-- read data
        dall[bcount].append([])
        nline = len(data)
        for j in range(nline):
            dall[bcount][dcount].append(float(data[j]))
        dcount += 1

        #-- for the case where the last line is not blank.
        if il == Nline-1:
            bcount += 1
    ifs.close()

    #--- for numpy style
    nblock = bcount
    ndeach = dcount
    all_data = np.zeros(nblock*ndeach*nline).reshape(nblock, ndeach, nline)
    for ib in range(nblock):
        for ii in range(ndeach):
            for il in range(nline):
                all_data[ib, ii, il] = dall[ib][ii][il]
    return all_data

#----- read peak file
'''
Output: 
- Iblock[nblock]: k-point ID
- pdata[nblock][ndata]
(0.peak id, 1.data id at raw data, 2-. dumped data)
'''
def read_SED_peaks(PFILE):
    ifs = open(PFILE, "r")
    LALL = ifs.readlines()
    ifs.close()

    Iblock = []; pdata = []; bcount = -1
    for line in LALL:
        data = line.split()
        if len(data) == 0:
            continue
        if line.find("Iblock:") != -1:
            Iblock.append(int(data[2]))
            pdata.append([])
            pcount = 0
            bcount += 1
            continue
        elif line[0] == "#":
            continue
        
        #-- get data
        pdata[bcount].append(data)
        pcount += 1

    return Iblock, pdata

#------ read block data: ./LOG_PEAK/block#.txt
def read_block_file(BFILE, IX, IY):
    ifs = open(BFILE, "r")
    xx = []; yy = []; count = 0
    while ifs:
        line = ifs.readline()
        data = line.split()
        if len(data) == 0: break
        if line[0] == "#": continue
        xx.append(float(data[IX]))
        yy.append(float(data[IY]))
        count += 1
    ifs.close()

    xdat = np.zeros(count)
    ydat = np.zeros(count)
    for i in range(count):
        xdat[i] = xx[i]
        ydat[i] = yy[i]
    
    return xdat, ydat



