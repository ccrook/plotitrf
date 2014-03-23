#!/usr/bin/python

import sys
import os
import json
from math import *
import numpy as np
from subprocess import call, STDOUT

# Subclass of VectorSet for displaying ITRF vectors


itrfs='ITRF96 ITRF97 ITRF2000 ITRF2005 ITRF2008'.split()
year0=1990
year1=2020

def genData():
    llh=[[0.0,90.0],[0.0,-90.0]];
    for i in (52,10,-27):
        for j in range(0,360,72):
            llh.append([float(j),float(i)])
            llh.append([float(j+36),float(-i)])

    llh.sort(key=lambda x:(-x[1],x[0]))
    pts=llh
    vec=np.zeros((len(llh),3))

    # Create concord input file

    tmpin='qqconcord.in'
    tmpout='qqconcord.out'
    devnull=open(os.devnull,'wb')

    with open(tmpin,"w") as f:
        for l in llh:
            f.write("{0:.1f}\t{1:.1f}\t0.0\n".format(l[0],l[1]))

    # Calculate ITRF XYZ at each epoch for each ITRF

    itrfdata={}
    for itrf in itrfs:
        itrfdata[itrf]=[]
        for epoch in (year0,year1):
            epochstr=str(epoch)+'0101'
            call(('concord',
                  '-iITRF96:enh:d',
                  '-o'+itrf+'_XYZ',
                  '-p5','-y'+epochstr,
                  tmpin,tmpout),
                  stdout=devnull,stderr=STDOUT)
            xyz=[]
            with open(tmpout) as f:
                for l in f:
                    parts=l.split()
                    if len(parts) == 3:
                        xyz.append([float(x) for x in parts])
            itrfdata[itrf].append(xyz)

    os.remove(tmpin)
    os.remove(tmpout)
    itrfdata=itrfdata
    data=dict(
        itrfs=itrfs,
        years=[year0,year1],
        points=pts,
        itrfdata=itrfdata
        )
    # print json.dumps(data,sort_keys=True, indent=4)
    with open("plotitrf.dat","w") as f:
        f.write("ITRFS: {0}\n".format("\t".join(itrfs)))
        f.write("YEARS:\t{0}\t{1}\n".format(year0,year1))
        for i in range(len(llh)):
            f.write("{0}\t{1}".format(llh[i][0],llh[i][1]))
            for itrf in itrfs:
                for iyear in range(2):
                    xyz=itrfdata[itrf][iyear][i]
                    f.write("\t{0}\t{1}\t{2}".format(*xyz))
            f.write("\n")

genData()
