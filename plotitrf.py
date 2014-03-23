#!/usr/bin/python

import sys
import os
import os.path
from math import *
import numpy as np
from subprocess import call, STDOUT

from PyQt4.QtCore import *
from PyQt4.QtGui import *
from QtWorld import QtWorld, VectorSet


# Subclass of VectorSet for displaying ITRF vectors

class ITRFVectors( VectorSet ):

    itrfs='ITRF96 ITRF97 ITRF2000 ITRF2005 ITRF2008'.split()
    year0=1990
    year1=2020

    def __init__( self ):
        self.loadData()
        VectorSet.__init__(self,self.pts,self.vec)
        self.setScale(20.0)
        self.setScaleUnits('mm',0.001)
        self.calcVectors( 'ITRF96','ITRF2008',2000.0)

    def loadData(self):
        filename=os.path.join(os.path.dirname(__file__),'plotitrf.dat')
        with open(filename) as f:
            itrfs=f.readline().split()
            if itrfs.pop(0) != 'ITRFS:':
                raise RuntimeError('Invalid ITRFS line in '+filename)
            self.itrfs=itrfs
            years=f.readline().split()
            if years.pop(0) != 'YEARS:' or len(years) != 2:
                raise RuntimeError('Invalid YEARS line in '+filename)
            self.year0=float(years[0])
            self.year1=float(years[1])

            llh=[]
            itrfdata={itrf:[[],[]] for itrf in itrfs}
            ndata=2+len(itrfs)*6

            for line in f:
                try:
                    ldata=[float(x) for x in line.split()]
                    if len(ldata) == 0:
                        continue
                    if len(ldata) != ndata:
                        raise RuntimeError("")
                except:
                    raise RuntimeError("Invalid data in "+filename+": "+line)

                llh.append(ldata[:2])
                nf=2
                for itrf in itrfs:
                    for y in (0,1):
                        itrfdata[itrf][y].append(ldata[nf:nf+3])
                        nf += 3

        for itrf in itrfs:
            for y in (0,1):
                itrfdata[itrf][y]=np.array(itrfdata[itrf][y])

        self.vec=np.zeros((len(llh),3))
        self.pts=np.array(llh)
        self.itrfdata=itrfdata
            
    def calcVectors( self, itrffrom, itrfto, year ):
        xyzf=self.itrfdata[itrffrom]
        xyzt=self.itrfdata[itrfto]
        f0=(self.year1-year)/(self.year1-self.year0)
        f1=1.0-f0
        self.itrffrom=itrffrom
        self.itrfto=itrfto
        self.vec=(xyzt[0]*f0+xyzt[1]*f1)-(xyzf[0]*f0+xyzf[1]*f1)
        self.setVectors( self.vec )


class ITRFWidget( QWidget ):

    def __init__( self, parent=None ):
        QWidget.__init__(self,parent)
        self.initGui()

    def initGui( self ):
        world = QtWorld()
        veclist=ITRFVectors()
        world.addVectors(veclist)

        grid=QGridLayout()
        itrfs=ITRFVectors.itrfs

        fslider=QSlider(Qt.Horizontal,self)
        flabel=QLabel()
        grid.addWidget(QLabel('ITRF from'),0,0)
        grid.addWidget(fslider,0,1)
        grid.addWidget(flabel,0,2)

        tslider=QSlider(Qt.Horizontal,self)
        tlabel=QLabel()
        grid.addWidget(QLabel('ITRF to'),1,0)
        grid.addWidget(tslider,1,1)
        grid.addWidget(tlabel,1,2)

        yslider=QSlider(Qt.Horizontal,self)
        ylabel=QLabel()
        grid.addWidget(QLabel('Epoch'),2,0)
        grid.addWidget(yslider,2,1)
        grid.addWidget(ylabel,2,2)

        sslider=QSlider(Qt.Horizontal,self)
        slabel=QLabel()
        grid.addWidget(QLabel('Scale'),3,0)
        grid.addWidget(sslider,3,1)
        grid.addWidget(slabel,3,2)

        layout=QVBoxLayout()
        layout.addWidget(world,1.0)
        layout.addWidget(world.controls())
        layout.addLayout(grid)

        self.setLayout( layout )

        #------------------------------

        def redraw():
            try:
                itrffrom=str(flabel.text())
                itrfto=str(tlabel.text())
                year=float(ylabel.text())
                veclist.calcVectors(itrffrom,itrfto,year)
            except:
                pass

        fslider.setMinimum(0)
        fslider.setMaximum(len(itrfs)-1)
        fslider.valueChanged.connect( lambda x: flabel.setText(itrfs[x]) is None and redraw())
        flabel.setText(itrfs[0])
        fslider.setValue(0)

        tslider.setMinimum(0)
        tslider.setMaximum(len(itrfs)-1)
        tslider.valueChanged.connect( lambda x: tlabel.setText(itrfs[x]) is None and redraw())
        tslider.setValue(len(itrfs)-1)

        year0=ITRFVectors.year0
        year1=ITRFVectors.year1
        nint=int((year1-year0)*5)
        yint=float(year1-year0)/nint

        yslider.setMinimum(0)
        yslider.setMaximum(nint)
        yslider.valueChanged.connect( lambda x: ylabel.setText("{0:.1f}".format(year0+x*yint)) is None and redraw())
        yslider.setValue(nint)
        
        sclmin=0.01
        sclmax=1000.0
        nscl=25
        lsclmin=log(sclmin)
        lsclmax=log(sclmax)
        lsclint=(lsclmax-lsclmin)/(nscl)

        sslider.setMinimum(0)
        sslider.setMaximum(nscl)
        sslider.valueChanged.connect( lambda x: slabel.setText("{0:.2f}".format(exp(lsclmin+x*lsclint))) is None and veclist.setScale(float(slabel.text())))
        sslider.setValue(15)


if __name__=='__main__':
    app = QApplication(sys.argv)
    window = ITRFWidget()
    window.show()
    sys.exit(app.exec_())
    main()
