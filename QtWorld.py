#!/usr/bin/python

import sys
import os.path
import math
import numpy as np
from PyQt4.QtGui import *
from PyQt4.QtCore import *

class QtWorld(QWidget):


    world_opacity=220
    
    def __init__(self):
        super(QtWorld, self).__init__()
        angle=np.concatenate(([0.0],np.radians(np.arange(0.0,180.0,5.0))))
        merid0=np.vstack((np.sin(angle),0.0*angle,np.cos(angle)))
        merid0[2,0]=1.2
        merid0=merid0.transpose()
        self.merid0 = merid0
        self.rot = np.identity(3)
        self.lon0 = 0
        self.radius=1.0
        self.mapwidth=5.0  # Map width as multiple of radius
        self.globe=False
        self.plotDarkSide=False
        self.vectors=[]
        self.loadCoast()
        self.initUI()
        
    def initUI(self):      
        self.setGeometry(300, 300, 600, 600)
        self.setMinimumSize(QSize(150,150))
        self.setWindowTitle('World plot')
        self.tracking=False
        self.lastPos=[0,0]
        self.cx=300
        self.cy=300
        self.show()

    def controls(self):
        control=QWidget()
        layout=QHBoxLayout()
        mapButton=QRadioButton("Map")
        mapButton.setChecked(not self.globe)
        globeButton=QRadioButton("Globe")
        globeButton.setChecked(self.globe)
        darkButton=QCheckBox("Plot dark side")
        darkButton.setChecked(self.plotDarkSide)
        layout.addWidget(mapButton)
        layout.addWidget(globeButton)
        layout.addWidget(darkButton)
        layout.addWidget(darkButton)
        layout.addItem(QSpacerItem(0,0,QSizePolicy.Ignored,QSizePolicy.Fixed))
        def globeClicked():
            self.globe=True
            darkButton.setEnabled(True)
            self.repaint()
        def mapClicked():
            self.globe=False
            darkButton.setEnabled(False)
            self.repaint()
        def darkSideChanged():
            self.plotDarkSide=darkButton.isChecked()
            self.repaint()
        globeButton.clicked.connect(globeClicked)
        mapButton.clicked.connect(mapClicked)
        darkButton.clicked.connect(darkSideChanged)
        control.setLayout(layout)
        return control

    def setGlobeMode( self, globe=True ):
        self.globe=globe
        self.repaint()

    def setMapMode( self, map=True ):
        self.setGlobeMode( not map )

    def mousePressEvent( self, event ):
        self.tracking=True
        self.lastPos=event.pos()

    def mouseMoveEvent( self, event ):
        if not self.tracking:
            return
        minshift=20
        rotfactor=40
        pos=event.pos()
        offset=pos-self.lastPos
        if not self.globe:
            offset.setY(0)
        if abs(offset.x())+abs(offset.y()) < minshift: 
            return

        lp=self.lastPos
        self.lastPos=pos

        if self.globe:
            self.rotateGlobe( lp, offset )
        else:
            self.rotateMap( lp, offset )
        self.repaint()

    def rotateGlobe( self, pos, offset ):
        dy=offset.x()
        dz=-offset.y()

        py=pos.x()-self.cx
        pz=-(pos.y()-self.cy)
        pr=math.sqrt(py*py+pz*pz)
        rf=pr/self.radius

        rotx=0
        if rf > 0.01:
            rf=(rf-0.01)**2
            if rf > 1.0:
                rf=1.0
            rotx=rf*(dz*py-dy*pz)/pr
            dy += rotx*pz/pr
            dz -= rotx*py/pr
            rotx /= pr

        crx=math.cos(rotx)
        srx=math.sin(rotx)
        rx=np.array([[1.0,0.0,0.0],[0.0,crx,srx],[0.0,-srx,crx]])

        r=math.sqrt(dy*dy+dz*dz)

        s1 = dy/r
        c1 = -dz/r
        r1=np.array([[1.0,0.0,0.0],[0.0,c1,s1],[0.0,-s1,c1]])

        c2=math.cos(r/self.radius)
        s2=math.sin(r/self.radius)
        r2=np.array([[c2,0.0,-s2],[0.0,1.0,0.0],[s2,0.0,c2]])

        rot=rx.dot(r1.transpose().dot(r2.dot(r1)))
        self.rot = self.rot.dot(rot)

    def rotateMap( self, pos, offset ):
        rot=360*offset.x()/(self.mapwidth*self.radius)
        self.lon0 -= rot

    def mouseReleaseEvent( self, event ):
        self.tracking=False


    def paintEvent(self, event):
        if self.globe:
            self.paintWorld()
        else:
            self.paintMap()

    def paintWorld(self):
        w=self.geometry().width()
        h=self.geometry().height()
	r=(w if w < h else h)/3
	cx=w/2
	cy=h/2
        self.radius=r
        self.cx=cx
        self.cy=cy

        offset=np.array([[0.0,cx,cy]])
        scale=np.array([[r,r,-r]])
        scalept=lambda p: p[:,0:3].dot(self.rot)*scale+offset
        scalevec=lambda p, x: x.dot(self.rot)*scale

        coastpts=scalept(self.cxyz)

        qp = QPainter()
        qp.begin(self)
        qp.setPen(QPen(QColor(0,0,0)))
        qp.setBrush(QBrush(QColor(255,255,255)))
        qp.drawEllipse(QRectF(cx-r, cy-r, r*2, r*2 ))

        qp.setPen(QPen(QColor(0,0,255)))

        if self.plotDarkSide:
            self.plotLines(qp,coastpts,breaks=self.cbreak,back=True)
            for v in self.vectors:
                v.paint(qp,scalept,scalevec,back=True)
            
        qp.setPen(QPen(QColor(0,0,0)))
        qp.setBrush(QBrush(QColor(255,255,255,self.world_opacity)))
        qp.drawEllipse(QRectF(cx-r, cy-r, r*2, r*2 ))

        qp.setPen(QPen(QColor(0,0,255)))
        self.plotLines(qp,coastpts,breaks=self.cbreak,back=False)
        for v in self.vectors:
            v.paint(qp,scalept,scalevec,back=False)
    
        self.paintScales(qp)
        qp.end()

    def paintMap(self):
        w=self.geometry().width()
        h=self.geometry().height()
        r1=w/(self.mapwidth+1.0)
        r2=h/3.0
        r=min(r1,r2)
	cx=w/2
	cy=h/2

        self.radius=r
        self.cx=cx
        self.cy=cy

        offset=np.array([[0.0,cx,cy]])
        
        def scalept(p):
            return np.vstack((
                p[:,0]*0.0,
                (np.mod(p[:,3]-self.lon0+180,360)-180)*(self.mapwidth/360.0),
                -p[:,2]
                )).transpose()*r+offset

        scale=np.array([[-r,r,-r]])

        def scalevec(p,x):
            evec=np.vstack((-np.sin(np.radians(p[:,3])),np.cos(np.radians(p[:,3])),p[:,3]*0.0)).transpose()
            nvec=np.cross(p[:,0:3],evec)
            # Oth vectors...
            uen=np.vstack((
                np.sum(p[:,0:3]*x,axis=1),
                np.sum(evec*x,axis=1),
                np.sum(nvec*x,axis=1),
            )).transpose()*scale
            return uen


        coastpts=scalept(self.cxyz)
        cbreak=np.logical_or(self.cbreak,False)
        cbreak[:-1]=np.logical_or(cbreak[:-1],np.abs(coastpts[:-1,1]-coastpts[1:,1]) > r*2)

        qp=QPainter()
        qp.begin(self)
        qp.setPen(QPen(QColor(0,0,0)))
        qp.setBrush(QBrush(QColor(255,255,255)))
        qp.drawRect(QRectF(cx-r*self.mapwidth/2.0, cy-r, r*self.mapwidth, r*2 ))

        qp.setPen(QPen(QColor(0,0,255)))
        self.plotLines(qp,coastpts,breaks=cbreak)
        for v in self.vectors:
            v.paint(qp,scalept,scalevec,hv=True)
    
        self.paintScales(qp)
        qp.end()

    def paintScales( self, painter ):
        fm=painter.fontMetrics()
        lh=fm.height()
        px=self.width()-fm.width('m');
        py=self.height()-lh/2.0;
        for v in self.vectors:
            v.paintScale(painter,px,py,self.width()/10.0,self.radius)

    def plotLines( self, painter, pts, breaks=None, back=None ):
        if breaks is None:
            breaks=points[:,0] != points[:,0]
        if back is not None:
            if back:
                breaks=np.logical_or(breaks,pts[:,0]>0)
            else:
                breaks=np.logical_or(breaks,pts[:,0]<0)
        
        for p0,p1,brk in zip(pts[:-1],pts[1:],breaks[:-1]):
            if not brk:
                painter.drawLine(p0[1],p0[2],p1[1],p1[2])

    @staticmethod
    def pointsFromLatLon( lonlat ):
        lon=np.radians(lonlat[:,0])
        lat=np.radians(lonlat[:,1])
        return np.vstack((np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat),lonlat[:,0])).transpose()

    def loadCoast(self):
        cstf=os.path.join(os.path.dirname(__file__),'coasts.dat')
        cpoints=np.loadtxt(cstf)
        cbreak=cpoints[:,0] < -360
        cbreak[:-1] = np.logical_or(cbreak[:-1],cbreak[1:])
        self.cxyz=QtWorld.pointsFromLatLon( cpoints )
        self.cbreak=cbreak

    def addVectors( self, vectorset ):
        if vectorset in self.vectors:
            return
        vectorset.changed.connect( self.repaint )
        self.vectors.append(vectorset)

class VectorSet( QObject ):

    changed=pyqtSignal( name='changed' )

    def __init__( self, points, vectors ):
        '''
        Expects an array of points [[lon1,lat1],[lon2,lat2]...]
        and an array of vectors [[dx1,dy1,dz1],[dx2,dy2,dz2]...]
        '''
        QObject.__init__( self )
        self.pts=QtWorld.pointsFromLatLon(points)
        self.vec=vectors
        self.scale=1.0
        self.colour=QColor('red')
        self.vcolour=QColor('blue')
        self.basesize=5.0
        self.linewidth=2.0
        self.units='m'
        self.unitFactor=1.0

    def setVectors( self, vectors ):
        self.vec=vectors
        self.changed.emit()

    def setStyle( self, colour=None, vcolour=None, basesize=None, linewidth=None ):
        self.colour=colour or self.colour
        self.vcolour=vcolour or self.vcolour
        self.basesize=basesize or self.basesize
        self.linewidth=linewidth or self.linewidth
        self.changed.emit()

    def setScale( self, scale ):
        self.scale = scale
        self.changed.emit()

    def setScaleUnits( self, units, unitFactor=1.0 ):
        self.units=units
        self.unitFactor=unitFactor

    def paint( self, painter, scalept, scalevec, back=None, hv=False ):
        pts=scalept(self.pts)
        vec=scalevec(self.pts,self.vec)*self.scale
        pen0=painter.pen()
        brush0=painter.brush()
        painter.setPen(QPen(self.colour))
        painter.setBrush(QBrush(self.colour))
        r=self.basesize/2.0
        for p in pts:
            if back is not None and (p[0] < 0) ^ back:
                continue
            painter.drawEllipse(QRectF(p[1]-r, p[2]-r, r*2, r*2 ))
        painter.setPen(QPen(self.colour,self.linewidth))
        for p0,p1 in zip(pts,vec):
            if back is not None and (p0[0] < 0) ^ back:
                continue
            painter.drawLine(p0[1],p0[2],p0[1]+p1[1],p0[2]+p1[2])
        if hv:
            painter.setPen(QPen(self.vcolour,self.linewidth))
            for p0,p1 in zip(pts,vec):
                if back is not None and (p0[0] < 0) ^ back:
                    continue
                painter.drawLine(p0[1],p0[2],p0[1],p0[2]+p1[0])
        painter.setPen(pen0)
        painter.setBrush(brush0)

    def paintScale( self, painter, xr, yr, width, r ):
        scllen=abs(width/(r*self.scale*self.unitFactor))
        ndp=int(math.floor(math.log10(scllen)))
        scllen1=10**ndp
        scllen2=scllen1
        ndp2=-ndp
        for x,ndpx in zip((1.5,2,3,5,10),(1,0,0,0,-1)):
            if x*scllen1 > scllen:
                break
            scllen2=x*scllen1
            ndp2=(-ndp)+ndpx
        ndp2=max(0,ndp2)
        plotlen=scllen2*r*self.scale*self.unitFactor
        sclstr=" {0:.{ndp}f}{1}".format(scllen2,self.units,ndp=ndp2)
        txtflags=Qt.TextSingleLine | Qt.TextDontClip | Qt.AlignRight | Qt.AlignBottom
        txtrect0=QRectF(xr,yr,0.0,0.0)
        txtrect=painter.boundingRect(txtrect0,txtflags,sclstr)
        painter.drawText(txtrect0,txtflags,sclstr)

        sclx=txtrect.left()-plotlen
        scly=txtrect.center().y()

        pen0=painter.pen()
        brush0=painter.brush()
        painter.setPen(QPen(self.colour))
        painter.setBrush(QBrush(self.colour))
        r=self.basesize/2.0
        painter.drawEllipse(QRectF(sclx-r, scly-r, r*2, r*2 ))
        painter.setPen(QPen(self.colour,self.linewidth))
        painter.drawLine(sclx,scly,sclx+plotlen,scly)
        painter.setPen(pen0)
        painter.setBrush(brush0)
        return txtrect.top()
        

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = QtWorld()
    sys.exit(app.exec_())
    main()
