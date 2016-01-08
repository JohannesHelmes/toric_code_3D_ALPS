#!/usr/bin/env python

import sys
import numpy as np
from PyQt4 import QtGui, QtCore
from IncMakerGuiTC3D import Ui_Dialog as Dlg

class MeinDialog(QtGui.QDialog, Dlg): 
    L=0
    Coords=[]
    Lines=[]
    Geometry=[]
    CurrIncStep=0
    MaxIncStep=0
    CurrPlane=0
    framesize=760
    diam=7
    RecentlyTouched=[]
    changed=False
    File=None
    def __init__(self): 
        QtGui.QDialog.__init__(self) 
        self.setupUi(self)
        self.Frame.setParent(None)
        self.setMouseTracking(True)
        self.RecentlyTouched={co: 0 for co in self.Coords}
        self.connect(self.Quit,QtCore.SIGNAL("clicked()"),self.onQuit)
        self.connect(self.Create,QtCore.SIGNAL("clicked()"),self.onCreate)
        self.connect(self.Change21,QtCore.SIGNAL("clicked()"),self.onChange21)
        self.connect(self.Save,QtCore.SIGNAL("clicked()"),self.onSave)
        self.connect(self.Load,QtCore.SIGNAL("clicked()"),self.onLoad)
        self.connect(self.CopyIncStep,QtCore.SIGNAL("clicked()"),self.onCopyToNextIncStep)
        self.connect(self.CopyPlane,QtCore.SIGNAL("clicked()"),self.onCopyToNextPlane)
        self.connect(self.IncStepNo,QtCore.SIGNAL("valueChanged(int)"),self.NoLChange)
        self.connect(self.PlaneNo,QtCore.SIGNAL("valueChanged(int)"),self.PlaneChange)
        self.Choose.clicked.connect(self.onChoose)
        self.pen = QtGui.QPen(QtGui.QColor(0,0,0))
        self.pen.setWidth(2)
        self.brush = QtGui.QBrush(QtGui.QColor(255,255,255))
        #self.Lines.append("")
        self.IncStepNo.setRange(0,0)


    def onQuit(self):
        self.close()

    def onChoose(self):
        self.Filename.setText(QtGui.QFileDialog.getSaveFileName())

    def PlaneChange(self):
        self.CurrPlane=self.PlaneNo.value()
        self.changed=True
        self.update()
        self.setFocus()

    def NoLChange(self):
        self.CurrIncStep=self.IncStepNo.value()
        self.changed=True
        self.update()
        self.setFocus()

    def onCopyToNextPlane(self):
        if self.CurrPlane==self.L-1:
            return
        incstep=self.CurrIncStep
        k=self.CurrPlane
        self.Geometry[incstep][k+1][:][:][:]=self.Geometry[incstep][k][:][:][:]
        self.PlaneNo.setValue(self.CurrPlane+1)
        self.PlaneChange()


    def onCopyToNextIncStep(self):
        if self.CurrIncStep==self.MaxIncStep:
            self.MaxIncStep+=1
            self.IncStepNo.setRange(0,self.MaxIncStep)
            self.Geometry.append(np.zeros((self.L,self.L,self.L,3),dtype=np.int8))
        self.Geometry[self.CurrIncStep+1]=np.copy(self.Geometry[self.CurrIncStep])
        self.IncStepNo.setValue(self.CurrIncStep+1)
        self.NoLChange()
            

    def onSave(self):
        self.File=open(self.Filename.text(),'w')
        for s in range(self.MaxIncStep+1):
            for k in range(self.L):
                for j in range(self.L):
                    for i in range(self.L):
                        for l in range(3):
                            sys.stdout.write(str(self.Geometry[s][k][j][i][l]))
                            self.File.write(str(self.Geometry[s][k][j][i][l]))
            sys.stdout.write('\n')
            self.File.write('\n')
        self.File.close()

    def onLoad(self):
        self.File=open(self.Filename.text(),'r')
        content=self.File.readlines()
        
        self.LValue.setValue(int(round((len(content[0])/3)**(1./3))))
        self.L=self.LValue.value()
        self.CurrIncStep=0
        self.MaxIncStep=len(content)-1
        self.IncStepNo.setRange(0,self.MaxIncStep)
        self.PlaneNo.setRange(0,self.L)
        self.onCreate()
        self.Geometry=[]
        for s,line in enumerate(content):
            self.Geometry.append(np.zeros((self.L,self.L,self.L,3),dtype=np.int8)) #4th coordinate: (0) right, (1) down, (2) back/diagonal
            for k in range(self.L):
                for j in range(self.L):
                    for i in range(self.L):
                        for l in range(3):
                            self.Geometry[s][k][j][i][l]=line[3*k*self.L**2+3*j*self.L+i*3+l]

        self.File.close()
        self.NoLChange()

       


    def onChange21(self):
        incstep=self.CurrIncStep
        k=self.CurrPlane
        for j in range(self.L):
            for i in range(self.L):
                for l in range(3):
                    if self.Geometry[incstep][k][j][i][l]==2:
                        self.Geometry[incstep][k][j][i][l]=1
        self.changed=True
        self.update()


    def Change01(self):
        incstep=self.CurrIncStep
        k=self.CurrPlane
        for j in range(self.L):
            for i in range(self.L):
                for l in range(3):
                    if self.Geometry[incstep][k][j][i][l]==0:
                        self.Geometry[incstep][k][j][i][l]=1
                    elif self.Geometry[incstep][k][j][i][l]==1:
                        self.Geometry[incstep][k][j][i][l]=0
        self.changed=True
        self.update()

    def keyPressEvent(self,event):
        if event.key()==QtCore.Qt.Key_A:
            self.onChange21()
        if event.key()==QtCore.Qt.Key_P:
            self.onCopyToNextPlane()
        if event.key()==QtCore.Qt.Key_I:
            self.onCopyToNextIncStep()
        if event.key()==QtCore.Qt.Key_D:
            self.Change01()

    def onCreate(self):
        self.L=self.LValue.value()
        xmin=10
        ymin=10
        xsize=self.framesize
        ysize=self.framesize
        plaqsize=self.framesize/self.L
        self.Geometry=[np.zeros((self.L,self.L,self.L,3),dtype=np.int8)]
        self.Coords=[]
        self.RecentlyTouched=np.zeros((self.L,self.L,3))
        self.PlaneNo.setRange(0,self.L-1)
        for j,y in enumerate(np.linspace(0,ysize,self.L+1)[:-1]):
            self.Coords.append([])
            for i,x in enumerate(np.linspace(0,xsize,self.L+1)[:-1]):
                self.Coords[j].append([])
                self.Coords[j][i].append(QtCore.QPoint(xmin+x+plaqsize/2,ymin+y))
                self.Coords[j][i].append(QtCore.QPoint(xmin+x,ymin+y+plaqsize/2))
                self.Coords[j][i].append(QtCore.QPoint(xmin+x+plaqsize/2,ymin+y+plaqsize/2))
                self.RecentlyTouched[j][i][0]=0
                self.RecentlyTouched[j][i][1]=0
                self.RecentlyTouched[j][i][2]=0
        self.changed=True
        self.update()

    def paintEvent(self,event):
        if (self.L>0)or(self.changed):
            self.changed=False

            painter = QtGui.QPainter(self)
            xmin=10
            ymin=10
            xsize=self.framesize
            ysize=self.framesize
            plaqsize=self.framesize/self.L
            incstep=self.CurrIncStep
            for x in np.linspace(0,xsize,self.L+1)[:-1]:
                for y in np.linspace(0,ysize,self.L+1)[:-1]:
                    self.brush = QtGui.QBrush(QtGui.QColor(255,255,255))
                    painter.setPen(self.pen)
                    painter.setBrush(self.brush)
                    painter.drawRect(xmin+x,ymin+y,plaqsize,plaqsize)
                    painter.setPen(3)
                    painter.drawLine(xmin+x,ymin+y,xmin+x+plaqsize/2,ymin+y+plaqsize/2)
                    painter.setPen(1)
            k=self.CurrPlane
            for j in range(self.L):
                for i in range(self.L):
                    for l in range(3):
                        state=self.Geometry[incstep][k][j][i][l]
                        coords=self.Coords[j][i][l]
                        if state==1:
                            self.brush = QtGui.QBrush(QtGui.QColor(0,0,255))
                        elif state==2:
                            self.brush = QtGui.QBrush(QtGui.QColor(255,0,0))
                        else:
                            self.brush = QtGui.QBrush(QtGui.QColor(255,255,255))
                        painter.setBrush(self.brush)
                        painter.drawEllipse(coords,self.diam,self.diam)


    beginCoordinates=[]
    def mousePressEvent(self,event):
        self.beginCoordinates=event.pos()
    
    def mouseReleaseEvent(self,event):
        endCoordinates=event.pos()
        for j in range(self.L):
            for i in range(self.L):
                for l in range(3):
                    spin=self.Coords[j][i][l]
                    distB=spin-self.beginCoordinates
                    distE=spin-endCoordinates
                    if (distB.x()**2+distB.y()**2<self.diam**2):
                        Bi=i
                        Bj=j
                    if (distE.x()**2+distE.y()**2<self.diam**2):
                        Ei=i
                        Ej=j

        incstep=self.CurrIncStep
        k=self.CurrPlane
        for j in range(Bj,Ej+1):
            for i in range(Bi,Ei+1):
                for l in range(3):
                    self.Geometry[incstep][k][j][i][l]=2
        self.changed=True
        self.update()


    def mouseMoveEvent(self,event):
        if (event.buttons()==QtCore.Qt.NoButton):
            k=self.CurrPlane
            incstep=self.CurrIncStep
            for j in range(self.L):
                for i in range(self.L):
                    for l in range(3):
                        spin=self.Coords[j][i][l]
                        dist=spin-event.pos()
                        if (dist.x()**2+dist.y()**2<self.diam**2):
                            if (self.RecentlyTouched[j][i][l]!=0):
                                return
                            else:
                                self.Geometry[incstep][k][j][i][l]+=2
                                self.Geometry[incstep][k][j][i][l]%=3
                                self.RecentlyTouched[j][i][l]=1
                                self.changed=True
                                self.update()
                                return
                        self.RecentlyTouched[j][i][l]=0

app = QtGui.QApplication(sys.argv) 
dialog = MeinDialog() 
dialog.show() 
sys.exit(app.exec_()) 
