#!/usr/bin/env python

import sys
import numpy as np
from PyQt4 import QtGui, QtCore
from IncMakerGuiTC2D import Ui_Dialog as Dlg

class MeinDialog(QtGui.QDialog, Dlg): 
    L=0
    Coords=[]
    Lines=[]
    NoL=0
    State={}
    RecentlyTouched={}
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
        self.connect(self.Next,QtCore.SIGNAL("clicked()"),self.onNext)
        self.connect(self.ILValue,QtCore.SIGNAL("valueChanged(int)"),self.NoLChange)
        self.Choose.clicked.connect(self.onChoose)
        self.pen = QtGui.QPen(QtGui.QColor(0,0,0))
        self.pen.setWidth(2)
        self.brush = QtGui.QBrush(QtGui.QColor(255,255,255))
        #self.Lines.append("")
        self.ILValue.setRange(0,0)


    def onQuit(self):
        self.close()

    def onChoose(self):
        self.Filename.setText(QtGui.QFileDialog.getSaveFileName())

    def NoLChange(self):
        self.NoL=self.ILValue.value()
        for char,coords in zip(self.Lines[self.NoL][:-1],self.Coords):
            self.State[coords]=int(char)

        self.changed=True
        self.update()
        self.setFocus()

    def onNext(self):
        if len(self.Lines)==0:
            self.Lines.append("")
        self.Lines[self.NoL]=""
        for coords in self.Coords:
            self.Lines[self.NoL]=self.Lines[self.NoL]+str(self.State[coords])
        self.Lines[self.NoL]=self.Lines[self.NoL]+'\n'
        self.Lines.append("")
        self.NoL+=1
        self.ILValue.setRange(0,self.NoL)
        self.ILValue.setValue(self.NoL)

    def onSave(self):
        self.File=open(self.Filename.text(),'w')
        for line in self.Lines:
            sys.stdout.write(line)
            self.File.write(line)
        self.File.close()

    def onLoad(self):
        self.File=open(self.Filename.text(),'r')
        content=self.File.readlines()
        self.LValue.setValue(int((len(content[0])/2)**0.5))
        self.NoL=0
        self.ILValue.setRange(0,len(content))
        for i,line in enumerate(content):
            self.Lines.append("")
            for char in line:
                self.Lines[i]=self.Lines[i]+char

        self.Lines.append("")
        for char in content[-1]:
            self.Lines[-1]=self.Lines[-1]+char

        self.File.close()
        self.onCreate()
        self.NoLChange()

       


    def onChange21(self):
        for coords in self.Coords:
            if self.State[coords]==2:
                self.State[coords]=1
        self.changed=True
        self.update()


    def Change01(self):
        for coords in self.Coords:
            if self.State[coords]==1:
                self.State[coords]=0
            elif self.State[coords]==0:
                self.State[coords]=1

        self.changed=True
        self.update()

    def keyPressEvent(self,event):
        if event.key()==QtCore.Qt.Key_A:
            self.onChange21()
        if event.key()==QtCore.Qt.Key_N:
            self.onNext()
        if event.key()==QtCore.Qt.Key_D:
            self.Change01()

    def onCreate(self):
        self.L=self.LValue.value()
        xmin=30
        ymin=120
        xsize=580
        ysize=580
        plaqsize=580/self.L
        self.Coords=[]
        for x in np.linspace(0,xsize,self.L+1)[:-1]:
            for y in np.linspace(0,ysize,self.L+1)[:-1]:
                self.Coords.append(QtCore.QPoint(xmin+x+plaqsize/2,ymin+y))
                self.Coords.append(QtCore.QPoint(xmin+x,ymin+y+plaqsize/2))
        for coords in self.Coords:
            self.State[coords]=0
        self.changed=True
        self.update()

    def paintEvent(self,event):
        if (self.L>0)or(self.changed):
            self.changed=False

            painter = QtGui.QPainter(self)
            xmin=30
            ymin=120
            xsize=580
            ysize=580
            plaqsize=580/self.L
            for x in np.linspace(0,xsize,self.L+1)[:-1]:
                for y in np.linspace(0,ysize,self.L+1)[:-1]:
                    self.brush = QtGui.QBrush(QtGui.QColor(255,255,255))
                    painter.setPen(self.pen)
                    painter.setBrush(self.brush)
                    painter.drawRect(xmin+x,ymin+y,plaqsize,plaqsize)
            for coords in self.Coords:
                    if self.State[coords]==1:
                        self.brush = QtGui.QBrush(QtGui.QColor(0,0,255))
                    elif self.State[coords]==2:
                        self.brush = QtGui.QBrush(QtGui.QColor(255,0,0))
                    else:
                        self.brush = QtGui.QBrush(QtGui.QColor(255,255,255))
                    painter.setBrush(self.brush)
                    painter.drawEllipse(coords,10,10)

    def mouseMoveEvent(self,event):
        for spin in self.Coords:
            dist=spin-event.pos()
            if (dist.x()**2+dist.y()**2<100):
                if (self.RecentlyTouched[spin]!=0):
                    return
                else:
                    self.State[spin]+=2
                    self.State[spin]%=3
                    self.RecentlyTouched[spin]=1
                    self.changed=True
                    self.update()
                    return
        self.RecentlyTouched={co: 0 for co in self.Coords}

app = QtGui.QApplication(sys.argv) 
dialog = MeinDialog() 
dialog.show() 
sys.exit(app.exec_()) 
