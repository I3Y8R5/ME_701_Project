from PyQt5.QtWidgets import (QApplication, QMainWindow, QLineEdit, QDialog, QLabel, 
                             QVBoxLayout, QHBoxLayout, QAction, QMessageBox, QComboBox, QFileDialog,
                             QPushButton, QWidget, QTabWidget)
from PyQt5.QtCore import QT_VERSION_STR, PYQT_VERSION_STR

import platform
import numpy as np
import sys

class MainWindow(QMainWindow):
    
    def __init__(self):
        QMainWindow.__init__(self)
        
        # Create the File menu
        self.menuFile = self.menuBar().addMenu("&File")
        self.actionLoad = QAction("&Load", self)
        self.actionLoad.triggered.connect(self.load)
        self.actionSaveAs = QAction("&Save As", self)
        self.actionSaveAs.triggered.connect(self.saveas)
        self.actionSaveTmp = QAction("&Save Template", self)
        self.actionSaveTmp.triggered.connect(self.SaveTmp)
        self.actionQuit = QAction("&Quit", self)
        self.actionQuit.triggered.connect(self.close)
        
        self.menuFile.addActions([self.actionLoad, self.actionSaveAs, self.actionSaveTmp, self.actionQuit])
    
        
        
        # Create the Help menu
        self.menuHelp = self.menuBar().addMenu("&Help")
        self.actionAbout = QAction("&About",self)
        self.actionAbout.triggered.connect(self.about)
        self.menuHelp.addActions([self.actionAbout])
        #______________________________________________________________
        
        
    
        # Initialize Tabs
        self.tab = QTabWidget()
        self.tab0 = QWidget()
        self.tab1 = QWidget()
        self.tab2 = QWidget()
        self.tab3 = QWidget()
        self.tab4 = QWidget()
        #______________________________________________________________
        
        # Setup main widget
        widget = QDialog()
        
        # Trying to set up tabs
        self.tab.addTab(self.tab0, "Input")
        self.tab.addTab(self.tab1, "Output")
        self.tab.addTab(self.tab2, "X-Axis")
        self.tab.addTab(self.tab3, "Y-Axis")
        self.tab.addTab(self.tab4, "Z-Axis")
        #______________________________________________________________
        
        layout = QVBoxLayout()
        layout.addWidget(self.tab)        

        # Creating tab0
        self.tab0.layout = QHBoxLayout()  
        
        # Creating Collumn One
        self.ocoll1 = QWidget()
        self.ocoll1.layout = QVBoxLayout()
        self.ol11 = QLabel("E Matrix")
        self.oc11 = QLineEdit()
        self.ol12 = QLabel("G Matrix")
        self.oc12 = QLineEdit()
        self.ol13 = QLabel("Volume Fraction of Matrix")
        self.oc13 = QLineEdit()
        self.ocoll1.layout.addWidget(self.ol11)
        self.ocoll1.layout.addWidget(self.oc11)
        self.ocoll1.layout.addWidget(self.ol12)
        self.ocoll1.layout.addWidget(self.oc12)
        self.ocoll1.layout.addWidget(self.ol13)
        self.ocoll1.layout.addWidget(self.oc13)
        self.ocoll1.setLayout(self.ocoll1.layout)

        # Creating Collumn Two
        self.ocoll2 = QWidget()
        self.ocoll2.layout = QVBoxLayout()
        self.ol21 = QLabel("E Fiber")
        self.oc21 = QLineEdit()
        self.ol22 = QLabel("G Fiber")
        self.oc22 = QLineEdit()
        self.ol23 = QLabel("Volume Fraction of Fiber")
        self.oc23 = QLineEdit()
        self.ocoll2.layout.addWidget(self.ol21)
        self.ocoll2.layout.addWidget(self.oc21)
        self.ocoll2.layout.addWidget(self.ol22)
        self.ocoll2.layout.addWidget(self.oc22)
        self.ocoll2.layout.addWidget(self.ol23)
        self.ocoll2.layout.addWidget(self.oc23)
        self.ocoll2.setLayout(self.ocoll2.layout)
        
        # Creating Collumn Three
        self.ocoll3 = QWidget()
        self.ocoll3.layout = QVBoxLayout()
        self.ol31 = QLabel("Poisson's Ratio of Matrix")
        self.oc31 = QLineEdit()
        self.ol32 = QLabel("Poisson's Ratio of Fiber")
        self.oc32 = QLineEdit()
        self.ol33 = QLabel("")
        self.oc33 = QPushButton("Calculate")
        self.oc33.pressed.connect(self.Calc)
        self.ocoll3.layout.addWidget(self.ol31)
        self.ocoll3.layout.addWidget(self.oc31)
        self.ocoll3.layout.addWidget(self.ol32)
        self.ocoll3.layout.addWidget(self.oc32)
        self.ocoll3.layout.addWidget(self.ol33)
        self.ocoll3.layout.addWidget(self.oc33)
        self.ocoll3.setLayout(self.ocoll3.layout)
        
        self.tab0.layout.addWidget(self.ocoll1)
        self.tab0.layout.addWidget(self.ocoll2)
        self.tab0.layout.addWidget(self.ocoll3)

        self.tab0.setLayout(self.tab0.layout)
        #_______________________________________________________
        
        # Creating tab1
        self.tab1.layout = QHBoxLayout()
        global dim
        dim = 3
        
        
        # Creating Collumn One
        self.coll1 = QWidget()
        self.coll1.layout = QVBoxLayout()
        self.l11 = QLabel("E 1-1")
        self.c11 = QLineEdit()
        self.l12 = QLabel("G 1-2")
        self.c12 = QLineEdit()
        self.l13 = QLabel("V 1-2")
        self.c13 = QLineEdit()
        self.calc = QPushButton('Calculate')
        self.calc.pressed.connect(self.Calc)
        self.coll1.layout.addWidget(self.l11)
        self.coll1.layout.addWidget(self.c11)
        self.coll1.layout.addWidget(self.l12)
        self.coll1.layout.addWidget(self.c12)
        self.coll1.layout.addWidget(self.l13)
        self.coll1.layout.addWidget(self.c13)
        self.coll1.layout.addWidget(self.calc)
        self.coll1.setLayout(self.coll1.layout)

        # Creating Collumn Two
        self.coll2 = QWidget()
        self.coll2.layout = QVBoxLayout()
        self.l21 = QLabel("E 2-2")
        self.c21 = QLineEdit()
        self.l22 = QLabel("G 2-3")
        self.c22 = QLineEdit()
        self.l23 = QLabel("V 1-3")
        self.c23 = QLineEdit()
        self.save = QPushButton('Save')
        self.save.pressed.connect(self.saveas)
        self.coll2.layout.addWidget(self.l21)
        self.coll2.layout.addWidget(self.c21)
        self.coll2.layout.addWidget(self.l22)
        self.coll2.layout.addWidget(self.c22)
        self.coll2.layout.addWidget(self.l23)
        self.coll2.layout.addWidget(self.c23)
        self.coll2.layout.addWidget(self.save)
        self.coll2.setLayout(self.coll2.layout)
        
        # Creating Collumn Three
        self.coll3 = QWidget()
        self.coll3.layout = QVBoxLayout()
        self.l31 = QLabel("E 2-3")
        self.c31 = QLineEdit()
        self.l32 = QLabel("G 1-3")
        self.c32 = QLineEdit()
        self.l33 = QLabel("V 2-3")
        self.c33 = QLineEdit()
        self.input = QPushButton('Input')
        self.coll3.layout.addWidget(self.l31)
        self.coll3.layout.addWidget(self.c31)
        self.coll3.layout.addWidget(self.l32)
        self.coll3.layout.addWidget(self.c32)
        self.coll3.layout.addWidget(self.l33)
        self.coll3.layout.addWidget(self.c33)
        self.coll3.layout.addWidget(self.input)
        self.coll3.setLayout(self.coll3.layout)
        
        # Creating Collumn Four
        self.coll4 = QWidget()
        self.coll4.layout = QVBoxLayout()
        self.plot1 = MatplotlibCanvas()
        self.coll4.layout.addWidget(self.plot1)
        self.coll4.setLayout(self.coll4.layout)        

        self.tab1.layout.addWidget(self.coll1)
        self.tab1.layout.addWidget(self.coll2)
        self.tab1.layout.addWidget(self.coll3)
        self.tab1.layout.addWidget(self.coll4)

        self.tab1.setLayout(self.tab1.layout)
        #_______________________________________________
        
        
        # Creating tab2
        global view
        view = 'x'
        self.tab2.layout = QHBoxLayout()
        dim = 2
        self.xplot = MatplotlibCanvas()
        
        self.spc1 = QLabel("                                                  ")
        self.spc2 = QLabel("                                                  ")
        self.tab2.layout.addWidget(self.spc1)
        self.tab2.layout.addWidget(self.xplot)
        self.tab2.layout.addWidget(self.spc2)
        self.tab2.setLayout(self.tab2.layout)
        #_______________________________________________
        
        # Creating tab3
        view = 'y'
        self.tab3.layout = QHBoxLayout()
        self.yplot = MatplotlibCanvas()
        self.spc3 = QLabel("                                                  ")
        self.spc4 = QLabel("                                                  ")
        self.tab3.layout.addWidget(self.spc3)
        self.tab3.layout.addWidget(self.yplot)
        self.tab3.layout.addWidget(self.spc4)
        self.tab3.setLayout(self.tab3.layout)
        #_______________________________________________
        
        # Creating tab4
        view = 'z'
        self.tab4.layout = QHBoxLayout()
        self.zplot = MatplotlibCanvas()
        self.spc5 = QLabel("                                                  ")
        self.spc6 = QLabel("                                                  ")
        self.tab4.layout.addWidget(self.spc5)
        self.tab4.layout.addWidget(self.zplot)
        self.tab4.layout.addWidget(self.spc6)
        self.tab4.setLayout(self.tab4.layout) 
        #_______________________________________________
        
        # Basically each tab is a widget in and of itself
        # For each tab you have to set the layout and add 
        # Widgets just like the main window
        

        widget.setLayout(layout)
        
        
        self.setCentralWidget(widget)
        
        
    def Calc(self):
        # Parameters Needed For Calculation
        Em = float(self.oc11.text())
        Ef = float(self.oc21.text())
        vm = float(self.oc13.text())
        vf = float(self.oc23.text())
        Gm = float(self.oc12.text())
        Gf = float(self.oc22.text())
        Vm = float(self.oc31.text())
        Vf = float(self.oc32.text())
        #____________________________________
        
        # Calculations
        E11 = (Em*vm + Ef*vf)
        E22 = ((vf/Ef) + (vm/Em))**(-1)
        E33 = E22
        eta = (3 -4*Vm + (Gm/Gf))/(4*(1 - Vm))
        G23 = (Gm*vf + eta*(1 - vf))/(eta*(1 - vf) + vf*(Gm/Gf))
        G12 = ((vf/Gf) + (vm/Gm))**(-1)
        G13 = G12
        Vxy = (vf*Vf + vm*Vm)
        Vzx = Vxy
        Vyz = Vm
        #_________________________________________________________
        self.c11.setText(str(E11))
        self.c21.setText(str(E22))
        self.c31.setText(str(E33))
        self.c12.setText(str(G12))
        self.c22.setText(str(G23))
        self.c32.setText(str(G13)) 
        self.c13.setText(str(Vxy))
        self.c23.setText(str(Vzx))
        self.c33.setText(str(Vyz))
        

    """
    Make modificaitions to saveas so that the material can be saved
    """

    def saveas(self) :
        options = QFileDialog.Options()
        fname, _ = QFileDialog.getSaveFileName(self, "Save as", options=options)
        #t0 = ['E 1-1', 'E 2-2', 'E 3-3']
        x0 = [eval(self.c11.text()),eval(self.c12.text()), eval(self.c13.text())]
        #t1 = ['G 1-2', 'G 2-3', 'G 1-3']
        x1 = [eval(self.c21.text()), eval(self.c22.text()), eval(self.c23.text())]
        #t2 = ['V 1-2', 'V 1-3', 'V 2-3']
        x2 = [eval(self.c31.text()), eval(self.c32.text()), eval(self.c33.text())]
        x = np.array((x0, x1, x2))
        x = x.transpose()
        np.savetxt(fname, x)
        
    """
    Use the above saveas to fill out the SaveTmp
    """
    
    def SaveTmp(self):
        options = QFileDialog.Options()
        fname, _ = QFileDialog.getSaveFileName(self, "Save Template", options=options)
        x0 = [eval(self.oc11.text()),eval(self.oc12.text()), eval(self.oc13.text())]
        x1 = [eval(self.oc21.text()), eval(self.oc22.text()), eval(self.oc23.text())]
        x2 = [eval(self.oc31.text()), eval(self.oc32.text()), 0.00]
        x = np.array((x0, x1, x2))
        np.savetxt(fname, x)
        
    
    def load(self, fopen):
        options = QFileDialog.Options()
        fname, _ = QFileDialog.getOpenFileName(self, "Load", options=options)
        i = 0
        with open(fname, 'r') as file:
            for line in file:
                x = line.split()
                print(i)
                print(x)
                if (i == 0):
                    self.oc11.setText(x[0])
                    self.oc12.setText(x[1])
                    self.oc13.setText(x[2])
                if (i == 1):
                    self.oc21.setText(x[0])
                    self.oc22.setText(x[1])
                    self.oc23.setText(x[2])
                if (i == 2):
                    self.oc31.setText(x[0])
                    self.oc32.setText(x[1])

                    
                i = i + 1
        fname.close()
                
    def about(self) :
        QMessageBox.about(self, 
            "About Composite Evaluator",
            """<b>Composite Evaluator</b>
               <p>Copyright &copy; 2017 Kord Byers, All Rights Reserved.
               <p>Python %s -- Qt %s -- PyQt %s on %s""" %
            (platform.python_version(),
             QT_VERSION_STR, PYQT_VERSION_STR, platform.system()))
        
from PyQt5.QtWidgets import QSizePolicy
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D

class MatplotlibCanvas(FigureCanvas):
    
    
    def __init__(self, parent=None, width=10, height=10, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        
        if (dim == 2):
            self.axes = self.fig.add_subplot(1,1,1)
            self.compute_initial_figure_2D()
            
        if (dim == 3):
            self.axes = self.fig.add_subplot(111, projection='3d')
            self.compute_initial_figure_3D()
            
        self.axes.hold(False)
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
    
    def compute_initial_figure_2D(self):
        initial = rectangle_p(1,1,1)
        fiber = cylinder(0.25)
        self.x = np.arange(0.0, 3.0, 0.01)
        if( view == 'x'):
            self.axes.plot(fiber.x, fiber.yup, 'k')
            self.axes.plot(fiber.x, fiber.ylo, 'k')
            self.axes.set_ylabel('z')
        if(view == 'y'):
            self.axes.plot(initial.l, initial.z, 'k')
            self.axes.plot(initial.l, initial.o, 'k')
            self.axes.plot(initial.o, initial.h, 'k')
            self.axes.plot(initial.z, initial.h, 'k')
            self.axes.plot(initial.l, np.ones(len(initial.o))*0.25, 'k')
            self.axes.plot(initial.l, np.ones(len(initial.o))*0.75, 'k')
            self.axes.set_ylabel('z')
        if(view == 'z'):
            self.axes.plot(initial.l, initial.z, 'k')
            self.axes.plot(initial.l, initial.o, 'k')
            self.axes.plot(initial.o, initial.h, 'k')
            self.axes.plot(initial.z, initial.h, 'k')
            self.axes.plot(np.ones(len(initial.o))*0.25, initial.l, 'k')
            self.axes.plot(np.ones(len(initial.o))*0.75, initial.l, 'k')
            self.axes.set_ylabel('y')
        self.axes.set_xlabel('x')
        
    def compute_initial_figure_3D(self):
        initial = rectangle_p(1,1,1)
        fiber = cylinder(0.25)
        #self.x = np.arange(0.0, 3.0, 0.01)
        #self.y = np.sin(2*np.pi*self.x)
        #self.z = np.cos(2*np.pi*self.x)
        #self.axes.plot(self.x, self.y, self.z)
        #self.axes.set_xlabel('x')
        #self.axes.set_ylabel('y(x)')
        self.axes.plot(initial.z, initial.z, initial.h, 'k')
        self.axes.plot(initial.z, initial.l, initial.z, 'k')
        self.axes.plot(initial.z, initial.l, initial.o, 'k')
        self.axes.plot(initial.w, initial.z, initial.z, 'k')
        self.axes.plot(initial.o, initial.l, initial.z, 'k')
        self.axes.plot(initial.w, initial.o, initial.z, 'k')
        self.axes.plot(initial.o, initial.o, initial.h, 'k')
        self.axes.plot(initial.o, initial.z, initial.h, 'k')
        self.axes.plot(initial.o, initial.l, initial.o, 'k')
        self.axes.plot(initial.w, initial.o, initial.o, 'k')
        self.axes.plot(initial.w, initial.z, initial.o, 'k')
        self.axes.plot(initial.z, initial.z,initial.h, 'k')
        self.axes.plot(initial.z, initial.o,initial.h, 'k')
        #for i in range(0, len(fiber.yup)):
        #    for j in range(0, len(fiber.y)):
        #        #self.axes.plot(fiber.x[i], fiber.y[j], fiber.yup[i])
        #        #self.axes.plot(fiber.x[i], fiber.y[j], fiber.ylo[i])
        for i in range(0, len(fiber.y)):
            self.axes.plot(np.ones(len(fiber.y))*fiber.x[i], fiber.y, np.ones(len(fiber.y))*fiber.yup[i], 'k')
            self.axes.plot(np.ones(len(fiber.y))*fiber.x[i], fiber.y, np.ones(len(fiber.y))*fiber.ylo[i], 'k')
        
    def redraw(self, x, y):
        self.axes.plot(x,y)
        self.draw() 

class rectangle_p():
    
    def __init__(self, length, height, width):
        self.z = [0,0,0,0,0]
        self.o = [1,1,1,1,1]
        self.l = [0,length/4,2*length/4,3*length/4,length]
        self.h = [0,height/4,2*height/4,3*height/4,height]
        self.w = [0,width/4,2*width/4,3*width/4,width]

class cylinder():
    
    def __init__(self, radius):
        self.x = np.linspace(0.25, 0.75)
        self.y = np.linspace(0,1)
        self.yup = []
        self.ylo = []
        for i in range(0, len(self.x)):
            temp = np.sqrt(radius**2 - (self.x[i] - 0.5)**2) + 0.5
            self.yup.append(temp)
            self.ylo.append(1 - temp)
            
        
        
            
        
    


app = QApplication(sys.argv)
widget = MainWindow()
widget.show()
app.exec_()