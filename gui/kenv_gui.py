
# A. Petrenko (Snezhinsk, 2019)
# petrenko@inp.nsk.su

import sys, time
import numpy as np

# Kapchinsky-Vladimirsky equation solver developped by V. Fedorov, D. Nikiforov, A. Petrenko
# https://github.com/fuodorov/kenv
import kenv as kv

from PyQt5.QtWidgets import QApplication, QPushButton, QVBoxLayout, \
    QHBoxLayout, QLabel, QSlider, QDoubleSpinBox, QWidget, QLineEdit, \
    QComboBox, QFileDialog
from PyQt5.QtCore import Qt

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
import json

print("KENV version = %s" % kv.__version__ )

mc2 = 0.511 # MeV

class Window(QWidget):
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)
        self.updating = False # needed in order to update everything for loaded config
        
        self.ConfigFile = "config.json"
        
        self.load_config(self.ConfigFile)

        self.SolNames = list(acc.Bz_beamline.keys())
        self.AccNames = list(acc.Ez_beamline.keys())
        
        self.Sol2Vary = self.SolNames[1]
        self.Acc2Vary = self.AccNames[0]

        # a figure instance to plot on
        self.figure = plt.figure()

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)
        
        self.populate_figure()

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)

        lay1 = QHBoxLayout()

        plt_lay = QVBoxLayout()
        plt_lay.addWidget(self.toolbar)
        plt_lay.addWidget(self.canvas)
        #plt_lay.addWidget(self.button)

        lay1.addLayout(plt_lay)

        ctr_widget = QWidget()

        ctrl_lay = QVBoxLayout()
        ctrl_lay.addStretch()
        
        ctr_widget.setLayout(ctrl_lay)
        ctr_widget.setFixedWidth(150)

        ctrl_lay.addWidget(QLabel("Beam current (A):"))
        self.IEdit = QLineEdit()
        ctrl_lay.addWidget(self.IEdit)
        self.IEdit.setText("%.0f" % beam.current) # A
        self.IEdit.editingFinished.connect(self.update_I)

        ctrl_lay.addWidget(QLabel("Solenoid field (Gs):"))
        
        self.solSelector = QComboBox()
        ctrl_lay.addWidget(self.solSelector)
        self.solSelector.addItems(self.SolNames)
        
        self.solSelector.setCurrentText(self.Sol2Vary)
        self.solSelector.currentIndexChanged.connect(self.sol_selected)       

        self.solEdit = QLineEdit()
        self.solEdit.setText("%.1f" % self.solB()) # Gs
        ctrl_lay.addWidget(self.solEdit)
        self.solEdit.editingFinished.connect(self.update_solenoid)

        ctrl_lay.addWidget(QLabel("Electric field (kV/m):"))
        
        self.AccSelector = QComboBox()
        ctrl_lay.addWidget(self.AccSelector)
        self.AccSelector.addItems(self.AccNames)
        
        self.AccSelector.setCurrentText(self.Acc2Vary)
        self.AccSelector.currentIndexChanged.connect(self.acc_selected)

        self.AccEdit = QLineEdit()
        self.AccEdit.setText("%.0f" % self.AccE()) # kV/m
        ctrl_lay.addWidget(self.AccEdit)
        self.AccEdit.editingFinished.connect(self.update_Acc_module)


        self.savebtn=QPushButton('Save config')
        ctrl_lay.addWidget(self.savebtn)
        #self.button2.resize(self.button2.sizeHint())
        self.savebtn.clicked.connect(self.save_config)

        self.saveasbtn=QPushButton('Save config as')
        ctrl_lay.addWidget(self.saveasbtn)
        #self.button2.resize(self.button2.sizeHint())
        self.saveasbtn.clicked.connect(self.save_config_as)

        self.loadbtn=QPushButton('Load config')
        ctrl_lay.addWidget(self.loadbtn)
        self.loadbtn.clicked.connect(self.load_config_clicked)
    
        #lay1.addLayout(ctrl_lay)
        lay1.addWidget(ctr_widget)

        self.setLayout(lay1)

        self.setGeometry(150, 150, 1200, 700)
        
        self.show()
        
    def populate_figure(self):
        
        self.figure.clear()
        self.axBz = plt.subplot2grid((7,1), (5,0), rowspan=2)
        self.axx = plt.subplot2grid((7,1), (0,0), sharex=self.axBz, rowspan=5)

        acc.compile()
        sim = kv.Simulation(beam, acc)
        sim.track()

        z  = acc.z
        
        SolNames = list(acc.Bz_beamline.keys())
        AccNames = list(acc.Ez_beamline.keys())
        
        for name in SolNames:
            z0 = acc.Bz_beamline[name].z0
            self.axBz.text(z0,-70.0, name, verticalalignment='top', horizontalalignment='center',
                      color='blue', fontsize=9, rotation='vertical', alpha=0.7, picker=5)
            
        for name in AccNames:
            z0 = acc.Ez_beamline[name].z0
            self.axBz.text(z0,+70.0, name, verticalalignment='bottom', horizontalalignment='center',
                      color='red', fontsize=9, rotation='vertical', alpha=0.7, picker=5)

        x_size = sim.envelope_x(z)*100.0 # cm
        self.x_size_line, = self.axx.plot(z, x_size, '-', linewidth=2, color='blue', alpha=0.6)

        Ekin = sim.gamma(z)*mc2 - mc2 # MeV
        self.Ekin_line, = self.axx.plot(z, Ekin, '-', linewidth=2, color='red', alpha=0.6)

        z0 = self.sol_z0()
        self.VLine_axx, = self.axx.plot((z0,z0), (-20,+20),     '-', linewidth=8, alpha=0.20, color='gray')
        self.VLine_axBz, = self.axBz.plot((z0,z0), (-5000,+5000), '-', linewidth=8, alpha=0.20, color='gray')

        self.axBz.plot((0,acc.z_stop), (0,0), linewidth=1, alpha=0.6, color='black') # zero line

#        z0 = self.Acc_z0()
#        self.AccVLine_axx, = self.axx.plot((z0,z0), (-20,+20),     '-', linewidth=8, alpha=0.15, color='red')
#        self.AccVLine_axBz, = self.axBz.plot((z0,z0), (-5000,+5000), '-', linewidth=8, alpha=0.15, color='red')

        self.Bz_line, = self.axBz.plot(z, 1e4*acc.Bz(z), '-', color='blue', alpha=0.5)
        self.Ez_line, = self.axBz.plot(z, 1e3*acc.Ez(z), '-', color='red', alpha=0.5)

        self.axBz.set_xlabel("z (m)")

        self.axBz.set_ylabel("$E_z$ (kV/m), $B_z$ (Gs)")
        self.axx.set_ylabel("$x$ (cm), $E_{\mathrm{kin}} (MeV)$")
        self.axx.grid(True);
        self.axx.set_ylim((0,10)) # cm
        self.axBz.set_ylim((-1200,1200)) # Gs

        self.figure.tight_layout()
        self.canvas.draw()
        self.figure.canvas.mpl_connect('pick_event', self.onpick)
    
    def onpick(self, event):
        name = event.artist.get_text()
        if name in self.SolNames:
            if self.Sol2Vary == name:
                self.sol_selected()
            else:
                self.solSelector.setCurrentText(name)
        elif name in self.AccNames:
            if self.Acc2Vary == name:
                self.acc_selected()
            else:
                self.AccSelector.setCurrentText(name)

    def solB(self):
        return acc.Bz_beamline[self.Sol2Vary].max_field*1e4 # Gs

    def AccE(self):
        return acc.Ez_beamline[self.Acc2Vary].max_field*1e3 # kV/m

    def sol_z0(self):
        return acc.Bz_beamline[self.Sol2Vary].z0 # m

    def Acc_z0(self):
        return acc.Ez_beamline[self.Acc2Vary].z0 # m

    def sol_selected(self):
        if self.updating: return
        self.Sol2Vary = self.solSelector.currentText()
        self.solEdit.setText("%.1f" % self.solB()) # Gs
        z0 = self.sol_z0()
        self.VLine_axx.set_xdata((z0,z0))
        self.VLine_axBz.set_xdata((z0,z0))
        self.canvas.draw()

    def acc_selected(self):
        if self.updating: return
        self.Acc2Vary = self.AccSelector.currentText()
        self.AccEdit.setText("%.0f" % self.AccE()) # kV/m
        z0 = self.Acc_z0()
        self.VLine_axx.set_xdata((z0,z0))
        self.VLine_axBz.set_xdata((z0,z0))
        self.canvas.draw()

    def update_solenoid(self):
        oldValue = self.solB()
        try:
            newValue = float(self.solEdit.text())
            if np.abs(newValue - oldValue) > 0.01:
                acc.Bz_beamline[self.Sol2Vary].max_field = newValue/1e4
                self.update_plot()
        except Exception as e:
            print(e)
            self.solEdit.setText("%.1f" % oldValue) # Gs

    def update_Acc_module(self):
        oldValue = self.AccE()
        try:
            newValue = float(self.AccEdit.text())
            if np.abs(newValue - oldValue) > 0.01:
                acc.Ez_beamline[self.Acc2Vary].max_field = newValue/1e3
                self.update_plot()
        except Exception as e:
            print(e)
            self.AccEdit.setText("%.0f" % oldValue) # kV/m


    def update_I(self):
        oldValue = beam.current
        try:
            newValue = float(self.IEdit.text())
            if np.abs(newValue - oldValue) > 0.01:
                beam.current = newValue # A
                self.update_plot()
        except Exception as e:
            print(e)
            self.IEdit.setText("%.0f" % oldValue) # I

    def update_plot(self):
        #print("update_plot()")

        z = acc.z

        acc.compile()
        sim = kv.Simulation(beam, acc)
        sim.track()

        x_size = sim.envelope_x(z)*100 #+ rnd_data*2 # cm
        self.x_size_line.set_ydata(x_size)

        self.Bz_line.set_ydata(1e4*acc.Bz(z))
        self.Ez_line.set_ydata(1e3*acc.Ez(z))

        Ekin = sim.gamma(z)*mc2 - mc2 # MeV
        self.Ekin_line.set_ydata(Ekin)

        # refresh canvas
        self.canvas.draw()
    
    def save_config(self):
        
        config = {
                "date": time.ctime(),
                "KENV version": kv.__version__,
                "Accelerator": {'z_start': acc.z_start, 'z_stop': acc.z_stop, 'dz': acc.dz},
                "Beam": {'current': beam.current, 'energy': beam.energy, 'radius': beam.radius,
                        'rp': beam.rp, 'normalized_emittance': beam.normalized_emittance,
                        'x': beam.x, 'y': beam.y, 'xp': beam.xp, 'yp': beam.yp},
        }
        
        arr = []
        for itm in acc.Bz_beamline.values():
            arr.append({'name': itm.name, 'z0': itm.z0, 'file_name': itm.file_name,
                        'max_field': itm.max_field})
        config['Bz_beamline'] = arr
        
        arr = []
        for itm in acc.Ez_beamline.values():
            arr.append({'name': itm.name, 'z0': itm.z0, 'file_name': itm.file_name,
                        'max_field': itm.max_field})
        config['Ez_beamline'] = arr
        
        with open(self.ConfigFile, "w") as f:
            json.dump(config, f, indent=1)

        self.setWindowTitle(self.ConfigFile)

    def save_config_as(self):
        path = QFileDialog.getSaveFileName(self, "Save As")[0]
        if path:
            self.ConfigFile = path
            self.save_config()

    def load_config_clicked(self):
        path = QFileDialog.getOpenFileName(self, "Open")[0]
        if path:
            self.load_config(path)
            self.apply_config()
    
    def load_config(self, fname):
        global acc, beam
        try:
            print('Loading %s' % fname)
            with open(fname, 'r') as f:
                cfg = json.load(f)
            
            itm = cfg['Accelerator']
            acc = kv.Accelerator(z_start=itm['z_start'],
                                 z_stop=itm['z_stop'],
                                 dz=itm['dz'])
            itm = cfg['Beam']

            beam = kv.Beam(
                    current=itm['current'], energy=itm['energy'], radius=itm['radius'],
                    rp=itm['rp'], normalized_emittance=itm['normalized_emittance'],
                    x=itm['x'], y=itm['y'], xp=itm['xp'], yp=itm['yp'],
            )
            
            print(beam)
            
            for itm in cfg['Bz_beamline']:
                acc.add_solenoid(name=itm['name'], center=itm['z0'],
                                 file_name=itm['file_name'],
                                 max_field=itm['max_field'])

            for itm in cfg['Ez_beamline']:
                acc.add_accel(name=itm['name'], center=itm['z0'],
                                 file_name=itm['file_name'],
                                 max_field=itm['max_field'])
            
            self.ConfigFile = fname
            self.setWindowTitle(self.ConfigFile)
            
        except Exception as e:
            print(e)

    def apply_config(self):
        self.updating = True
        
        print("Applying config.")
        self.populate_figure()
        self.IEdit.setText("%.0f" % beam.current) # A

        self.solSelector.clear()
        self.solSelector.addItems(self.SolNames)

        self.AccSelector.clear()
        self.AccSelector.addItems(self.AccNames)

        if not self.Sol2Vary in self.SolNames:
            self.Sol2Vary = self.SolNames[0]
        
        self.solSelector.setCurrentText(self.Sol2Vary)
        self.solEdit.setText("%.1f" % self.solB()) # Gs

        if not self.Acc2Vary in self.AccNames:
            self.Acc2Vary = self.AccNames[0]
        
        self.AccSelector.setCurrentText(self.Acc2Vary)
        self.AccEdit.setText("%.1f" % self.AccE()) # kV/m

        self.updating = False
        
app = QApplication(sys.argv)
main = Window()

app.exec_()
