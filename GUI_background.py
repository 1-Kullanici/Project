from numeric_sim import N_simulation as nSim
from numeric_sim import C_profiles as cProf
from numeric_sim import Impurity

import matplotlib
matplotlib.use('QtAgg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

import numpy as np

from PyQt6.QtGui import QValidator
from PyQt6.QtCore import QSize, Qt
from PyQt6.QtWidgets import (
    QComboBox,
    QSpinBox,
    QDoubleSpinBox,
    QLabel,
    QMainWindow,
    QProgressBar,
    QPushButton,
    QVBoxLayout,
    QHBoxLayout,
    QGridLayout,
    QFormLayout,
    QWidget,
    QFileDialog,
    QMenuBar,
    QLineEdit,
    QApplication,
    QSizePolicy,
    QDialog,
    QMessageBox
)


class ScientificDoubleSpinBox(QDoubleSpinBox):
    def __init__(self, *args, **kwargs):
        super(ScientificDoubleSpinBox, self).__init__(*args, **kwargs)

    def textFromValue(self, value):
        return "{:.3e}".format(value)

    def validate(self, text, pos):
        try:
            if text == "" or text[-1] in 'e.-+':
                return (QValidator.State.Intermediate, text, pos)
            float(text)
            return (QValidator.State.Acceptable, text, pos)
        except ValueError:
            return (QValidator.State.Invalid, text, pos)
        

class MainWindow(QMainWindow):
    def __init__(self):                                 # Initialize the main window
        super().__init__()

        # Define default parameters
        global _Cb_default, _Cth_default, _Dopant_default, _T0_default, _T1_default, _xL_default, _xL_unit_default, _t0_default, _t1_default, _prgrss_default, _prgrss_lgnd_default, _prgrss_val_default, _xJun1_default, _xJun2_default
        _Cb_default          = 0     # in atoms/cm^3
        _Cth_default         = 1e15  # in atoms/cm^3
        _Dopant_default      = 2     # Boron
        _T0_default          = 900   # in °C
        _T1_default          = 900   # in °C
        _xL_default          = 600   # in µm or nm -> should be converted to cm before simulation
        _xL_unit_default     = 1     # 0: µm, 1: nm
        _t0_default          = 3000  # in s
        _t1_default          = 3000  # in s
        _prgrss_default      = 0
        _prgrss_lgnd_default = "Available"
        _prgrss_val_default  = "..."
        _xJun1_default       = "..."
        _xJun2_default       = "..."

        # Create a container widget
        widget = QWidget()
        self.setCentralWidget(widget)

        #################################
        ## Label/Widget for extra part ##
        #################################

        # Create a butoons for "About", "Help", and "Exit" parts
        aboutB = QPushButton("About")
        helpB = QPushButton("Help")
        exitB = QPushButton("Exit")

        # Change the button colors
        aboutB.setStyleSheet("color: #C0C0C0")
        helpB.setStyleSheet("color:  #C0C0C0")
        exitB.setStyleSheet("color:  #FF6347")

        # "About" button functionality
        aboutB.clicked.connect(lambda: self.show_popup(dialog_text  = "Numerical Simulation Tool for Diffusion\n\n Author: Ibrahim Furkan Tezcan\n\n Version: 1.1-stable\n\n Release: 09.Jun.2024\n\n Last Update: 10.Jun.2024",
                                                       dialog_title = "About"))
        

        # "Help" button functionality
        helpB.clicked.connect(lambda: self.show_popup(dialog_text  = "The simulation parameters are: the type of dopant, the temperature for predeposition, the temperature for drive-in, the spatial length, the predeposition time, and the drive-in time.\n\n The simulation results are the concentration profiles of the dopant in the silicon wafer and the junction depths for predeposition and drive-in. The concentration profiles are plotted in linear and logarithmic scales. The junction depths are given in the unit of the spatial length.\n\n The simulation is terminated by clicking the \"Simulate!\" button. The simulation can be terminated by clicking the \"Terminate!\" button during simulation. The parameters and plots can be reset by clicking the \"Reset!\" button.", 
                                                      dialog_title = "Help")
                                                      )
        
        # Exit button functionality
        exitB.clicked.connect(self.close)

        ######################################
        ## Label/Widget for plotting graphs ##
        ######################################

        # Create two figures
        self.figure_lin = Figure()
        self.figure_log = Figure()

        # Linear Scale Plot Labels/Canvas
        linearPlotCanvasTitle = QLabel("Concentration Profile (Linear Scale)")
        linearPlotCanvasTitle.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.linearPlotCanvas = FigureCanvasQTAgg(self.figure_lin)

        # Log Scale Plot Labels/Canvas
        logPlotCanvasTitle = QLabel("Concentration Profile (Logarithmic Scale)")
        logPlotCanvasTitle.setAlignment(Qt.AlignmentFlag.AlignCenter)                                  
        self.logPlotCanvas = FigureCanvasQTAgg(self.figure_log)

        # Create a toolbar for the canvas
        self.linPC_tb = NavigationToolbar(self.linearPlotCanvas, self)
        self.logPC_tb = NavigationToolbar(self.logPlotCanvas, self)

        ################################################
        ## Label/Widget for numerical simulation tool ##
        ################################################

        # Dopant selection panel
        Dopant_label = QLabel("Select impurity: ")
        self.Dopant_in = QComboBox()
        self.Dopant_in.addItem("Antimony (Sb)")
        self.Dopant_in.addItem("Arsenic (As)")
        self.Dopant_in.addItem("Boron (B)")
        self.Dopant_in.addItem("Phosphorus (P)")
        self.Dopant_in.setCurrentIndex(_Dopant_default)

        self.dopantProfile_list = self.createDopantProfile()
        self.Dopant = self.dopantProfile_list[self.Dopant_in.currentIndex()]
        
        # User input panel
        Cb_in = QLabel("Clipping conc. (Cb - Optional): ")    # Previously: Cb_in = QLabel("Choose minimum clipping concentration (Cb)")
        self.Cb  = ScientificDoubleSpinBox()
        self.Cb.setSuffix(" atoms/cm^3")
        self.Cb.setMinimum(0)
        self.Cb.setMaximum(1e18)
        self.Cb.setValue(_Cb_default)
        # print(Cb)

        Cth_in = QLabel("Background conc. (Cth): ")             # Previously: Cth_in = QLabel("Choose threshold concentration for junction depth calculation (Cth)")
        self.Cth  = ScientificDoubleSpinBox()
        self.Cth.setSuffix(" atoms/cm^3")
        self.Cth.setMinimum(0)
        self.Cth.setMaximum(1e18)
        self.Cth.setValue(_Cth_default)
        # print(Cth)

        T0_in = QLabel("Temperature for predep. (T0): ")
        self.T0  = QSpinBox()
        self.T0.setSuffix(" °C")
        self.T0.setMinimum(900)
        self.T0.setMaximum(1200)
        self.T0.setValue(_T0_default)
        # print(T0)

        T1_in = QLabel("Temperature for drive-in (T1): ")
        self.T1  = QSpinBox()
        self.T1.setSuffix(" °C")
        self.T1.setMinimum(900)
        self.T1.setMaximum(1200)
        self.T1.setValue(_T1_default)
        # print(T1)

        xL_in = QLabel("Spatial length (xd): ")
        self.xL_unit  = QComboBox()
        self.xL_unit.addItem("µm")
        self.xL_unit.addItem("nm")
        self.xL_unit.setCurrentIndex(_xL_unit_default)
        self.xL_inUnit  = QSpinBox()
        self.xL_inUnit.setMinimum(0)
        self.xL_inUnit.setMaximum(1000)
        self.xL_inUnit.setValue(_xL_default)
        self.xL = self.xL_unitConverter(self.xL_inUnit.value())
        # print(xL)
                
        t0_in = QLabel("Predeposition time (t0): ")
        self.t0  = QSpinBox()
        self.t0.setSuffix(" s")
        self.t0.setMinimum(0)
        self.t0.setMaximum(86400)          # 8.64e4 s = 24 hours
        self.t0.setValue(_t0_default)
        # print(t0)

        t1_in = QLabel("Drive-in time (t1): ")
        self.t1  = QSpinBox()
        self.t1.setSuffix(" s")
        self.t1.setMinimum(0)
        self.t1.setMaximum(86400)          # 8.64e4 s = 24 hours
        self.t1.setValue(_t1_default)
        # print(t1)

        # Clear button
        Clear = QPushButton("Reset!")

        # Clear button functionality
        Clear.clicked.connect(lambda: self.clear()) # Clear the canvas

        # "Simulate!" button
        self.startSimulation = QPushButton("Simulate!")
        self.startSimulation.setStyleSheet("color: #FFD700")

        # "Simulate!" button functionality 
        self.startSimulation.clicked.connect(lambda: self.toggle_simulation())

        # Progress bar
        self.progress = QProgressBar()
        self.progress.setMinimum(0)
        self.progress.setMaximum(100)
        self.progress.setValue(0)

        # Printed parts
        xJunc_1_label  = QLabel("Junction depth for predep.")
        xJunc_2_label  = QLabel("Junction depth for drive-in")

        self.xJunc_1   = QLabel("...")
        self.xJunc_2   = QLabel("...")

        self.xJunc_unit1    = QLabel(self.xL_unit.currentText())
        self.xJunc_unit2    = QLabel(self.xL_unit.currentText())
        
        self.progress_label = QLabel("Available")                          
        self.progress_val   = QLabel("...")

        """"""
        #################################
        ## Place widgets on the window ##
        #################################

        # layout = QVBoxLayout()
        layout = QGridLayout()

        layout.addWidget(Dopant_label,          0, 0, 1, 1)
        layout.addWidget(self.Dopant_in,        0, 1, 1, 2)

        layout.addWidget(Cb_in,                 1, 0, 1, 1)
        layout.addWidget(self.Cb,               1, 1, 1, 2)

        layout.addWidget(Cth_in,                2, 0, 1, 1)
        layout.addWidget(self.Cth,              2, 1, 1, 2)

        layout.addWidget(xL_in,                 3, 0, 1, 1)
        layout.addWidget(self.xL_inUnit,        3, 1, 1, 1)
        layout.addWidget(self.xL_unit,          3, 2, 1, 1)

        layout.addWidget(T0_in,                 0, 3, 1, 1)
        layout.addWidget(self.T0,               0, 4, 1, 1)
        layout.addWidget(T1_in,                 0, 5, 1, 1)
        layout.addWidget(self.T1,               0, 6, 1, 1)

        layout.addWidget(t0_in,                 1, 3, 1, 1)
        layout.addWidget(self.t0,               1, 4, 1, 1)
        layout.addWidget(t1_in,                 1, 5, 1, 1)
        layout.addWidget(self.t1,               1, 6, 1, 1)

        layout.addWidget(self.startSimulation,  2, 3, 1, 2)
        layout.addWidget(Clear,                 2, 5, 1, 2)

        layout.addWidget(self.progress,         3, 3, 1, 2)
        layout.addWidget(self.progress_val,     3, 5, 1, 1)
        layout.addWidget(self.progress_label,   3, 6, 1, 1)

        layout.addWidget(xJunc_1_label,         0, 7, 1, 2)
        layout.addWidget(self.xJunc_1,          1, 7, 1, 1)
        layout.addWidget(self.xJunc_unit1,      1, 8, 1, 1)

        layout.addWidget(xJunc_2_label,         2, 7, 1, 2)
        layout.addWidget(self.xJunc_2,          3, 7, 1, 1)
        layout.addWidget(self.xJunc_unit2,      3, 8, 1, 1)

        layout.addWidget(aboutB,                0, 9, 1, 2)
        layout.addWidget(helpB,                 1, 9, 1, 2)
        layout.addWidget(exitB,                 2, 9, 2, 2)

        layout.addWidget(linearPlotCanvasTitle, 4, 0, 1, 5)
        layout.addWidget(self.linearPlotCanvas, 5, 0, 1, 5)
        layout.addWidget(self.linPC_tb,         6, 0, 1, 5)

        layout.addWidget(logPlotCanvasTitle,    4, 5, 1, 6)
        layout.addWidget(self.logPlotCanvas,    5, 5, 1, 6)
        layout.addWidget(self.logPC_tb,         6, 5, 1, 6)

        # Create the menu bar
        # self._createMenuBar(layout)

        # Set the layout on the application's window
        widget.setLayout(layout)

        # Set the central widget of the Window. Widget will expand to take up all the space in the window by default.
        self.setCentralWidget(widget)

    def xL_unitConverter(self, xL_inUnit):              # Convert the unit of xL to cm  
        if self.xL_unit.currentIndex():
            return xL_inUnit*1e-7 # convert nm to cm
        else:
            return xL_inUnit*1e-4 # convert µm to cm

    def xL_unitConverter_inv(self, xL):                 # Convert the unit of xL to µm or nm
        if self.xL_unit.currentIndex():
            return xL*1e7         # convert cm to nm
        else:
            return xL*1e4         # convert cm to µm

    def plot(self, Cp_1:cProf, Cp_2:cProf, Cth_profile:cProf, x_ax:list, xJunc_1:float, xJunc_2:float): # Plot the results
        # Clear the previous plot (if any)
        self.figure_lin.clear()
        self.figure_log.clear()

        # Get Cth from the Cth_profile
        Cth = Cth_profile.arr[0]

        # Convert the junction depths to specified unit
        xJunc_1_inUnit = self.xL_unitConverter_inv(xJunc_1)
        xJunc_2_inUnit = self.xL_unitConverter_inv(xJunc_2)

        # Set the progress bar: 30% complete
        self.updateProgress(40)

        # Linear plot
        ax1 = self.figure_lin.add_subplot()
        ax1.set_xlabel('Position ({})'.format(self.xL_unit.currentText()))
        ax1.set_ylabel('Concentration (atoms/cm^3)')
        ax1.plot(x_ax, Cp_1.arr, color='C0', label='Predep.',  zorder=1)
        ax1.plot(x_ax, Cp_2.arr, color='C1', label='Drive-in', zorder=0)
        ax1.plot(x_ax, Cth_profile.arr, color='C2', label='Background Conc.', linestyle='dashed', zorder=2)
        ax1.legend(loc='upper right', prop={'size': 6})                             # Set the legend location
        
        # Set the progress bar: 60% complete
        self.updateProgress(70)

        # Logarithmic plot
        ax2 = self.figure_log.add_subplot()
        ax2.set_yscale('log')                                                       # Set the y-axis to log scale
        ax2.set_xlabel('Position ({})'.format(self.xL_unit.currentText()))
        ax2.set_ylabel('Concentration (atoms/cm^3)')
        ax2.plot(x_ax, Cp_1.arr, color='C0', label='Predep.',  zorder=1)
        ax2.plot(x_ax, Cp_2.arr, color='C1', label='Drive-in', zorder=0)
        ax2.plot(x_ax, Cth_profile.arr, color='C2', label='Background Conc.', linestyle='dashed', zorder=2)
        ax2.scatter(xJunc_1_inUnit, Cth, color='brown', marker='o', label='Junction Depth for Predep.',  zorder=3)
        ax2.scatter(xJunc_2_inUnit, Cth, color='black', marker='o', label='Junction Depth for Drive-in', zorder=3)
        ax2.legend(loc='upper right', prop={'size': 6})                             # Set the legend location
        ax2.set_ylim(bottom=1)                                                      # Set the y-axis limits

        # Set the progress bar: 90% complete
        self.updateProgress(90)

        # Enable tight layout
        self.figure_lin.tight_layout()
        self.figure_log.tight_layout()

        # Set the progress bar: 95% complete
        self.updateProgress(95)

        # Draw the plot
        self.linearPlotCanvas.draw()
        self.logPlotCanvas.draw()

    def createDopantProfile(self):                      # Create a list of dopant profiles
        imp_Sb = Impurity(4.58, 3.88, 1e20)
        imp_As = Impurity(9.17, 3.99, 2e21)
        imp_B  = Impurity(1.0,  3.5,  3e20)
        imp_P  = Impurity(4.7,  3.68, 1e21)

        self.dopantProfile_list = [imp_Sb, imp_As, imp_B, imp_P]
        return self.dopantProfile_list

    def updateDopantProfile(self, dopant_idx:int):      # Update the dopant profile to selected one
        self.Dopant=self.dopantProfile_list[dopant_idx]
        # self.Do, self.Ea, self.Co = self.Dopant.get_attr()

    def updateJuncDepth(self, label:QLabel, val:float, unit:QLabel): # Update the junction depths
        val_inUnit = self.xL_unitConverter_inv(val)
        label.setText('{:.2f}'.format(val_inUnit))
        unit.setText(self.xL_unit.currentText())

    def updateProgressLabel(self, text:str):            # Update the progress label and the GUI
        self.progress_label.setText(text)
        QApplication.processEvents()  # Process all pending events

    def updateProgress(self, value:int, unit:str='%'):  # Update the progress bar
        self.progress.setValue(value)
        self.progress_val.setText("{}{}".format(value, unit))
        QApplication.processEvents()  # Process all pending events

    def updateParameters(self):                         # Update input parameters
        self.updateDopantProfile(self.Dopant_in.currentIndex())
        self.xL = self.xL_unitConverter(self.xL_inUnit.value())
        kwarg = {

            "Cb"        : self.Cb.value(),
            "Cth"       : self.Cth.value(),
            "Dopant"    : self.Dopant,
            "T0"        : self.T0.value(),
            "T1"        : self.T1.value(),
            "xL"        : self.xL,
            "t0"        : self.t0.value(),
            "t1"        : self.t1.value(),

        }

        print("Updated parameters are: ", kwarg)
        return kwarg
    
    def kwargParser(self, **kwargs):                    # Parse the keyword arguments
        values = []
        for _, value in kwargs.items():
            values.append(value)
        return tuple(values)

    def resetParameters(self):                          # Reset parameters to default values
        self.Dopant_in.setCurrentIndex(_Dopant_default)
        self.Cb.setValue(_Cb_default)
        self.Cth.setValue(_Cth_default)
        self.T0.setValue(_T0_default)
        self.T1.setValue(_T1_default)
        self.xL_unit.setCurrentIndex(_xL_unit_default)
        self.xL_inUnit.setValue(_xL_default)
        self.t0.setValue(_t0_default)
        self.t1.setValue(_t1_default)
        self.resetProgress()
    
    def resetProgress(self):                            # Reset the progress bar
        self.progress.setValue(_prgrss_default)
        self.progress_label.setText(_prgrss_lgnd_default)
        self.progress_val.setText(_prgrss_val_default)

    def resetJuncDepth(self):                           # Reset the junction depths
        self.xJunc_1.setText(_xJun1_default)
        self.xJunc_2.setText(_xJun2_default)

    def clearCanvas(self):                              # Reset everything on the canvas
        # Clear the previous plot (if any)
        self.figure_lin.clear()
        self.figure_log.clear()
        # Draw the blank figure
        self.linearPlotCanvas.draw()
        self.logPlotCanvas.draw()

    def clear(self):                                    # Reset everything to default values
        # Progress bar notificaiton
        self.updateProgressLabel("Clearing...")
        self.updateProgress(0)

        # Clear the transformation parameters
        self.resetParameters()
        self.updateProgress(33)

        # Clear the junction depths
        self.resetJuncDepth()
        self.updateProgress(66)

        # Clear canvas
        self.clearCanvas()
        self.updateProgress(100)

        # Clear the progress bar
        self.updateProgressLabel("Done!")
        self.resetProgress()
        print("---------CLEARED---------")

    def toggle_simulation(self):                        # Toggle the Simulate! button
        if self.startSimulation.text() == "Simulate!":
            self.startSimulation.setText("Terminate!")
            self.startSimulation.setStyleSheet("color: #DC143C")
            self.simulate()
            self.startSimulation.setText("Simulate!")
            self.startSimulation.setStyleSheet("color: #FFD700")
        else:
            self.startSimulation.setText("Terminating...")
            self.startSimulation.setStyleSheet("color: #800000")
            self.startSimulation.setEnabled(False)
            self.terminate_simulation()
            self.startSimulation.setText("Simulate!")
            self.startSimulation.setStyleSheet("color: #FFD700")
            self.startSimulation.setEnabled(True)

    def simulate(self):                                 # Run the simulation
        self.clearCanvas()
        self.resetJuncDepth()
        self.resetProgress()
        self.updateProgressLabel("Preparing...")

        # Set initial parameters
        param =  self.updateParameters()
        Cb, Cth, Dopant, T0, T1, xD, t0, t1 = self.kwargParser(**param)

        # Create instances of the N_simulation class
        self.preDep  = nSim(Dopant, T0)
        self.driveIn = nSim(Dopant, T1)

        # Set the progress bar: 33% complete
        self.updateProgress(33)

        # Set the simulation parameters
        x_i  = int(xD/self.preDep.x_step)+1    #x_step is set to either 1e-8 cm (1 Angstrom) or 1e-7 (1 nm) inside nSim - constant for both simulations
        t_j0 = int(t0/self.preDep.t_step)+1    #Number of time iterations for predep.
        t_j1 = int(t1/self.driveIn.t_step)+1   #Number of time iterations for drive-in
        # print("# of iterations for predep.: ", t_j0, "\n# of iterations for drive-in: ", t_j1)
        
        # Set the progress bar: 66% complete
        self.updateProgress(66)

        # Create an instance of the C_profiles class
        C_1 = cProf(x_i=x_i, Cb=Cb)
        C_2 = cProf(x_i=x_i, Cb=Cb)
        Cth_profile = cProf().create_empty_profile(x_i=x_i, Cb=Cth)
        
        # Set the progress bar: 100% complete
        self.updateProgress(100)
        self.updateProgressLabel("Done!")

        # Set the progress bar: 0% complete
        self.updateProgress(0)
        self.updateProgressLabel("Running Simulation: 1/2...")

        # Run the simulation for predep.
        Cp_1, xjunc_1 = self.preDep.lumerical_on_budget(C_1, Cb, Cth, t_j=t_j0, 
                                                        progressPercentageOutput=self.updateProgress, 
                                                        progressOutput=self.updateProgressLabel,
                                                        )
        self.updateProgressLabel("Predep. sim. is completed. Saving profile...")
        
        # Print the junction depths:
        self.updateJuncDepth(self.xJunc_1, xjunc_1, self.xJunc_unit1)

        # Set the progress bar: 0% complete
        self.updateProgress(0)
        self.updateProgressLabel("Running Simulation: 2/2...")

        # Set the initial profile for the next simulation
        C_2.Cold = Cp_1                   

        # Run the simulation for drive-in
        Cp_2, xjunc_2 = self.driveIn.lumerical_on_budget(C_2, Cb, Cth, t_j=t_j1,
                                                         process=1, 
                                                         progressPercentageOutput=self.updateProgress, 
                                                         progressOutput=self.updateProgressLabel,
                                                         )
        self.updateProgressLabel("Drive-in simulation completed. Saving profile...")

        # Print the junction depths:
        self.updateJuncDepth(self.xJunc_2, xjunc_2, self.xJunc_unit2)

        # Set the progress bar: 0% complete
        self.updateProgress(0)
        self.updateProgressLabel("Plotting the profiles...")

        # Cut the first element of the profile to avoid the initial condition
        Cp_1.cut_initial()
        Cp_2.cut_initial()
        Cth_profile.cut_initial()

        # Construct the x-axis array
        x_step_inUnit = self.xL_unitConverter_inv(self.preDep.x_step) # Convert the x_step to the specified unit
        range_arr = np.array(range(Cth_profile.size()))
        x_ax = range_arr*x_step_inUnit                                # Construct the x-axis array in specified unit

        # Set the progress bar: 10% complete
        self.updateProgress(20)

        # Plot the results
        self.plot(Cp_1, Cp_2, Cth_profile, x_ax, xjunc_1, xjunc_2)

        # Set the progress bar: 100% complete
        self.updateProgress(100)
        self.updateProgressLabel("Done!")

        # # Clear the junk
        # del self.preDep, self.driveIn
        # del C_1, C_2, Cp_1, Cp_2, Cth_profile

    def terminate_simulation(self):                     # Terminate the simulation - Not working properly
        self.preDep.terminate()
        self.driveIn.terminate()
    
    def show_popup(self, dialog_text:str="You've found me!", dialog_title:str="Easter Egg"): # Show a popup window
        msg = QMessageBox()
        msg.setWindowTitle(dialog_title)
        msg.setText(dialog_text)
        msg.setIcon(QMessageBox.Icon.Information)
        msg.setStandardButtons(QMessageBox.StandardButton.Ok)
        msg.exec()


# TODO: (Optional) Try to handle the boundary conditions in a more elegant way (indx 0 and negative indexes) for both processes
 