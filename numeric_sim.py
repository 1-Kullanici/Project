###########################################
## DEVELOPER: Ibrahim Furkan Tezcan      ##
## Version: 1.0                          ##  
## Date: 06.Jun.2024                     ##
###########################################

import sys
import numpy as np
import matplotlib.pyplot as plt
# import argparse
# from numpy.lib.stride_tricks import as_strided as ast

class Impurity:     # Dopant properties - TODO: CHECK "Co". It may depend on temperature.
    """
        This class contains attributes and properties of an impurity object.
    """
    def __init__(self, Do:float, Ea:float, Co:float):
        self.Do = Do
        self.Ea = Ea
        self.Co = Co 
    
    def get_attr(self) -> tuple:
        return (self.Do, self.Ea, self.Co)


class _C_profile:   #Concentration profile
    """
        This class contains attributes and properties of a concentration profile object.\n
        Extends C_profiles class.
    """
    def __init__(self, x_i:int=1, Cb:float=0):
        self.arr = np.ones(x_i)*Cb   # Creates a blank concentration profile when initialized

    def size(self) -> int:
        return self.arr.size
    
    def get_val(self, i:int) -> float:
        return self.arr[i]

    def set_val(self, val, i:int):
        self.arr[i] = val   

    def get_profile(self) -> np.array:
        return self.arr
    
    def cut_initial(self):
        self.arr = self.arr[1:]


class C_profiles:   # Conc. profiles at two different times
    """
        This object contains two concentration profiles at two different times.\n
        Concentration profile is a 1D array that contains the concentration values of a dopant at sampled positions.
    """

    def __init__(self, x_i:int=1, Cb:float=0):
        self.Cnew = _C_profile(x_i, Cb)   # Creates a blank concentration profile for t_j when initialized
        self.Cold = _C_profile(x_i, Cb)   # Creates a blank concentration profile for t_j-1 when initialized
    
    def create_empty_profile(self, x_i:int=1, Cb:float=0) -> _C_profile:
        return _C_profile(x_i, Cb)
    
    def make_profiles(self, Cnew:_C_profile, Cold:_C_profile):
        self.Cnew = Cnew
        self.Cold = Cold

    def get_profiles(self) -> tuple:
        return (self.Cnew, self.Cold)
    
    def set_profile(self, C:_C_profile, i:int):
        if i == 0:
            self.Cnew = C
        elif i == 1:
            self.Cold = C
        else:
            print("Invalid index.")

    def size(self) -> int:
        return self.Cnew.size()

    def update_profiles(self, Cij:_C_profile):
        if self.size() == Cij.size():
            self.Cold = self.Cnew
            self.Cnew = Cij
        else:
            print("Size mismatch. Profiles not updated.")


class N_simulation: # Simulation class
    """
        This class performs numerical simulations using difference equation derived from diffusion equation.
    """
    def __init__(self, dopant:Impurity=None, T:int=900):                                                                #DONE!
        # Termination flag
        self.terminateFlag = False

        # Impurity parameters
        self.Ea = dopant.Ea                             #eV
        self.D0 = dopant.Do                             #cm^2/s
        self.C0 = dopant.Co                             #cm^-3
        
        # Constants
        self.Boltzmann = 8.617e-5                       #eV/K
        # self.x_step = 1e-8                              #cm (1 Angstrom)
        self.x_step = 1e-7                              #cm (1 nm)
        print("Position step: ", self.x_step)

        # Calculate diffusivity based on dopant
        self.T = T + 273.15                             #convert °C to K 
        self.D = self.diffusivity(self.T)    

        # Calculate time step for convergence
        self.t_step = (self.x_step**2) / (2*self.D)     #seconds
        print("Time step: ", self.t_step)    

    def diffusivity(self, T=900) -> float:                                                                              #DONE!
        """
            diffusivity(T=900)
            
        Calculates diffusivity given the parameters.
        
        Parameters:
        --------------------------------
        T   -   Temperature for process (degree C)  : int
        """
        return (self.D0 * np.exp(-self.Ea/(self.Boltzmann) * (1/T)))

    def lumerical_on_budget(self, C:C_profiles, Cb:float=0, Cth:float=1e15, t_j:int=1, process:bool=0, progressPercentageOutput=print, progressOutput=print) -> _C_profile:    #DONE!
        """
            lumerical_on_budget(C, Cb=0, Cth=100, t_j=1, process=0, progressPercentageOutput=print, progressOutput=print)
            
        Numerically calculates the concentration profile of given dopant.\n
        Cb must be smaller than Cth.\n
        If process is 0 (set by default), calculation will be done for predeposition.\n
        If process is 1, calculation will be done for drive-in.
        
        Parameters:
        --------------------------------
        C                        -   Concentration profile provided                  : _C_profile
        Cb                       -   Bottom concentration clip (atoms/cm^3)          : float
        Cth                      -   Threshold (backgrnd) concentration (atoms/cm^3) : float
        t_j                      -   Time iteration for process (seconds)            : int
        process                  -   Selected process (predep./drive-in)             : bool
        progressPercentageOutput -   Function to print the progress percentage       : function
        progressOutput           -   Function to print the progress                  : function
        """

        xjunc=0
        coef = self.D*self.t_step/(self.x_step**2)

        if process == 0:
            # j-1 iteration of time
            for j in range(1, t_j):
                C.Cold.set_val(self.C0, 0) #set initial condition and boundary condition
                # progressOutput(100*j/t_j, "%", "completed.", end="\r") 
                progressPercentageOutput(int(100*j/t_j))
                for i in range(1, C.Cold.size()-1):
                    Cij=C.Cold.get_val(i) + coef * (C.Cold.get_val(i+1) - 2*C.Cold.get_val(i) + C.Cold.get_val(i-1))
                    C.Cnew.set_val(Cij, i)
                C.update_profiles(C.Cnew)
                if self.terminateFlag:
                    progressOutput("Simulation is terminated.")
                    break
            
            # Find junction depth
            xjunc_idx=0
            min_diff = float('inf')
            for i in range(1, C.Cnew.size()-1):
                current_diff =  abs(C.Cnew.get_val(i) - Cth)
                if current_diff < min_diff:
                    min_diff = current_diff
                    xjunc_idx = i
            xjunc = xjunc_idx*self.x_step

        elif process == 1:
            xjunc_idx=0
            min_diff = float('inf')
            # j-1 iteration of time
            for j in range(1, t_j):
                C.Cold.set_val(Cb, -1) #set initial condition and boundary condition
                # progressOutput(100*j/t_j, "%", "completed.", end="\r")     
                progressPercentageOutput(int(100*j/t_j))
                for i in range(1, C.size()-1):
                    Cij=C.Cold.get_val(i) + coef * (C.Cold.get_val(i+1) - 2*C.Cold.get_val(i) + C.Cold.get_val(i-1))
                    C.Cnew.set_val(Cij, i)
                C.update_profiles(C.Cnew)
                if self.terminateFlag:
                    progressOutput("Simulation is terminated.")
                    break

            # Find junction depth
            xjunc_idx=0
            min_diff = float('inf')
            for i in range(1, C.Cnew.size()-1):
                current_diff =  abs(C.Cnew.get_val(i) - Cth)
                if current_diff < min_diff:
                    min_diff = current_diff
                    xjunc_idx = i
            xjunc = xjunc_idx*self.x_step

        else:
            progressOutput("Process not selected properly. Returning given profile.")

        Cn = C.get_profiles()[0]
        return Cn, xjunc

    def terminate(self):
        self.terminateFlag = True
         

class plot: #TESTING...
    """
        This class contains methods to plot the data.
    """
    def __init__(self, Cp_1:_C_profile=None, Cp_2:_C_profile=None, Cth_profile:_C_profile=None, xJunc_1:float=0, xJunc_2:float=0):
        fig, (ax1, ax2) = plt.subplots(2)
        fig.suptitle('Concentration Profile of the Dopant')
        # Get Cth from the Cth_profile
        Cth = Cth_profile.arr[0]

        # Linear plot
        ax1.plot(0,0)
        ax1.set_title('Linear Plot')
        ax1.set_xlabel('Position ({})'.format('cm'))
        ax1.set_ylabel('Concentration (atoms/cm^3)')
        ax1.plot(Cp_1.arr, color='blue', label='Predep.')
        ax1.plot(Cp_2.arr, color='red',  label='Drive-in')
        ax1.plot(Cth_profile.arr, color='green', label='Background Conc.', linestyle='dashed')
        ax1.legend(loc='upper right', prop={'size': 5})                             # Set the legend location

        # Logarithmic plot
        ax2.plot(0,1)
        ax2.set_yscale('log') 
        ax2.set_title('Logarithmic Plot')                                                      # Set the y-axis to log scale
        ax2.set_xlabel('Position ({})'.format('cm'))
        ax2.set_ylabel('Concentration (atoms/cm^3)')
        ax2.plot(Cp_1.arr, color='blue', label='Predep.')
        ax2.plot(Cp_2.arr, color='red',  label='Drive-in')
        ax2.plot(Cth_profile.arr, color='green', label='Background Conc.', linestyle='dashed')
        ax2.scatter(xJunc_1, Cth, color='cyan', marker='o', label='Junction Depth for Predep.')
        ax2.scatter(xJunc_2, Cth, color='magenta', marker='o', label='Junction Depth for Drive-in')
        ax2.legend(loc='upper right', prop={'size': 5})                             # Set the legend location
        ax2.set_ylim(bottom=1) 

    def plot_all(self, C:_C_profile=None):
        plt.tight_layout()
        plt.show()


def createDopantProfile(dopant_idx:int):        # DONE!
    imp_Sb = Impurity(4.58, 3.88, 1e20)
    imp_As = Impurity(9.17, 3.99, 2e21)
    imp_B  = Impurity(1.0,  3.5,  3e20)
    imp_P  = Impurity(4.7,  3.68, 1e21)

    dopantProfile_list = [imp_Sb, imp_As, imp_B, imp_P]
    return dopantProfile_list[dopant_idx]

def main(Cb, Cth, T0, T1, x, t0, t1):           #TESTING...

    # Create an instance of the Impurity class
    impurity = createDopantProfile(1)

    # Create instances of the N_simulation class
    preDep  = N_simulation(impurity, T0)
    driveIn = N_simulation(impurity, T1)

    # Set the simulation parameters
    x_i  = int(x/preDep.x_step)+1    #x_step is 1e-8 cm (1 Angstrom) - constant for both simulations
    t_j0 = int(t0/preDep.t_step)+1   #Number of time iterations for predep.
    t_j1 = int(t1/driveIn.t_step)+1  #Number of time iterations for drive-in
    print("# of iterations for predep.: ", t_j0, "\n# of iterations for drive-in: ", t_j1)

    # Create an instance of the C_profiles class
    C_1 = C_profiles(x_i=x_i, Cb=Cb)
    C_2 = C_profiles(x_i=x_i, Cb=Cb)
    Cth_profile = C_profiles().create_empty_profile(x_i, Cth)

    # Run the simulation
    print("Running predep simulation...")
    Cp_1,xjunc_1 = preDep.lumerical_on_budget(C_1, Cb, Cth, t_j=t_j0)               #predep.
    print("Junction depth for predep.: ", xjunc_1)

    print("Predep simulation completed. Saving profile...")
    C_2.Cold = Cp_1                                         #set the initial profile for the next simulation

    print("Running drive-in simulation...")
    Cp_2,xjunc_2 = driveIn.lumerical_on_budget(C_2, Cb, Cth, t_j=t_j1, process=1)   #drive-in
    print("Drive-in simulation completed. Saving profile...")
    print("Junction depth for drive-in: ", xjunc_2)

    # Cut the first element of the profile to avoid the initial condition
    Cp_1.cut_initial()
    Cp_2.cut_initial()

    # Create an instance of the plot class
    plotter = plot(Cp_1=Cp_1, Cp_2=Cp_2, Cth_profile=Cth_profile, xJunc_1=xjunc_1, xJunc_2=xjunc_2)

    # Plot the results
    plotter.plot_all()

    # Return the result
    print("Exiting...")
    # return Cp_1, Cp_2


if __name__ == "__main__": #RUN BOY RUN!
    main(0, 1e15, 1000, 1150, 3e-4, 600, 30)
    sys.exit()


# TODO: (Optional) Try to handle the boundary conditions in a more elegant way (indx 0) for both processes