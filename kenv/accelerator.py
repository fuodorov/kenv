# accelerator.py

'''Create an accelerator.
'''

import numpy as np
from scipy import interpolate

__all__ = ['Element',
           'Accelerator',
           'read_elements']


class Element:
    '''Sets one accelerator element.

    Sets one accelerator element with parameters:
    z0 [m] --- element position,
    max_field --- maximum field, 
    file_name --- field profile,
    name --- unique element name.
    '''
    def __init__(self,
                 z0: float,
                 max_field: float,
                 file_name: str,
                 name: str):
        self.z0 = z0
        self.max_field = max_field
        self.file_name = file_name
        self.name = name


field_files = {} # buffer for field files

def read_elements(z:np.arange,
                  beamline: dict) -> interpolate.interp1d:
    '''Sews elements into a function of z.

    Sews elements into a function of z with parameters:
    beamline --- set of accelerator elements,
    z [m] --- coordinate.
    '''
    global field_files
    F = 0
    if not beamline:
        z_data = [i/1000 for i in range(1000)]
        F_data = [0 for i in range(1000)]
        f = interpolate.interp1d(
            z_data, F_data,
            fill_value=(0, 0), bounds_error=False
        )
        F = F + f(z)
    else:
        for element in beamline.values():
            if not (element.file_name in field_files):
                field_files[element.file_name] = np.loadtxt(element.file_name)
            M = field_files[element.file_name]
            z_data = M[:,0]
            F_data = M[:,1]
            f = interpolate.interp1d(
                element.z0+z_data, element.max_field*F_data,
                fill_value=(0, 0), bounds_error=False
            )
            F = F + f(z)
    F = interpolate.interp1d(z, F, fill_value=(0, 0), bounds_error=False)
    return F

        
class Accelerator:
    '''Create an accelerator.
    
    Create an accelerator with parameters:
    start [m] --- the beginning of the accelerator,
    stop [m] --- end of the accelerator,
    step [m] --- step method along the accelerator.
    
    Can add a solenoids, quadrupoles and accelerating modules.
    '''
    def __init__(self, 
                 start: float,
                 stop: float,
                 step: float):
        self.start = start
        self.stop = stop
        self.step = step
        self.parameter = np.arange(start, stop, step)
        self.Bz_beamline = {}
        self.Ez_beamline = {}
        self.Gz_beamline = {}
        self.Bz = interpolate.interp1d
        self.Ez = interpolate.interp1d
        self.Gz = interpolate.interp1d
    
    def add_solenoid(self,
                     name: str,
                     center: float,
                     max_field: float,
                     file_name: str ) -> None:
        '''
        Creates a solenoid in the accelerator.
        
        Creates a solenoid in the accelerator with parameters:
        name --- solenoid's id,
        centert [m] --- solenoid's center,
        max field [T] --- solenoid's maximum field,
        file name --- experimental profile of the Bz field,
        '''
        self.Bz_beamline[name] = Element(center, max_field, file_name, name)
    
    def add_accel(self,
                     name: str,
                     center: float,
                     max_field: float,
                     file_name: str ) -> None:
        '''
        Creates an accelerating module in the accelerator.
        
        Creates an accelerating module in the accelerator with parameters:
        name --- accelerating module's id,
        centert [m] --- accelerating module's center,
        max field [MV/m] --- accelerating module's maximum field,
        file name --- experimental profile of the Ez field,
        '''
        self.Ez_beamline[name] = Element(center, max_field, file_name, name)
    
    def add_quadrupole(self,
                     name: str,
                     center: float,
                     max_field: float,
                     file_name: str ) -> None:
        '''
        Creates a quadrupole in the accelerator.
        
        Creates a quadrupole in the accelerator with parameters:
        name --- quadrupole's id,
        centert [m] --- quadrupole's center,
        max field [T/m] --- quadrupole's maximum field,
        file name --- experimental profile of the Gz field,
        '''
        self.Gz_beamline[name] = Element(center, max_field, file_name, name)
    
    def delete_solenoid(self,
                     name: str='all') -> None:
        '''
        Delete a solenoid in the accelerator.
        
        Delete a solenoid in the accelerator with parameters:
        name --- solenoid's id
        
        *if there is no name that will be removed all
        '''
        if name == 'all':
            self.Bz_beamline = {}
        else:
            self.Bz_beamline.pop(name)
    
    def delete_accel(self,
                     name: str='all') -> None:
        '''
        Delete a accelerating module in the accelerator.
        
        Delete a accelerating module in the accelerator with parameters:
        name --- quadrupole's id
        
        *if there is no name that will be removed all
        '''
        if name == 'all':
            self.Ez_beamline = {}
        else:
            self.Ez_beamline.pop(name)
    
    def delete_quadrupole(self,
                     name: str='all') -> None:
        '''
        Delete a quadrupole in the accelerator.
        
        Delete a quadrupole in the accelerator with parameters:
        name --- quadrupole's id
        
        *if there is no name that will be removed all
        '''
        if name == 'all':
            self.Gz_beamline = {}
        else:
            self.Gz_beamline.pop(name)
    
    
    def compile(self) -> None:
        '''Compilation of the accelerator.
        '''
        self.Bz = read_elements(self.parameter, self.Bz_beamline)
        self.Ez = read_elements(self.parameter, self.Ez_beamline)
        self.Gz = read_elements(self.parameter, self.Gz_beamline)
        
        print('Accelerator compiled.')
     
    def __str__(self):
        string = 'Accelerator structure.\n'
        string += '\tSolenoids:\n'
        for element in self.Bz_beamline.values():
            string +="\t[ %.5f m, %.5f T, '%s', '%s'],\n" % (element.z0, element.max_field, element.file_name, element.name)
        string += '\tAccelerating modules:\n'
        for element in self.Ez_beamline.values():
            string +="\t[ %.5f m, %.5f Mv/m, '%s', '%s'],\n" % (element.z0, element.max_field, element.file_name, element.name)
        string += '\tQuadrupoles:\n'
        for element in self.Gz_beamline.values():
            string +="\t[ %.5f m, %.5f T/m, '%s', '%s'],\n" % (element.z0, element.max_field, element.file_name, element.name)
        return string
    
    add_sol = add_new_solenoid = add_solenoid
    add_acc = add_new_accel = add_accel
    add_quad = add_new_quadrupole = add_quadrupole
    
    del_sol = del_solenoid = delete_solenoid
    del_acc = del_accel = delete_accel
    del_quad = del_quadrupole = delete_quadrupole
