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

def read_elements(beamline: dict,
                  z:np.arange) -> interpolate.interp1d:
    '''Sews elements into a function of z.

    Sews elements into a function of z with parameters:
    beamline --- set of accelerator elements,
    z [m] --- coordinate.
    '''
    field_files = {}
    F = 0
    F_dev = 0
    if not beamline:
        z_data = [i/1000 for i in range(1000)]
        F_data = [0 for i in range(1000)]
        f = interpolate.interp1d(
            z_data, F_data,
            fill_value=(0, 0), bounds_error=False
        )
        F = F + f(z)
        F_dev = F
    else:
        for element in beamline.values():
            if not (element.file_name in field_files):
                field_files[element.file_name] = np.loadtxt(element.file_name)
            M = field_files[element.file_name]
            z_data = M[:,0]
            F_data = M[:,1]
            dz = 2*(z_data[-1] - z_data[0])/len(z_data)

            f = interpolate.interp1d(
                element.z0+z_data, element.max_field*F_data, kind='cubic',
                fill_value=(0, 0), bounds_error=False
            )
            F = F + f(z)

            z_data_dev = z_data
            F_data_dev = np.gradient(F_data, dz)
            f_dev = interpolate.interp1d(
                element.z0+z_data_dev, -element.max_field*F_data_dev, kind='cubic',
                fill_value=(0, 0), bounds_error=False
            )
            F_dev = F_dev + f_dev(z)
    F = interpolate.interp1d(z, F, kind='cubic', fill_value=(0, 0), bounds_error=False)
    F_dev = interpolate.interp1d(z, F_dev, kind='cubic', fill_value=(0, 0), bounds_error=False)
    return F, F_dev


class Accelerator:
    '''Create an accelerator.

    Create an accelerator with parameters:
    z_start [m] --- the beginning of the accelerator,
    z_stop [m] --- end of the accelerator,
    dz [m] --- step method along the accelerator.

    Can add a solenoids, quadrupoles and accelerating modules.

    Accelerator's parameters after compile:
    beamline:
    Bz_beamline, Ez_beamline, Gz_beamline
    function:
    Ez, Bz, Gz
    and
    dEzdz, dBzdz, dGzdz

    '''
    Bz_beamline = {}
    Ez_beamline = {}
    Gz_beamline = {}
    Bz = interpolate.interp1d
    Ez = interpolate.interp1d
    Gz = interpolate.interp1d
    dBzdz = interpolate.interp1d
    dEzdz = interpolate.interp1d
    dGzdz = interpolate.interp1d

    def __init__(self,
                 z_start: float,
                 z_stop: float,
                 dz: float):
        self.z_start = self.start = z_start
        self.z_stop = self.stop = z_stop
        self.dz = self.step = dz

        self.z = self.parameter = np.arange(z_start, z_stop, dz)

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
        self.Bz, self.dBzdz = read_elements(self.Bz_beamline, self.parameter)
        self.Ez, self.dEzdz = read_elements(self.Ez_beamline, self.parameter)
        self.Gz, self.dGzdz = read_elements(self.Gz_beamline, self.parameter)



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
