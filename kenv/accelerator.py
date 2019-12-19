# accelerator.py
'''Create an accelerator.

'''

import numpy as np
from scipy import interpolate, integrate

__all__ = ['Element',
           'Accelerator',
           'read_elements']


class Element:
    '''Sets one accelerator element.

    Sets one accelerator element with parameters:
    z0 [m] --- element position,
    max_field --- maximum field,
    file_name --- field profile,
    name --- unique element name,

    and shifted to [x, xp, y, yp],

    and z_start, z_stop, length

    '''
    z_start = .0e0
    z_stop = .0e0
    length = .0e0
    field = .0e0
    def __init__(self,
                 z0: float,
                 max_field: float,
                 file_name: str,
                 name: str,
                 *,
                 x: float=0.0,
                 xp: float=0.0,
                 y: float=0.0,
                 yp: float=0.0):
        self.z0 = z0
        self.max_field = max_field
        self.file_name = file_name
        self.name = name
        self.x = x
        self.xp = xp
        self.y = y
        self.yp = yp

def read_elements(beamline: dict,
                  z:np.arange,*, n=1000) -> interpolate.interp1d:
    '''Sews elements into a function of z.

    Sews elements into a function of z with parameters:
    beamline --- set of accelerator elements,
    z [m] --- coordinate.

    '''
    field_files = {}
    F = 0
    F_prime = 0
    F_int = 0
    offset_correct_x = 0
    offset_correct_xp = 0
    offset_correct_y = 0
    offset_correct_yp = 0
    if not beamline:
        z_data = [i/n for i in range(n)]
        F_data = [0 for i in range(n)]
        f = interpolate.interp1d(
            z_data, F_data,
            fill_value=(0, 0), bounds_error=False
        )
        F = F + f(z)
        F_prime = F
        F_int = F
        offset_correct_x = F
        offset_correct_xp = F
        offset_correct_y = F
        offset_correct_yp = F
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

            element.length = (integrate.cumtrapz(f(z), z)**2)[-1]/(integrate.cumtrapz(f(z)**2, z))[-1]
            element.field = (integrate.cumtrapz(f(z)**2, z)[-1])/(integrate.cumtrapz(f(z), z)[-1])
            element.z_start = element.z0 - element.length/2
            element.z_stop = element.z0 + element.length/2
            #derivative
            z_data_prime = z_data
            F_data_prime = np.gradient(F_data, dz)
            f_prime = interpolate.interp1d(
                element.z0+z_data_prime, element.max_field*F_data_prime, kind='cubic',
                fill_value=(0, 0), bounds_error=False
            )
            F_prime = F_prime + f_prime(z)

            #offset correction
            z_data = np.linspace(element.z_start, element.z_stop, n)
            f_x = interpolate.interp1d(
                z_data, [element.x for i in range(n)],
                fill_value=(0,0), bounds_error=False
            )
            f_xp = interpolate.interp1d(
                z_data, [element.xp for i in range(n)],
                fill_value=(0, 0), bounds_error=False
            )
            f_y = interpolate.interp1d(
                z_data, [element.y for i in range(n)],
                fill_value=(0, 0), bounds_error=False
            )
            f_yp = interpolate.interp1d(
                z_data, [element.yp for i in range(n)],
                fill_value=(0, 0), bounds_error=False
            )

            offset_correct_x = offset_correct_x + f_x(z)
            offset_correct_xp = offset_correct_xp + f_xp(z)
            offset_correct_y = offset_correct_y + f_y(z)
            offset_correct_yp = offset_correct_yp + f_yp(z)

    F = interpolate.interp1d(z, F, kind='cubic', fill_value=(0, 0), bounds_error=False)
    F_prime = interpolate.interp1d(z, F_prime, kind='cubic', fill_value=(0, 0), bounds_error=False)

    F_int = integrate.cumtrapz(F(z), z)
    F_int = interpolate.interp1d(z[1:], F_int, kind='cubic', fill_value=(F_int[0], F_int[-1]), bounds_error=False)

    offset_correct_x = interpolate.interp1d(z, offset_correct_x, kind='linear', fill_value=(0, 0), bounds_error=False)
    offset_correct_y = interpolate.interp1d(z, offset_correct_y, kind='linear', fill_value=(0, 0), bounds_error=False)
    offset_correct_xp = interpolate.interp1d(z, offset_correct_xp, kind='linear', fill_value=(0, 0), bounds_error=False)
    offset_correct_yp = interpolate.interp1d(z, offset_correct_yp, kind='linear', fill_value=(0, 0), bounds_error=False)

    return F, F_prime, F_int, offset_correct_x, offset_correct_xp, offset_correct_y, offset_correct_yp

class Accelerator:
    '''Create an accelerator.

    Create an accelerator with parameters:
    z_start [m] --- the beginning of the accelerator,
    z_stop [m] --- end of the accelerator,
    dz [m] --- step method along the accelerator.

    Can add a solenoids, quadrupoles and accelerating modules.

    Accelerator's parameters after compile:
    beamline:
    Bz_beamline, Ez_beamline, Gz_beamline,
    Bx_beamline, By_beamline,
    function:
    Ez, Bz, Gz, Bx, By,
    and
    dEzdz, dBzdz, dGzdz, dBxdz, dBydz,
    and
    Ezdz, Bzdz, Gzdz, Bxdz, Bydz

    '''
    Bz_beamline = {}
    Bx_beamline = {}
    By_beamline = {}
    Ez_beamline = {}
    Gz_beamline = {}
    Bz = interpolate.interp1d
    Bx = interpolate.interp1d
    By = interpolate.interp1d
    Ez = interpolate.interp1d
    Gz = interpolate.interp1d
    dBzdz = interpolate.interp1d
    dBxdz = interpolate.interp1d
    dBydz = interpolate.interp1d
    dEzdz = interpolate.interp1d
    dGzdz = interpolate.interp1d
    Bzdz = interpolate.interp1d
    Bxdz = interpolate.interp1d
    Bydz = interpolate.interp1d
    Ezdz = interpolate.interp1d
    Gzdz = interpolate.interp1d
    Dx = interpolate.interp1d
    Dxp = interpolate.interp1d
    Dy = interpolate.interp1d
    Dyp = interpolate.interp1d

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
                     file_name: str,
                     *,
                     x: float=0.0,
                     xp: float=0.0,
                     y: float=0.0,
                     yp: float=0.0) -> None:
        '''Creates a solenoid in the accelerator.

        Creates a solenoid in the accelerator with parameters:
        name --- solenoid's id,
        center [m] --- solenoid's center,
        max_field [T] --- solenoid's maximum field,
        file_name --- experimental profile of the Bz field

        '''
        self.Bz_beamline[name] = Element(center, max_field, file_name, name,
                                         x=x, xp=xp, y=y, yp=yp)

    def add_accel(self,
                  name: str,
                  center: float,
                  max_field: float,
                  file_name: str,
                  *,
                  x: float=0.0,
                  xp: float=0.0,
                  y: float=0.0,
                  yp: float=0.0) -> None:
        '''Creates an accelerating module in the accelerator.

        Creates an accelerating module in the accelerator with parameters:
        name --- accelerating module's id,
        center [m] --- accelerating module's center,
        max_field [MV/m] --- accelerating module's maximum field,
        file_name --- experimental profile of the Ez field

        '''
        self.Ez_beamline[name] = Element(center, max_field, file_name, name)

    def add_quadrupole(self,
                       name: str,
                       center: float,
                       max_field: float,
                       file_name: str) -> None:
        '''Creates a quadrupole in the accelerator.

        Creates a quadrupole in the accelerator with parameters:
        name --- quadrupole's id,
        center [m] --- quadrupole's center,
        max_field [T/m] --- quadrupole's maximum field,
        file_name --- experimental profile of the Gz field

        '''
        self.Gz_beamline[name] = Element(center, max_field, file_name, name)

    def add_corrector_x(self,
                        name: str,
                        center: float,
                        max_field: float,
                        file_name: str) -> None:
        '''Creates a corrector in the accelerator.

        Creates a corrector in the accelerator with parameters:
        name --- corrector's id,
        center [m] --- corrector's center,
        max_field [T] --- corrector's maximum field,
        file_name --- experimental profile of the By field

        '''
        self.By_beamline[name] = Element(center, max_field, file_name, name)

    def add_corrector_y(self,
                        name: str,
                        center: float,
                        max_field: float,
                        file_name: str) -> None:
        '''Creates a corrector in the accelerator.

        Creates a corrector in the accelerator with parameters:
        name --- corrector's id,
        center [m] --- corrector's center,
        max_field [T] --- corrector's maximum field,
        file_name --- experimental profile of the Bx field

        '''
        self.Bx_beamline[name] = Element(center, max_field, file_name, name)

    def delete_solenoid(self,
                        name: str='all') -> None:
        '''Delete a solenoid in the accelerator.

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
        '''Delete a accelerating module in the accelerator.

        Delete a accelerating module in the accelerator with parameters:
        name --- accel's id

        *if there is no name that will be removed all

        '''
        if name == 'all':
            self.Ez_beamline = {}
        else:
            self.Ez_beamline.pop(name)

    def delete_quadrupole(self,
                          name: str='all') -> None:
        '''Delete a quadrupole in the accelerator.

        Delete a quadrupole in the accelerator with parameters:
        name --- quadrupole's id

        *if there is no name that will be removed all

        '''
        if name == 'all':
            self.Gz_beamline = {}
        else:
            self.Gz_beamline.pop(name)

    def delete_corrector_x(self,
                           name: str='all') -> None:
        '''Delete a corrector in the accelerator.

        Delete a corrector in the accelerator with parameters:
        name --- corrector's id

        *if there is no name that will be removed all

        '''
        if name == 'all':
            self.By_beamline = {}
        else:
            self.By_beamline.pop(name)

    def delete_corrector_y(self,
                           name: str='all') -> None:
        '''Delete a corrector in the accelerator.

        Delete a corrector in the accelerator with parameters:
        name --- corrector's id

        *if there is no name that will be removed all

        '''
        if name == 'all':
            self.Bx_beamline = {}
        else:
            self.Bx_beamline.pop(name)

    def compile(self) -> None:
        '''Compilation of the accelerator.'''

        zero_box = 0
        self.Bz, self.dBzdz, self.Bzdz, self.Dx, self.Dxp, self.Dy, self.Dyp  = read_elements(self.Bz_beamline, self.parameter)
        self.Bx, self.dBxdz, self.Bxdz, *zero_box = read_elements(self.Bx_beamline, self.parameter)
        self.By, self.dBydz, self.Bydz, *zero_box = read_elements(self.By_beamline, self.parameter)
        self.Ez, self.dEzdz, self.Ezdz, *zero_box = read_elements(self.Ez_beamline, self.parameter)
        self.Gz, self.dGzdz, self.Gzdz, *zero_box = read_elements(self.Gz_beamline, self.parameter)

    def __str__(self):
        string = 'Accelerator structure.\n'
        string += '\tSolenoids:\n'
        for element in self.Bz_beamline.values():
            string +="\t[ %.5f m, %.5f T, '%s', '%s', %.5f m, %.5f rad, %.5f m, %.5f rad],\n"\
             % (element.z0, element.max_field, element.file_name, element.name,
                element.x, element.xp, element.y, element.yp)
        string += '\tAccelerating modules:\n'
        for element in self.Ez_beamline.values():
            string +="\t[ %.5f m, %.5f T, '%s', '%s'],\n"\
             % (element.z0, element.max_field, element.file_name, element.name)
        string += '\tQuadrupoles:\n'
        for element in self.Gz_beamline.values():
            string +="\t[ %.5f m, %.5f T, '%s', '%s'],\n"\
             % (element.z0, element.max_field, element.file_name, element.name)
        string += '\tCorrectors x:\n'
        for element in self.By_beamline.values():
            string +="\t[ %.5f m, %.5f T, '%s', '%s'],\n"\
             % (element.z0, element.max_field, element.file_name, element.name)
        string += '\tCorrectors y:\n'
        for element in self.Bx_beamline.values():
            string +="\t[ %.5f m, %.5f T, '%s', '%s'],\n"\
             % (element.z0, element.max_field, element.file_name, element.name)
        return string

    add_sol = add_new_solenoid = add_solenoid
    add_acc = add_new_accel = add_accel
    add_quad = add_new_quadrupole = add_quadrupole
    add_corr_x = add_new_corrector_x = add_corrector_x
    add_corr_y = add_new_corrector_y = add_corrector_y

    del_sol = del_solenoid = delete_solenoid
    del_acc = del_accel = delete_accel
    del_quad = del_quadrupole = delete_quadrupole
    del_corr_x = del_corrector_x = delete_corrector_x
    del_corr_y = del_corrector_y = delete_corrector_y
