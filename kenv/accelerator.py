# accelerator.py
'''Create an accelerator.

'''

import numpy as np
from scipy import interpolate, integrate, misc

__all__ = ['Element',
           'Accelerator',
           'read_fields',
           'read_offsets']


class Element:
    '''Sets one accelerator element.

    Sets one accelerator element with parameters:
    z0 [m] --- element position,
    max_field --- maximum field,
    file_name --- field profile,
    name --- unique element name,

    and shifted by x, xp, y, yp,

    and z_start, z_stop, length

    '''

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

        self.z_start = .0e0
        self.z_stop = .0e0
        self.length = .0e0
        self.field = .0e0

def read_fields(beamline: dict,
                z: np.arange) -> interpolate.interp1d:
    '''Sews elements into a function of z.

    Sews elements into a function of z with parameters:
    beamline --- set of accelerator elements,
    z [m] --- coordinate.

    '''

    field_files = {}
    F, F_prime, F_int = 0, 0, 0

    if not beamline:
        z_data = [i/len(z) for i in range(len(z))]
        F_data = [0 for i in range(len(z))]
        f = interpolate.interp1d(
            z_data, F_data,
            fill_value=(0, 0), bounds_error=False
        )
        F = F + f(z)
        F_prime = F
        F_int = F
    else:
        for element in beamline.values():
            if not (element.file_name in field_files):
                field_files[element.file_name] = np.loadtxt(element.file_name)
            M = field_files[element.file_name]
            z_data = M[:,0] + element.z0
            F_data = M[:,1]
            dz = (z_data[-1] - z_data[0])/len(z_data)

            f = interpolate.interp1d(
                z_data, element.max_field*F_data, kind='cubic',
                fill_value=(0, 0), bounds_error=False
            )

            f_prime = interpolate.interp1d(
                z_data, element.max_field*np.gradient(F_data, dz), kind='cubic',
                fill_value=(0, 0), bounds_error=False
            )

            F = F + f(z)
            F_prime = F_prime + f_prime(z)

            if not element.max_field == 0:
                element.length = ((integrate.cumtrapz(f(z_data), z_data)[-1])**2)/((integrate.cumtrapz(f(z_data)**2, z_data))[-1])
                element.field = (integrate.cumtrapz(f(z_data)**2, z_data)[-1])/(integrate.cumtrapz(f(z_data), z_data)[-1])
                element.z_start = element.z0 - element.length
                element.z_stop = element.z0 + element.length
            else:
                element.length = 0
                element.field = 0
                element.z_start = element.z0
                element.z_stop = element.z0

    F = interpolate.interp1d(z, F, kind='cubic', fill_value=(0, 0), bounds_error=False)
    F_int = integrate.cumtrapz(F(z), z)
    F_prime = interpolate.interp1d(z, F_prime, kind='cubic', fill_value=(0, 0), bounds_error=False)
    F_int = interpolate.interp1d(z[1:], F_int, kind='cubic', fill_value=(F_int[0], F_int[-1]), bounds_error=False)

    return F, F_prime, F_int

def read_offsets(beamline: dict,
                z: np.arange) -> interpolate.interp1d:
    '''Sews elements into a function of z.

    Sews elements into a function of z with parameters:
    beamline --- set of accelerator elements,
    z [m] --- coordinate.

    '''

    field_files = {}
    F = 0
    offset_correct_x, offset_correct_y, offset_correct_xp, offset_correct_yp  = 0, 0, 0, 0

    if not beamline:
        z_data = [i/len(z) for i in range(len(z))]
        F_data = [0 for i in range(len(z))]
        f = interpolate.interp1d(
            z_data, F_data,
            fill_value=(0, 0), bounds_error=False
        )
        F = F + f(z)
        offset_correct_x, offset_correct_y, offset_correct_xp, offset_correct_yp = F, F, F, F
    else:
        for element in beamline.values():
            if not (element.file_name in field_files):
                field_files[element.file_name] = np.loadtxt(element.file_name)
            M = field_files[element.file_name]
            z_data = M[:,0] + element.z0
            F_data = M[:,1]

            if not element.length == 0:
                z_data = np.linspace(element.z_start, element.z_stop, len(z_data))

                f_x = interpolate.interp1d(
                    z_data, [element.x + (z_data[i]-element.z0)*element.xp for i in range(len(z_data))],
                    fill_value=(0,0), bounds_error=False
                )
                f_xp = interpolate.interp1d(
                    z_data, [element.xp for i in range(len(z_data))],
                    fill_value=(0, 0), bounds_error=False
                )
                f_y = interpolate.interp1d(
                    z_data, [element.y + (z_data[i]-element.z0)*element.yp for i in range(len(z_data))],
                    fill_value=(0, 0), bounds_error=False
                )
                f_yp = interpolate.interp1d(
                    z_data, [element.yp for i in range(len(z_data))],
                    fill_value=(0, 0), bounds_error=False
                )
            else:
                f = interpolate.interp1d(
                    z_data, 0*F_data, kind='cubic',
                    fill_value=(0, 0), bounds_error=False
                )
                f_x, f_xp, f_y, f_yp = f, f, f, f

            offset_correct_x = offset_correct_x + f_x(z)
            offset_correct_xp = offset_correct_xp + f_xp(z)
            offset_correct_y = offset_correct_y + f_y(z)
            offset_correct_yp = offset_correct_yp + f_yp(z)

    offset_correct_x = interpolate.interp1d(z, offset_correct_x, kind='linear', fill_value=(0, 0), bounds_error=False)
    offset_correct_y = interpolate.interp1d(z, offset_correct_y, kind='linear', fill_value=(0, 0), bounds_error=False)
    offset_correct_xp = interpolate.interp1d(z, offset_correct_xp, kind='linear', fill_value=(0, 0), bounds_error=False)
    offset_correct_yp = interpolate.interp1d(z, offset_correct_yp, kind='linear', fill_value=(0, 0), bounds_error=False)

    return offset_correct_x, offset_correct_xp, offset_correct_y, offset_correct_yp

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
    Dx, Dxp, Dy, Dyp,
    and
    dEzdz, dBzdz, dGzdz, dBxdz, dBydz,
    and
    Ezdz, Bzdz, Gzdz, Bxdz, Bydz

    '''

    def __init__(self,
                 z_start: float,
                 z_stop: float,
                 dz: float):
        self.z_start = self.start = z_start
        self.z_stop = self.stop = z_stop
        self.dz = self.step = dz
        self.z = self.parameter = np.arange(z_start, z_stop, dz)

        self.Bx_beamline, self.By_beamline, self.Bz_beamline = {}, {}, {}
        self.Ez_beamline = {}
        self.Gz_beamline = {}
        self.Bx, self.By, self.Bz = interpolate.interp1d, interpolate.interp1d, interpolate.interp1d
        self.Ez = interpolate.interp1d
        self.Gz = interpolate.interp1d
        self.dBxdz, self.dBydz, self.dBzdz = interpolate.interp1d, interpolate.interp1d, interpolate.interp1d
        self.dEzdz = interpolate.interp1d
        self.dGzdz = interpolate.interp1d
        self.Bxdz, self.Bydz, self.Bzdz = interpolate.interp1d, interpolate.interp1d, interpolate.interp1d
        self.Ezdz = interpolate.interp1d
        self.Gzdz = interpolate.interp1d
        self.Dx, self.Dxp, self.Dy, self.Dyp = interpolate.interp1d, interpolate.interp1d, interpolate.interp1d, interpolate.interp1d

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
        and shifted by x [m], xp [rad], y [m], yp [rad]

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
        and shifted by x [m], xp [rad], y [m], yp [rad]

        '''
        self.Ez_beamline[name] = Element(center, max_field, file_name, name,
                                         x=x, xp=xp, y=y, yp=yp)

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

        self.Bx, self.dBxdz, self.Bxdz = read_fields(self.Bx_beamline, self.parameter)
        self.By, self.dBydz, self.Bydz = read_fields(self.By_beamline, self.parameter)
        self.Bz, self.dBzdz, self.Bzdz = read_fields(self.Bz_beamline, self.parameter)
        self.Ez, self.dEzdz, self.Ezdz = read_fields(self.Ez_beamline, self.parameter)
        self.Gz, self.dGzdz, self.Gzdz = read_fields(self.Gz_beamline, self.parameter)

        Dx_Bz, Dxp_Bz, Dy_Bz, Dyp_Bz = read_offsets(self.Bz_beamline, self.parameter)
        Dx_Ez, Dxp_Ez, Dy_Ez, Dyp_Ez = read_offsets(self.Ez_beamline, self.parameter)

        self.Dx = interpolate.interp1d(self.z, Dx_Bz(self.z) + Dx_Ez(self.z), kind='linear', fill_value=(0, 0), bounds_error=False)
        self.Dxp = interpolate.interp1d(self.z, Dxp_Bz(self.z) + Dxp_Ez(self.z), kind='linear', fill_value=(0, 0), bounds_error=False)
        self.Dy = interpolate.interp1d(self.z, Dy_Bz(self.z) + Dy_Ez(self.z), kind='linear', fill_value=(0, 0), bounds_error=False)
        self.Dyp = interpolate.interp1d(self.z, Dyp_Bz(self.z) + Dyp_Ez(self.z), kind='linear', fill_value=(0, 0), bounds_error=False)

    def __str__(self):
        string = 'Accelerator structure.\n'
        string += '\tSolenoids:\n'
        for element in self.Bz_beamline.values():
            string +="\t[ %.5f m, %.5f T, '%s', '%s', %.5f m, %.5f rad, %.5f m, %.5f rad],\n"\
             % (element.z0, element.max_field, element.file_name, element.name,
                element.x, element.xp, element.y, element.yp)
        string += '\tAccelerating modules:\n'
        for element in self.Ez_beamline.values():
            string +="\t[ %.5f m, %.5f T, '%s', '%s', %.5f m, %.5f rad, %.5f m, %.5f rad],\n"\
             % (element.z0, element.max_field, element.file_name, element.name,
                element.x, element.xp, element.y, element.yp)
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
