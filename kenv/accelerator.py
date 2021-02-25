import numpy as np
from scipy import interpolate, integrate

__all__ = ['Element', 'Accelerator', 'read_fields', 'read_offsets']


class Element:
    """
    Sets one accelerator element.

    """

    def __init__(self,
                 z0: float,
                 max_field: float,
                 file_name: str,
                 name: str,
                 *,
                 x: float = 0.0,
                 xp: float = 0.0,
                 y: float = 0.0,
                 yp: float = 0.0):
        """
        Initialization of an accelerator element.

        Parameters
        ----------
        z0: float
            Position of the element center
        max_field: float
            Maximum value of the field on the element axis
        file_name: string
            Name of the file
        name: string
            Name of the element

        x: float, optional
            Offset of the element along the x-axis, [m]
        xp: float, optional
            Element rotation in the z-x plane, [rad]
        y: float, optional
            Offset of the element along the y-axis, [m]
        yp: float, optional
            Element rotation in the z-y plane, [rad]
        """
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
    """
    Field reading function.

    Parameters
    ----------
    beamline: dict
        Accelerator`s beamline
    z: an instance of the `np.arange` class
        Axis z

    Returns
    -------
    F: an instance of the `interpolate.interp1d` class
        Field function
    F_prime: an instance of the `interpolate.interp1d` class
        Field prime function
    F_int: an instance of the `interpolate.interp1d` class
        Field integrate function
    """

    field_files = {}
    F, F_prime, F_int = 0, 0, 0

    if not beamline:
        z_data = [i / len(z) for i in range(len(z))]
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
            z_data = M[:, 0] + element.z0
            F_data = M[:, 1]
            dz = (z_data[-1] - z_data[0]) / len(z_data)

            f = interpolate.interp1d(
                z_data,
                element.max_field * F_data,
                kind='cubic',
                fill_value=(0, 0),
                bounds_error=False
            )

            f_prime = interpolate.interp1d(
                z_data,
                element.max_field * np.gradient(F_data, dz),
                kind='cubic',
                fill_value=(0, 0),
                bounds_error=False
            )

            F = F + f(z)
            F_prime = F_prime + f_prime(z)

            if not element.max_field == 0:
                int_f = integrate.cumtrapz(f(z_data), z_data)[-1]
                int_ff = integrate.cumtrapz(f(z_data)**2, z_data)[-1]
                element.length = int_f**2 / int_ff
                element.field = int_ff / int_f
                element.z_start = element.z0 - element.length
                element.z_stop = element.z0 + element.length
            else:
                element.length = 0
                element.field = 0
                element.z_start = element.z0
                element.z_stop = element.z0

    F = interpolate.interp1d(z, F, kind='cubic',
                             fill_value=(0, 0), bounds_error=False)
    F_int = integrate.cumtrapz(F(z), z)
    F_prime = interpolate.interp1d(z, F_prime, kind='cubic',
                                   fill_value=(0, 0), bounds_error=False)
    F_int = interpolate.interp1d(z[1:], F_int, kind='cubic',
                                 fill_value=(F_int[0], F_int[-1]),
                                 bounds_error=False)

    return F, F_prime, F_int


def read_offsets(beamline: dict,
                 z: np.arange) -> interpolate.interp1d:
    """
    Offsets and rotations reading function.

    Parameters
    ----------
    beamline: dict
        Accelerator`s beamline
    z: an instance of the `np.arange` class
        Axis z

    Returns
    -------
    offset_correct_x: an instance of the `interpolate.interp1d` class
        Offset function
    offset_correct_xp: an instance of the `interpolate.interp1d` class
        Rotation function
    offset_correct_y: an instance of the `interpolate.interp1d` class
        Offset function
    offset_correct_yp: an instance of the `interpolate.interp1d` class
        Rotation function
    """

    field_files = {}
    F = 0
    offset_correct_x, offset_correct_y = 0, 0
    offset_correct_xp, offset_correct_yp = 0, 0

    if not beamline:
        z_data = [i / len(z) for i in range(len(z))]
        F_data = [0 for i in range(len(z))]
        f = interpolate.interp1d(
            z_data, F_data,
            fill_value=(0, 0), bounds_error=False
        )
        F = F + f(z)
        offset_correct_x, offset_correct_y = F, F
        offset_correct_xp, offset_correct_yp = F, F
    else:
        for element in beamline.values():
            if not (element.file_name in field_files):
                field_files[element.file_name] = np.loadtxt(element.file_name)
            M = field_files[element.file_name]
            z_data = M[:, 0] + element.z0
            F_data = M[:, 1]

            if not element.length == 0:
                z_data = np.linspace(element.z_start, element.z_stop,
                                     len(z_data))
                x, xp = element.x, element.xp
                y, yp = element.y, element.yp
                z0 = element.z0
                f_x = interpolate.interp1d(
                    z_data,
                    [x + (z_data[i] - z0) * xp for i in range(len(z_data))],
                    fill_value=(0, 0),
                    bounds_error=False
                )
                f_xp = interpolate.interp1d(
                    z_data,
                    [xp for i in range(len(z_data))],
                    fill_value=(0, 0),
                    bounds_error=False
                )
                f_y = interpolate.interp1d(
                    z_data,
                    [y + (z_data[i] - z0) * yp for i in range(len(z_data))],
                    fill_value=(0, 0),
                    bounds_error=False
                )
                f_yp = interpolate.interp1d(
                    z_data,
                    [yp for i in range(len(z_data))],
                    fill_value=(0, 0),
                    bounds_error=False
                )
            else:
                f = interpolate.interp1d(
                    z_data,
                    0 * F_data,
                    kind='cubic',
                    fill_value=(0, 0),
                    bounds_error=False
                )
                f_x, f_xp, f_y, f_yp = f, f, f, f

            offset_correct_x = offset_correct_x + f_x(z)
            offset_correct_xp = offset_correct_xp + f_xp(z)
            offset_correct_y = offset_correct_y + f_y(z)
            offset_correct_yp = offset_correct_yp + f_yp(z)

    offset_correct_x = interpolate.interp1d(
        z,
        offset_correct_x,
        kind='linear',
        fill_value=(0, 0),
        bounds_error=False
    )
    offset_correct_y = interpolate.interp1d(
        z,
        offset_correct_y,
        kind='linear',
        fill_value=(0, 0),
        bounds_error=False
    )
    offset_correct_xp = interpolate.interp1d(
        z,
        offset_correct_xp,
        kind='linear',
        fill_value=(0, 0),
        bounds_error=False
    )
    offset_correct_yp = interpolate.interp1d(
        z,
        offset_correct_yp,
        kind='linear',
        fill_value=(0, 0),
        bounds_error=False
    )
    return (offset_correct_x, offset_correct_xp,
            offset_correct_y, offset_correct_yp)


class Accelerator:
    """
    Top-level accelerator class that contains all the accelerator data.

    The `Accelerator` class has several important attributes:

    - `Bx_beamline`, a beamline which contains the solenoids information
    - `By_beamline`, a beamline which contains the solenoids information
    - `Bz_beamline`, a beamline which contains the solenoids information
    - `Ez_beamline`, a beamline which contains the accels information
    - `Gz_beamline`, a beamline which contains the quadrupoles information
    """

    def __init__(self,
                 z_start: float,
                 z_stop: float,
                 dz: float):
        """
        Initializes an accelerator.

        Parameters
        ----------
        z_start: float
            The position of the start of the accelerator
        z_stop: float
            The position of the end of the accelerator
        dz: float
            The step of the accelerator
        """
        self.z_start = self.start = z_start
        self.z_stop = self.stop = z_stop
        self.dz = self.step = dz
        self.z = self.parameter = np.arange(z_start, z_stop, dz)
        self.Bx_beamline, self.By_beamline, self.Bz_beamline = {}, {}, {}
        self.Ez_beamline = {}
        self.Gz_beamline = {}
        self.Bx, self.By = interpolate.interp1d, interpolate.interp1d
        self.Bz = interpolate.interp1d
        self.Ez = interpolate.interp1d
        self.Gz = interpolate.interp1d
        self.dBxdz, self.dBydz = interpolate.interp1d, interpolate.interp1d
        self.dBzdz = interpolate.interp1d
        self.dEzdz = interpolate.interp1d
        self.dGzdz = interpolate.interp1d
        self.Bxdz, self.Bydz = interpolate.interp1d, interpolate.interp1d
        self.Bzdz = interpolate.interp1d
        self.Ezdz = interpolate.interp1d
        self.Gzdz = interpolate.interp1d
        self.Dx, self.Dxp = interpolate.interp1d, interpolate.interp1d
        self.Dy, self.Dyp = interpolate.interp1d, interpolate.interp1d

    def add_solenoid(self,
                     name: str,
                     center: float,
                     max_field: float,
                     file_name: str,
                     *,
                     x: float = 0.0,
                     xp: float = 0.0,
                     y: float = 0.0,
                     yp: float = 0.0) -> None:
        """
        Creates a solenoid lenses in the accelerator.

        Parameters
        ----------
        name: string
            Solenoid's id
        center: float
            Solenoid's center
        max_field: float
            Solenoid's maximum field [T]
        file_name: string
            Experimental profile of the Bz field

        x: float, optional
            Offset of the element along the x-axis
        xp: float, optional
            Element rotation in the z-x plane
        y: float, optional
            Offset of the element along the y-axis
        yp: float, optional
            Element rotation in the z-y plane
        """
        self.Bz_beamline[name] = Element(center, max_field, file_name, name,
                                         x=x, xp=xp, y=y, yp=yp)

    def add_accel(self,
                  name: str,
                  center: float,
                  max_field: float,
                  file_name: str,
                  *,
                  x: float = 0.0,
                  xp: float = 0.0,
                  y: float = 0.0,
                  yp: float = 0.0) -> None:
        """
        Creates an accelerating module in the accelerator.

        Parameters
        ----------
        name: string
            Accel's id
        center: float
            Accel's center
        max_field: float
            Accel's maximum field [MV/m]
        file_name: string
            Experimental profile of the Ez field

        x: float, optional
            Offset of the element along the x-axis
        xp: float, optional
            Element rotation in the z-x plane
        y: float, optional
            Offset of the element along the y-axis
        yp: float, optional
            Element rotation in the z-y plane
        """
        self.Ez_beamline[name] = Element(center, max_field, file_name, name,
                                         x=x, xp=xp, y=y, yp=yp)

    def add_quadrupole(self,
                       name: str,
                       center: float,
                       max_field: float,
                       file_name: str) -> None:
        """
        Creates a quadrupole in the accelerator.

        Parameters
        ----------
        name: string
            Quad's id
        center: float
            Quad's center
        max_field: float
            Quad's maximum field [T/m]
        file_name: string
            Experimental profile of the Gz field

        x: float, optional
            Offset of the element along the x-axis
        xp: float, optional
            Element rotation in the z-x plane
        y: float, optional
            Offset of the element along the y-axis
        yp: float, optional
            Element rotation in the z-y plane
        """
        self.Gz_beamline[name] = Element(center, max_field, file_name, name)

    def add_corrector_x(self,
                        name: str,
                        center: float,
                        max_field: float,
                        file_name: str) -> None:
        """
        Creates a corrector x in the accelerator.

        Parameters
        ----------
        name: string
            Corrector's id
        center: float
            Corrector's center
        max_field: float
            Corrector's maximum field [T]
        file_name: string
            Experimental profile of the By field

        x: float, optional
            Offset of the element along the x-axis
        xp: float, optional
            Element rotation in the z-x plane
        y: float, optional
            Offset of the element along the y-axis
        yp: float, optional
            Element rotation in the z-y plane
        """
        self.By_beamline[name] = Element(center, max_field, file_name, name)

    def add_corrector_y(self,
                        name: str,
                        center: float,
                        max_field: float,
                        file_name: str) -> None:
        """
        Creates a corrector y in the accelerator.

        Parameters
        ----------
        name: string
            Corrector's id
        center: float
            Corrector's center
        max_field: float
            Corrector's maximum field [T]
        file_name: string
            Experimental profile of the Bx field

        x: float, optional
            Offset of the element along the x-axis
        xp: float, optional
            Element rotation in the z-x plane
        y: float, optional
            Offset of the element along the y-axis
        yp: float, optional
            Element rotation in the z-y plane
        """
        self.Bx_beamline[name] = Element(center, max_field, file_name, name)

    def delete_solenoid(self, name: str = 'all') -> None:
        """
        Delete a solenoid in the accelerator.

        Parameters
        ----------
        name: string
            Solenoid`s id
            (If there is no name that will be removed all)
        """
        if name == 'all':
            self.Bz_beamline = {}
        else:
            self.Bz_beamline.pop(name)

    def delete_accel(self, name: str = 'all') -> None:
        """
        Delete an accelerating module in the accelerator.

        Parameters
        ----------
        name: string
            Accel`s id
            (If there is no name that will be removed all)
        """
        if name == 'all':
            self.Ez_beamline = {}
        else:
            self.Ez_beamline.pop(name)

    def delete_quadrupole(self, name: str = 'all') -> None:
        """
        Delete a quadrupole in the accelerator.

        Parameters
        ----------
        name: string
            Quad`s id
            (If there is no name that will be removed all)
        """
        if name == 'all':
            self.Gz_beamline = {}
        else:
            self.Gz_beamline.pop(name)

    def delete_corrector_x(self, name: str = 'all') -> None:
        """
        Delete a corrector x in the accelerator.

        Parameters
        ----------
        name: string
            Corrector`s id
            (If there is no name that will be removed all)
        """
        if name == 'all':
            self.By_beamline = {}
        else:
            self.By_beamline.pop(name)

    def delete_corrector_y(self, name: str = 'all') -> None:
        """
        Delete a corrector y in the accelerator.

        Parameters
        ----------
        name: string
           Corrector`s id
           (If there is no name that will be removed all)
        """
        if name == 'all':
            self.Bx_beamline = {}
        else:
            self.Bx_beamline.pop(name)

    def compile(self) -> None:
        """
        Compilation of the accelerator.
        """

        self.Bx, self.dBxdz, self.Bxdz = read_fields(self.Bx_beamline,
                                                     self.parameter)
        self.By, self.dBydz, self.Bydz = read_fields(self.By_beamline,
                                                     self.parameter)
        self.Bz, self.dBzdz, self.Bzdz = read_fields(self.Bz_beamline,
                                                     self.parameter)
        self.Ez, self.dEzdz, self.Ezdz = read_fields(self.Ez_beamline,
                                                     self.parameter)
        self.Gz, self.dGzdz, self.Gzdz = read_fields(self.Gz_beamline,
                                                     self.parameter)
        Dx_Bz, Dxp_Bz, Dy_Bz, Dyp_Bz = read_offsets(self.Bz_beamline,
                                                    self.parameter)
        Dx_Ez, Dxp_Ez, Dy_Ez, Dyp_Ez = read_offsets(self.Ez_beamline,
                                                    self.parameter)

        self.Dx = interpolate.interp1d(
            self.z,
            Dx_Bz(self.z) + Dx_Ez(self.z),
            kind='linear',
            fill_value=(0, 0),
            bounds_error=False
        )
        self.Dxp = interpolate.interp1d(
            self.z,
            Dxp_Bz(self.z) + Dxp_Ez(self.z),
            kind='linear',
            fill_value=(0, 0),
            bounds_error=False
        )
        self.Dy = interpolate.interp1d(
            self.z,
            Dy_Bz(self.z) + Dy_Ez(self.z),
            kind='linear',
            fill_value=(0, 0),
            bounds_error=False
        )
        self.Dyp = interpolate.interp1d(
            self.z,
            Dyp_Bz(self.z) + Dyp_Ez(self.z),
            kind='linear',
            fill_value=(0, 0),
            bounds_error=False
        )

    def __str__(self):
        string = 'Accelerator structure.\n'
        string += '\tSolenoids:\n'
        for element in self.Bz_beamline.values():
            string += (f"\t[ {element.z0} m, {element.max_field} T, "
                       f"{element.file_name}, {element.name}, "
                       f"{element.x} m, {element.xp} rad, "
                       f"{element.y} m, {element.yp} rad] \n")
        string += '\tAccelerating modules:\n'
        for element in self.Ez_beamline.values():
            string += (f"\t[ {element.z0} m, {element.max_field} T, "
                       f"{element.file_name}, {element.name}, "
                       f"{element.x} m, {element.xp} rad, "
                       f"{element.y} m, {element.yp} rad] \n")
        string += '\tQuadrupoles:\n'
        for element in self.Gz_beamline.values():
            string += (f"\t[ {element.z0} m, {element.max_field} T, "
                       f"{element.file_name}, {element.name} \n")
        string += '\tCorrectors x:\n'
        for element in self.By_beamline.values():
            string += (f"\t[ {element.z0} m, {element.max_field} T, "
                       f"{element.file_name}, {element.name} \n")
        string += '\tCorrectors y:\n'
        for element in self.Bx_beamline.values():
            string += (f"\t[ {element.z0} m, {element.max_field} T, "
                       f"{element.file_name}, {element.name} \n")
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
