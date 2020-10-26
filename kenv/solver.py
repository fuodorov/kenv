# solver.py
'''Simulation of the envelope beam in the accelerator.'''

import numpy as np
from scipy import interpolate, integrate, misc
from scipy.integrate import solve_ivp
from .constants import *
from .beam import *

__all__ = ['Sim',
           'Simulation',
           'Equations']

class Equations:
    '''Located derivative for further integration.'''

    def __init__(self, beam, accelerator,
                 particle=Particle()):
        self.beam = beam
        self.accelerator = accelerator
        self.particle = particle

    def envelope_prime(self,
                       z:np.arange,
                       X:list) -> list:
        '''Located derivative for further integration
         Kapchinscky equation for envelope beam.

        '''

        x = X[0]
        xp = X[1]
        y = X[2]
        yp = X[3]

        g = self.beam.gamma + self.beam.charge*self.accelerator.Ezdz(z)/mass_rest_electron
        dgdz = self.beam.charge*self.accelerator.Ez(z)/mass_rest_electron
        d2gdz2 = self.beam.charge*self.accelerator.dEzdz(z)/mass_rest_electron
        beta = np.sqrt(1 - 1 / (g*g))
        p = g*beta*mass_rest_electron

        K_s = (self.beam.charge*speed_light*self.accelerator.Bz(z) / (2*p*MeV))**2
        K_q = (self.beam.charge*speed_light*self.accelerator.Gz(z) / (p*MeV))
        K_x = K_s - K_q
        K_y = K_s + K_q

        emitt_x = self.beam.normalized_emittance_x / (g*beta)
        emitt_y = self.beam.normalized_emittance_y / (g*beta)

        P = 2*self.beam.current / (alfven_current * (g*beta)**3)

        dxdz = xp
        dxpdz = 2*P / (x + y) + emitt_x*emitt_x / x**3 - K_x*x - \
                dgdz*xp / (beta*beta*g) - d2gdz2*x / (2*beta*beta*g)
        dydz = yp
        dypdz = 2*P / (x + y) + emitt_y*emitt_y / y**3 - K_y*y - \
                dgdz*yp / (beta*beta*g) - d2gdz2*y / (2*beta*beta*g)

        return [dxdz, dxpdz, dydz, dypdz]

    def centroid_prime(self,
                       z:np.arange,
                       X:list) -> list:
        '''Located derivative for further integration
         Kapchinscky equation for centroid trajectory.

        '''

        x = X[0]
        xp = X[1]
        y = X[2]
        yp = X[3]
        phi = X[4]

        g = self.beam.gamma + self.beam.charge*self.accelerator.Ezdz(z)/mass_rest_electron
        beta = np.sqrt(1 - 1 / (g*g))
        p = g*beta*mass_rest_electron*MeV/speed_light

        offset_x = self.accelerator.Dx(z)
        offset_xp = self.accelerator.Dxp(z)
        offset_y = self.accelerator.Dy(z)
        offset_yp = self.accelerator.Dyp(z)
        x_corr = x - offset_x
        y_corr = y - offset_y
        r_corr = np.sqrt((x_corr)**2 + (y_corr)**2)

        Bx = self.accelerator.Bx(z)
        By = self.accelerator.By(z)
        Bz = self.accelerator.Bz(z)
        dBzdz = self.accelerator.dBzdz(z)
        d2Bzdz2 = misc.derivative(self.accelerator.dBzdz, z, dx=self.accelerator.dz, n=1)
        Gz = self.accelerator.Gz(z)
        Bz = Bz - d2Bzdz2*r_corr**2/4                          # row remainder
        Bx = Bx + Gz*y_corr - dBzdz*x_corr/2 + Bz*offset_xp    # row remainder
        By = By + Gz*x_corr - dBzdz*y_corr/2 + Bz*offset_yp    # row remainder
        Brho = p/self.beam.charge

        Ez = self.accelerator.Ez(z)*MeV
        dEzdz = self.accelerator.dEzdz(z)*MeV
        d2Ezdz2 = misc.derivative(self.accelerator.dEzdz, z, dx=self.accelerator.dz, n=1)*MeV
        Ez = Ez - d2Ezdz2*r_corr**2/4                    # row remainder
        Ex = - dEzdz*x_corr/2 + Ez*offset_xp             # row remainder
        Ey = - dEzdz*y_corr/2 + Ez*offset_yp             # row remainder

        dxdz = xp
        dxpdz = (Ex - Ez*xp) / (beta*speed_light*Brho) - (By - yp*Bz) / Brho
        dydz = yp
        dypdz = (Ey - Ez*yp) / (beta*speed_light*Brho) + (Bx - xp*Bz) / Brho
        dphidz = Bz / (2*Brho)

        return [dxdz, dxpdz, dydz, dypdz, dphidz]

    def particle_prime(self,
                       z:np.arange,
                       X:list, envelope_x, envelope_y) -> list:
        '''Located derivative for further integration
         Kapchinscky equation for envelope beam.

        '''

        x = X[0]
        xp = X[1]
        y = X[2]
        yp = X[3]
        phi = X[4]

        if self.particle.energy == 0.0:
            g = self.beam.gamma + self.beam.charge*self.accelerator.Ezdz(z)/mass_rest_electron
        else:
            g = self.particle.gamma + self.beam.charge*self.accelerator.Ezdz(z)/mass_rest_electron
        dgdz = self.beam.charge*self.accelerator.Ez(z)/mass_rest_electron
        d2gdz2 = self.beam.charge*self.accelerator.dEzdz(z)/mass_rest_electron
        beta = np.sqrt(1 - 1 / (g*g))
        p = g*beta*mass_rest_electron

        K_s = (self.beam.charge*speed_light*self.accelerator.Bz(z) / (2*p*MeV))**2
        K_q = (self.beam.charge*speed_light*self.accelerator.Gz(z) / (p*MeV))
        K_x = K_s - K_q
        K_y = K_s + K_q

        Brho = g*beta*mass_rest_electron*MeV/speed_light/self.beam.charge

        P = 2*self.beam.current / (alfven_current * (g*beta)**3)
        if (x*x + y*y > envelope_x(z)**2 + envelope_y(z)**2) and (x != 0.0  and y !=0.0):
            dxdz = xp
            dxpdz = 2*P * ((envelope_x(z)+envelope_y(z))/2*envelope_x(z)) / x - K_x*x - \
                    dgdz*xp / (beta*beta*g) - d2gdz2*x / (2*beta*beta*g)
            dydz = yp
            dypdz = 2*P * ((envelope_x(z)+envelope_y(z))/2*envelope_y(z)) / y - K_y*y - \
                    dgdz*yp / (beta*beta*g) - d2gdz2*y / (2*beta*beta*g)
        else:
            dxdz = xp
            dxpdz = 2*P / ((envelope_x(z)+envelope_y(z))/2*envelope_x(z)) * x - K_x*x - \
                    dgdz*xp / (beta*beta*g) - d2gdz2*x / (2*beta*beta*g)
            dydz = yp
            dypdz = 2*P / ((envelope_x(z)+envelope_y(z))/2*envelope_y(z)) * y - K_y*y - \
                    dgdz*yp / (beta*beta*g) - d2gdz2*y / (2*beta*beta*g)

        dphidz = self.accelerator.Bz(z) / (2*Brho)

        return [dxdz, dxpdz, dydz, dypdz, dphidz]

class Simulation:
    '''Simulation of the envelope beam in the accelerator.

    Basic parameters after track:
    gamma,
    envelope_x and envelope_y,
    envelope_xp and envelope_yp,
    centroid_x and centroid_y,
    centroid_xp and centroid_yp,
    larmor_angle

    '''

    def __init__(self,
                 beam,
                 accelerator,
                 particle=Particle()):

        self.beam = beam
        self.accelerator = accelerator
        self.particle = particle
        self.envelope_x = []
        self.envelope_xp = []
        self.envelope_y = []
        self.envelope_yp = []

        self.centroid_x = []
        self.centroid_xp = []
        self.centroid_y = []
        self.centroid_yp = []

        self.particle_x = []
        self.particle_xp = []
        self.particle_y = []
        self.particle_yp = []

        self.larmor_angle = []
        self.particle_larmor_angle = []

    def track(self,
              particle: bool=False,
              rtol: float=1e-5,
              atol: float=1e-8,
              method: str='RK23'):
        '''Tracking!'''

        # initial conditions
        equations = Equations(self.beam, self.accelerator, self.particle)

        X0_beam = np.array([self.beam.radius_x, self.beam.radius_xp,
                            self.beam.radius_y, self.beam.radius_yp])
        X0_centroid = np.array([self.beam.x, self.beam.xp,
                                self.beam.y, self.beam.yp,
                                self.beam.larmor_angle])
        X0_particle = np.array([self.particle.x, self.particle.xp,
                                self.particle.y, self.particle.yp,
                                self.particle.larmor_angle])
        # solver
        beam_envelope = solve_ivp(equations.envelope_prime,
        t_span=[self.accelerator.parameter[0], self.accelerator.parameter[-1]],
        y0=X0_beam, t_eval=self.accelerator.parameter, method=method, rtol=rtol, atol=atol).y

        self.gamma = self.beam.gamma + self.beam.charge*self.accelerator.Ezdz(self.accelerator.parameter)/mass_rest_electron

        self.envelope_x = beam_envelope[0,:]
        self.envelope_xp = beam_envelope[1,:]
        self.envelope_y = beam_envelope[2,:]
        self.envelope_yp = beam_envelope[3,:]

        self.envelope_x = interpolate.interp1d(self.accelerator.parameter, self.envelope_x, kind='cubic', fill_value=(0, 0), bounds_error=False)
        self.envelope_xp = interpolate.interp1d(self.accelerator.parameter, self.envelope_xp, kind='cubic', fill_value=(0, 0), bounds_error=False)
        self.envelope_y = interpolate.interp1d(self.accelerator.parameter, self.envelope_y, kind='cubic', fill_value=(0, 0), bounds_error=False)
        self.envelope_yp = interpolate.interp1d(self.accelerator.parameter, self.envelope_yp, kind='cubic', fill_value=(0, 0), bounds_error=False)

        centroid_trajectory = solve_ivp(equations.centroid_prime,
        t_span=[self.accelerator.parameter[0], self.accelerator.parameter[-1]],
        y0=X0_centroid, t_eval=self.accelerator.parameter, method=method, rtol=rtol, atol=atol).y

        self.centroid_x = centroid_trajectory[0,:]
        self.centroid_xp = centroid_trajectory[1,:]
        self.centroid_y = centroid_trajectory[2,:]
        self.centroid_yp = centroid_trajectory[3,:]
        phi = self.larmor_angle = centroid_trajectory[4,:]

        self.centroid_x = interpolate.interp1d(self.accelerator.parameter, self.centroid_x, kind='cubic', fill_value=(0, 0), bounds_error=False)
        self.centroid_xp = interpolate.interp1d(self.accelerator.parameter, self.centroid_xp, kind='cubic', fill_value=(0, 0), bounds_error=False)
        self.centroid_y = interpolate.interp1d(self.accelerator.parameter, self.centroid_y, kind='cubic', fill_value=(0, 0), bounds_error=False)
        self.centroid_yp = interpolate.interp1d(self.accelerator.parameter, self.centroid_yp, kind='cubic', fill_value=(0, 0), bounds_error=False)
        self.larmor_angle = interpolate.interp1d(self.accelerator.parameter, self.larmor_angle, kind='cubic', fill_value=(0, 0), bounds_error=False)

        self.gamma = interpolate.interp1d(self.accelerator.parameter, self.gamma, kind='cubic', fill_value=(0, 0), bounds_error=False)

        if particle == True:
            def wrapper(t, y):
                return equations.particle_prime(t,y, self.envelope_x, self.envelope_y)

            particle_trajectory = solve_ivp(wrapper,
            t_span=[self.accelerator.parameter[0], self.accelerator.parameter[-1]],
            y0=X0_particle, t_eval=self.accelerator.parameter, method=method, rtol=rtol, atol=atol).y

            psi = self.particle_larmor_angle = particle_trajectory[4,:]
            self.particle_x = particle_trajectory[0,:]*np.cos(psi) - particle_trajectory[2,:]*np.sin(psi)
            self.particle_y = particle_trajectory[0,:]*np.sin(psi) + particle_trajectory[2,:]*np.cos(psi)
            self.particle_xp = particle_trajectory[1,:]*np.cos(psi) - particle_trajectory[3,:]*np.sin(psi)
            self.particle_yp = particle_trajectory[1,:]*np.sin(psi) + particle_trajectory[3,:]*np.cos(psi)

            self.particle_x = interpolate.interp1d(self.accelerator.parameter, self.particle_x, kind='cubic', fill_value=(0, 0), bounds_error=False)
            self.particle_xp = interpolate.interp1d(self.accelerator.parameter, self.particle_xp, kind='cubic', fill_value=(0, 0), bounds_error=False)
            self.particle_y = interpolate.interp1d(self.accelerator.parameter, self.particle_y, kind='cubic', fill_value=(0, 0), bounds_error=False)
            self.particle_yp = interpolate.interp1d(self.accelerator.parameter, self.particle_yp, kind='cubic', fill_value=(0, 0), bounds_error=False)
            self.particle_larmor_angle = interpolate.interp1d(self.accelerator.parameter, self.particle_larmor_angle, kind='cubic', fill_value=(0, 0), bounds_error=False)

Sim = Simulation
