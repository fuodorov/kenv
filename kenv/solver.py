# solver.py
'''Simulation of the envelope beam in the accelerator.'''

import numpy as np
from scipy import interpolate, integrate
from scipy.integrate import solve_ivp
from .constants import *

__all__ = ['Sim',
           'Simulation',
           'KapchinskyEquations']

class KapchinskyEquations:
    '''Located derivative for further integration.'''

    def __init__(self, beam, accelerator):
        self.beam = beam
        self.accelerator = accelerator

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
        phi = X[4]

        g = self.beam.gamma + self.beam.charge*self.accelerator.Ezdz(z)/mass_rest_electron
        dgdz = self.beam.charge*self.accelerator.Ez(z)/mass_rest_electron
        d2gdz2 = self.beam.charge*self.accelerator.dEzdz(z)/mass_rest_electron
        beta = np.sqrt(1 - 1 / (g*g))
        p = g*beta*mass_rest_electron

        K_s = (speed_light*self.accelerator.Bz(z) / (2*p*1e6))**2
        K_q = (speed_light*self.accelerator.Gz(z) / (p*1e6))
        K_corr_x = (speed_light*self.accelerator.By(z) / (p*1e6))
        K_corr_y = (speed_light*self.accelerator.Bx(z) / (p*1e6))
        K_x = K_s + K_q
        K_y = K_s - K_q

        emitt_x = self.beam.normalized_emittance_x / (p/mass_rest_electron)
        emitt_y = self.beam.normalized_emittance_y / (p/mass_rest_electron)

        P = 2*self.beam.current / (alfven_current*(p/mass_rest_electron)**3)

        dxdz = xp
        dxpdz = 2*P / (x + y) + emitt_x*emitt_x / x**3 - K_x*x - \
                dgdz*xp / (beta*beta*g) - d2gdz2*x / (2*beta*beta*g) - K_corr_x
        dydz = yp
        dypdz = 2*P / (x + y) + emitt_y*emitt_y / y**3 - K_y*y - \
                dgdz*yp / (beta*beta*g) - d2gdz2*y / (2*beta*beta*g) + K_corr_y
        dphidz = -K_s**0.5

        return [dxdz, dxpdz, dydz, dypdz, dphidz]


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

        dgdz = self.beam.charge*self.accelerator.Ez(z)/mass_rest_electron
        d2gdz2 = self.beam.charge*self.accelerator.dEzdz(z)/mass_rest_electron

        g = self.beam.gamma + self.beam.charge*self.accelerator.Ezdz(z)/mass_rest_electron
        beta = np.sqrt(1 - 1 / (g*g))
        p = g*beta*mass_rest_electron

        K_s = (speed_light*self.accelerator.Bz(z) / (2*p*1e6))**2
        K_q = (speed_light*self.accelerator.Gz(z) / (p*1e6))
        K_corr_x = (speed_light*self.accelerator.By(z) / (p*1e6))
        K_corr_y = (speed_light*self.accelerator.Bx(z) / (p*1e6))
        K_x = K_s + K_q
        K_y = K_s - K_q

        dxdz = xp
        dxpdz = - K_x*x - dgdz*xp / (beta*beta*g) - d2gdz2*x / (2*beta*beta*g) - K_corr_x
        dydz = yp
        dypdz = - K_y*y - dgdz*yp / (beta*beta*g) - d2gdz2*y / (2*beta*beta*g) + K_corr_y
        dphidz = -K_s**0.5

        return [dxdz, dxpdz, dydz, dypdz, dphidz]

class Simulation:
    '''Simulation of the envelope beam in the accelerator.

    Basic parameters after track:
    gamma,
    envelope_x and envelope_y,
    centroid_x and centroid_y,
    larmor_angle

    '''

    envelope_x = []
    envelope_xp = []
    envelope_y = []
    envelope_yp = []

    centroid_x = []
    centroid_xp = []
    centroid_y = []
    centroid_yp = []

    larmor_angle = []

    def __init__(self,
                 beam,
                 accelerator):

        self.beam = beam
        self.accelerator = accelerator

    def track(self,
              rtol:float = 1e-6):
        '''Tracking!'''

        Equations = KapchinskyEquations(self.beam, self.accelerator)

        X0_beam = np.array([self.beam.radius_x, self.beam.radius_xp,
                            self.beam.radius_y, self.beam.radius_yp,
                            self.beam.larmor_angle])
        X0_centroid = np.array([self.beam.x, self.beam.xp,
                                self.beam.y, self.beam.yp,
                                self.beam.larmor_angle])

        beam = solve_ivp(Equations.envelope_prime, t_span=[self.accelerator.parameter[0],self.accelerator.parameter[-1]], y0=X0_beam, t_eval=self.accelerator.parameter, rtol=rtol).y
        centroid = solve_ivp(Equations.centroid_prime, t_span=[self.accelerator.parameter[0],self.accelerator.parameter[-1]], y0=X0_centroid, t_eval=self.accelerator.parameter, rtol=rtol).y

        self.envelope_x = beam[0,:]
        self.envelope_xp = beam[1,:]
        self.envelope_y = beam[2,:]
        self.envelope_yp = beam[3,:]

        self.centroid_x = centroid[0,:]*np.cos(centroid[4,:]) + centroid[2,:]*np.sin(centroid[4,:])
        self.centroid_y = -centroid[0,:]*np.sin(centroid[4,:]) + centroid[2,:]*np.cos(centroid[4,:])
        self.centroid_xp = centroid[1,:]*np.cos(centroid[4,:]) + centroid[3,:]*np.sin(centroid[4,:])
        self.centroid_yp = -centroid[1,:]*np.sin(centroid[4,:]) + centroid[3,:]*np.cos(centroid[4,:])
        self.larmor_angle = centroid[4,:]

        self.gamma = self.beam.gamma + self.beam.charge*self.accelerator.Ezdz(self.accelerator.parameter)/mass_rest_electron

        #TO function
        self.envelope_x = interpolate.interp1d(self.accelerator.parameter, self.envelope_x, kind='cubic', fill_value=(0, 0), bounds_error=False)
        self.envelope_xp = interpolate.interp1d(self.accelerator.parameter, self.envelope_xp, kind='cubic', fill_value=(0, 0), bounds_error=False)
        self.envelope_y = interpolate.interp1d(self.accelerator.parameter, self.envelope_y, kind='cubic', fill_value=(0, 0), bounds_error=False)
        self.envelope_yp = interpolate.interp1d(self.accelerator.parameter, self.envelope_yp, kind='cubic', fill_value=(0, 0), bounds_error=False)

        self.centroid_x = interpolate.interp1d(self.accelerator.parameter, self.centroid_x, kind='cubic', fill_value=(0, 0), bounds_error=False)
        self.centroid_xp = interpolate.interp1d(self.accelerator.parameter, self.centroid_xp, kind='cubic', fill_value=(0, 0), bounds_error=False)
        self.centroid_y = interpolate.interp1d(self.accelerator.parameter, self.centroid_y, kind='cubic', fill_value=(0, 0), bounds_error=False)
        self.centroid_yp = interpolate.interp1d(self.accelerator.parameter, self.centroid_yp, kind='cubic', fill_value=(0, 0), bounds_error=False)
        self.larmor_angle = interpolate.interp1d(self.accelerator.parameter, self.larmor_angle, kind='cubic', fill_value=(0, 0), bounds_error=False)

        self.gamma = interpolate.interp1d(self.accelerator.parameter, self.gamma, kind='cubic', fill_value=(0, 0), bounds_error=False)

Sim = Simulation
