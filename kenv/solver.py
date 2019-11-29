# solver.py
'''Simulation of the envelope beam in the accelerator.
'''
import numpy as np
from scipy import interpolate, integrate
from scipy.integrate import solve_ivp

from .beam import *
from .accelerator import *

__all__ = ['Sim',
           'Simulation']

class Simulation:
    '''Simulation of the envelope beam in the accelerator.

    basic parameters after track:
    gamma,
    envelope_x and envelope_y,
    centroid_x and centroid_y,
    larmor_angle
    '''

    mass_rest_electron = 0.5109989461
    clight = 299792458
    alfven_current = 17000

    gamma = 0

    envelope_x = []
    envelope_y = []

    centroid_x = []
    centroid_y = []
    larmor_angle = []

    def __init__(self,
             beam,
             accelerator):

        self.beam = beam
        self.accelerator = accelerator

    def track(self, rtol:float = 1e-6):
        '''Tracking!
        '''

        self.gamma = (self.beam.energy/self.mass_rest_electron + 1) + integrate.cumtrapz(-self.accelerator.Ez(self.accelerator.parameter)/self.mass_rest_electron, self.accelerator.parameter)
        self.gamma = interpolate.interp1d(self.accelerator.parameter[1:], self.gamma, fill_value=(self.gamma[0], self.gamma[-1]), bounds_error=False)

        def dXdz_beam(z:np.arange, X:list,
                 dz:float=self.accelerator.step,
                 gamma:interpolate.interp1d=self.gamma,
                 Bz:interpolate.interp1d=self.accelerator.Bz,
                 Ez:interpolate.interp1d=self.accelerator.Ez,
                 Gz:interpolate.interp1d=self.accelerator.Gz,
                 dBzdz:interpolate.interp1d=self.accelerator.dBzdz,
                 dEzdz:interpolate.interp1d=self.accelerator.dEzdz,
                 dGzdz:interpolate.interp1d=self.accelerator.dGzdz,
                 beam_current:float=self.beam.current,
                 emitt_n_x:float=self.beam.normalized_emittance_x,
                 emitt_n_y:float=self.beam.normalized_emittance_y,
                 c:float=self.clight,
                 mc:float=self.mass_rest_electron,
                 alfven_current:float=self.alfven_current) -> list:
            '''Located derivative for further integration.
            '''
            radius_x = X[0]
            x = X[1]
            radius_y = X[2]
            y = X[3]

            g = gamma(z)
            dgdz = Ez(z) / mc

            d2gdz2 = dEzdz(z) / mc

            beta = np.sqrt(1 - 1 / (g * g))

            K_s = ( c * Bz(z) / (2*mc*1e6) / (g * beta))**2
            K_q = ( c * Gz(z) / (mc*1e6)) / (g * beta)
            K_x = K_s + K_q
            K_y = K_s - K_q

            emitt_x = emitt_n_x / (g * beta)
            emitt_y = emitt_n_y / (g * beta)

            P = 2 * beam_current / (alfven_current * beta * beta * beta * g * g * g)

            dradius_xdz = x
            dxdz = 2 * P / (radius_x + radius_y) + emitt_x * emitt_x / (radius_x * radius_x * radius_x) - K_x * radius_x - \
                   dgdz * x / (beta * beta * g) - d2gdz2 * radius_x/(2 * beta * beta * g)
            dradius_ydz = y
            dydz = 2 * P/(radius_x + radius_y) + emitt_y * emitt_y / (radius_y * radius_y * radius_y) - K_y * radius_y - \
                   dgdz * y/(beta * beta * g) - d2gdz2 * radius_y/(2 * beta * beta * g)

            return [dradius_xdz, dxdz, dradius_ydz, dydz]

        def dXdz_centroid(z:np.arange, X:list,
                 dz:float=self.accelerator.step,
                 gamma:interpolate.interp1d=self.gamma,
                 Bz:interpolate.interp1d=self.accelerator.Bz,
                 Ez:interpolate.interp1d=self.accelerator.Ez,
                 Gz:interpolate.interp1d=self.accelerator.Gz,
                 dBzdz:interpolate.interp1d=self.accelerator.dBzdz,
                 dEzdz:interpolate.interp1d=self.accelerator.dEzdz,
                 dGzdz:interpolate.interp1d=self.accelerator.dGzdz,
                 c:float=self.clight,
                 mc:float=self.mass_rest_electron,
                 alfven_current:float=self.alfven_current,
                 rtol:float=rtol) -> list:
            '''Located derivative for further integration.
            '''
            radius_x = X[0]
            x = X[1]
            radius_y = X[2]
            y = X[3]
            phi = X[4]

            g = gamma(z)
            dgdz = Ez(z) / mc
            d2gdz2 = dEzdz(z) / mc

            beta = np.sqrt(1 - 1 / (g * g))

            K_s = ( c * Bz(z) / (2*mc*1e6) / (g * beta))**2
            K_q = ( c * Gz(z) / (mc*1e6)) / (g * beta)
            K_x = K_s + K_q
            K_y = K_s - K_q

            dradius_xdz = x
            dxdz =  - K_x * radius_x - dgdz * x / (beta * beta * g) - d2gdz2 * radius_x/(2 * beta * beta * g)
            dradius_ydz = y
            dydz = - K_y * radius_y - dgdz * y/(beta * beta * g) - d2gdz2 * radius_y/(2 * beta * beta * g)
            dphidz = -K_s**0.5

            return [dradius_xdz, dxdz, dradius_ydz, dydz, dphidz]

        X0_beam = np.array([self.beam.radius_x, self.beam.radius_xp, self.beam.radius_y, self.beam.radius_yp])
        X0_centroid = np.array([self.beam.x, self.beam.xp, self.beam.y, self.beam.yp, self.beam.larmor_angle])

        beam = solve_ivp(dXdz_beam, t_span=[self.accelerator.parameter[0],self.accelerator.parameter[-1]], y0=X0_beam, t_eval=self.accelerator.parameter, rtol=rtol).y
        centroid = solve_ivp(dXdz_centroid, t_span=[self.accelerator.parameter[0],self.accelerator.parameter[-1]], y0=X0_centroid, t_eval=self.accelerator.parameter, rtol=rtol).y

        self.envelope_x = beam[0,:]
        self.envelope_y = beam[2,:]

        self.centroid_x = centroid[0,:]*np.cos(centroid[4,:]) + centroid[2,:]*np.sin(centroid[4,:])
        self.centroid_y = -centroid[0,:]*np.sin(centroid[4,:]) + centroid[2,:]*np.cos(centroid[4,:])
        self.larmor_angle = centroid[4,:]

        #TO function
        self.envelope_x = interpolate.interp1d(self.accelerator.parameter, self.envelope_x, kind='cubic', fill_value=(0, 0), bounds_error=False)
        self.envelope_y = interpolate.interp1d(self.accelerator.parameter, self.envelope_y, kind='cubic', fill_value=(0, 0), bounds_error=False)

        self.centroid_x = interpolate.interp1d(self.accelerator.parameter, self.centroid_x, kind='cubic', fill_value=(0, 0), bounds_error=False)
        self.centroid_y = interpolate.interp1d(self.accelerator.parameter, self.centroid_y, kind='cubic', fill_value=(0, 0), bounds_error=False)
        self.larmor_angle = interpolate.interp1d(self.accelerator.parameter, self.larmor_angle, kind='cubic', fill_value=(0, 0), bounds_error=False)




Sim = Simulation
