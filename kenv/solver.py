# solver.py
'''Simulation of the envelope beam in the accelerator.
'''
import numpy as np

from .beam import *
from .accelerator import *
from scipy import interpolate, integrate
from scipy.integrate import solve_ivp

__all__ = ['Sim',
           'Simulation']

class Simulation:
    '''Simulation of the envelope beam in the accelerator.

    basic parameters after track:
    envelope_x and envelope_y,
    centroid_x and centroid_y,
    larmor_angle
    '''
    mass_rest_electron = 0.511
    clight = 299792458
    alfven_current = 17000

    def __init__(self,
             beam,
             accelerator):

        self.current = beam.current
        self.energy = beam.energy
        self.x = beam.x
        self.y = beam.y
        self.radius_x = beam.radius_x
        self.radius_y = beam.radius_y
        self.xp = beam.xp
        self.yp = beam.yp
        self.normalized_emittance_x = beam.normalized_emittance_x
        self.normalized_emittance_y = beam.normalized_emittance_y
        self.phase = beam.larmor_angle

        self.start = accelerator.start
        self.stop = accelerator.stop
        self.step = accelerator.step

        self.z = self.parameter = accelerator.parameter
        self.Bz = accelerator.Bz
        self.Ez = accelerator.Ez
        self.Gz = accelerator.Gz
        self.dBzdz = accelerator.dBzdz
        self.dEzdz = accelerator.dEzdz
        self.dGzdz = accelerator.dGzdz
        self.gamma_0 = self.energy/self.mass_rest_electron + 1
        self.gamma = self.gamma_0 + integrate.cumtrapz(-self.Ez(self.parameter)/self.mass_rest_electron, self.parameter)
        self.gamma = interpolate.interp1d(self.parameter[1:], self.gamma, fill_value=(self.gamma[0], self.gamma[-1]), bounds_error=False)

        self.envelope_x = []
        self.envelope_y = []

        self.centroid_x = []
        self.centroid_y = []
        self.larmor_angle = []

    def track(self, rtol:float = 1e-6):
        '''Tracking!
        '''

        def dXdz_beam(z:np.arange, X:list,
                 dz:float=self.step,
                 gamma:interpolate.interp1d=self.gamma,
                 Bz:interpolate.interp1d=self.Bz,
                 Ez:interpolate.interp1d=self.Ez,
                 Gz:interpolate.interp1d=self.Gz,
                 dBzdz:interpolate.interp1d=self.dBzdz,
                 dEzdz:interpolate.interp1d=self.dEzdz,
                 dGzdz:interpolate.interp1d=self.dGzdz,
                 beam_current:float=self.current,
                 emitt_n_x:float=self.normalized_emittance_x,
                 emitt_n_y:float=self.normalized_emittance_y,
                 c:float=self.clight,
                 mc:float=self.mass_rest_electron,
                 alfven_current:float=self.alfven_current) -> list:
            '''Located derivative for further integration.
            '''
            sigma_x = X[0]
            x = X[1]
            sigma_y = X[2]
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

            dsigma_xdz = x
            dxdz = 2 * P / (sigma_x + sigma_y) + emitt_x * emitt_x / (sigma_x * sigma_x * sigma_x) - K_x * sigma_x - \
                   dgdz * x / (beta * beta * g) - d2gdz2 * sigma_x/(2 * beta * beta * g)
            dsigma_ydz = y
            dydz = 2 * P/(sigma_x + sigma_y) + emitt_y * emitt_y / (sigma_y * sigma_y * sigma_y) - K_y * sigma_y - \
                   dgdz * y/(beta * beta * g) - d2gdz2 * sigma_y/(2 * beta * beta * g)


            return [dsigma_xdz, dxdz, dsigma_ydz, dydz]

        def dXdz_centroid(z:np.arange, X:list,
                 dz:float=self.step,
                 gamma:interpolate.interp1d=self.gamma,
                 Bz:interpolate.interp1d=self.Bz,
                 Ez:interpolate.interp1d=self.Ez,
                 Gz:interpolate.interp1d=self.Gz,
                 dBzdz:interpolate.interp1d=self.dBzdz,
                 dEzdz:interpolate.interp1d=self.dEzdz,
                 dGzdz:interpolate.interp1d=self.dGzdz,
                 beam_current:float=.0e0,
                 emitt_n_x:float=.0e0,
                 emitt_n_y:float=.0e0,
                 c:float=self.clight,
                 mc:float=self.mass_rest_electron,
                 alfven_current:float=self.alfven_current,
                 rtol:float=rtol) -> list:
            '''Located derivative for further integration.
            '''
            sigma_x = X[0]
            x = X[1]
            sigma_y = X[2]
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

            dsigma_xdz = x
            dxdz =  - K_x * sigma_x - dgdz * x / (beta * beta * g) - d2gdz2 * sigma_x/(2 * beta * beta * g)
            dsigma_ydz = y
            dydz = - K_y * sigma_y - dgdz * y/(beta * beta * g) - d2gdz2 * sigma_y/(2 * beta * beta * g)
            dphidz = -K_s**0.5

            return [dsigma_xdz, dxdz, dsigma_ydz, dydz, dphidz]


        X0_beam = np.array([self.radius_x, self.xp, self.radius_y, self.yp])
        X0_centroid = np.array([self.x, 0, self.y, 0, self.phase])

        beam = solve_ivp(dXdz_beam, t_span=[self.parameter[0],self.parameter[-1]], y0=X0_beam, t_eval=self.parameter, rtol=rtol).y
        centroid = solve_ivp(dXdz_centroid, t_span=[self.parameter[0],self.parameter[-1]], y0=X0_centroid, t_eval=self.parameter, rtol=rtol).y

        self.envelope_x = beam[0,:]
        self.envelope_y = beam[2,:]

        self.centroid_x = centroid[0,:]*np.cos(centroid[4,:]) + centroid[2,:]*np.sin(centroid[4,:])
        self.centroid_y = -centroid[0,:]*np.sin(centroid[4,:]) + centroid[2,:]*np.cos(centroid[4,:])
        self.larmor_angle = centroid[4,:]

Sim = Simulation
