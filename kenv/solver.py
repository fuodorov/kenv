# solver.py
'''Simulation of the envelope beam in the accelerator.
'''
import numpy as np

from .beam import *
from .accelerator import *
from scipy import interpolate, integrate
from scipy.integrate import odeint, solve_ivp

__all__ = ['Sim',
           'Simulation']

class Simulation:
    '''Simulation of the envelope beam in the accelerator.
    '''
    mass_rest_electron = 0.511
    clight = 299792458
    alfven_current = 17000

    def __init__(self,
             beam,
             accelerator):

        self.current = beam.current
        self.energy = beam.energy
        self.radius_x = beam.radius_x
        self.radius_y = beam.radius_y
        self.angular_x = beam.angular_x
        self.angular_y = beam.angular_y
        self.normalized_emittans_x = beam.normalized_emittans_x
        self.normalized_emittans_y = beam.normalized_emittans_y

        self.particles = beam.Particles

        self.start = accelerator.start
        self.stop = accelerator.stop
        self.step = accelerator.step
        self.parameter = accelerator.parameter
        self.Bz = accelerator.Bz
        self.Ez = accelerator.Ez
        self.Gz = accelerator.Gz

        self.gamma_0 = self.energy/self.mass_rest_electron + 1
        self.gamma = self.gamma_0 + integrate.cumtrapz(-self.Ez(self.parameter)/self.mass_rest_electron, self.parameter)
        self.gamma = interpolate.interp1d(self.parameter[1:], self.gamma, fill_value=(self.gamma[0], self.gamma[-1]), bounds_error=False)

        self.particle_x = {}
        self.particle_y = {}

        self.envelope_x = []
        self.envelope_y = []

    def track(self, solver = 'odeint', rtol:float = 1e-5):
        '''Tracking!
        '''

        if solver == 'odeint':
            def dXdz(X:list, z:np.arange,
                     dz:float=self.step,
                     gamma:interpolate.interp1d=self.gamma,
                     Bz:interpolate.interp1d=self.Bz,
                     Ez:interpolate.interp1d=self.Ez,
                     Gz:interpolate.interp1d=self.Gz,
                     beam_current:float=self.current,
                     emitt_n_x:float=self.normalized_emittans_x,
                     emitt_n_y:float=self.normalized_emittans_y,
                     c:float=self.clight,
                     mc:float=self.mass_rest_electron,
                     alfven_current:float=self.alfven_current) ->list:
                '''Located derivative for further integration.
                '''
                sigma_x = X[0]
                x = X[1]
                sigma_y = X[2]
                y = X[3]

                g = gamma(z)
                dgdz = Ez(z) / mc
                d2gdz2 = (Ez(z+dz) - Ez(z) / dz / mc)

                beta = np.sqrt(1 - 1 / (g * g))

                K_s = ( c * Bz(z) / (2*mc*1e6) / (g * beta))**2
                K_q = ( c * Gz(z) / (mc*1e6)) / (g * beta)
                K_x = K_s + K_q
                K_y = K_s - K_q

                emitt_x = emitt_n_x / (g * beta)
                emitt_y = emitt_n_y / (g * beta)

                P = 2 * beam_current / (alfven_current * beta * beta * beta * g * g * g)

                dxdz = 2 * P / (sigma_x + sigma_y) + emitt_x * emitt_x / (sigma_x * sigma_x * sigma_x) - K_x * sigma_x - \
                       dgdz * x / (beta * beta * g) - d2gdz2 * sigma_x/(2 * beta * beta * g)
                dsigma_xdz = x

                dydz = 2 * P/(sigma_x + sigma_y) + emitt_y * emitt_y / (sigma_y * sigma_y * sigma_y) - K_y * sigma_y - \
                       dgdz * y/(beta * beta * g) - d2gdz2 * sigma_y/(2 * beta * beta * g)
                dsigma_ydz = y

                return [dsigma_xdz, dxdz, dsigma_ydz, dydz]

            X0 = [self.radius_x, self.angular_x, self.radius_y, self.angular_y]
            self.envelope_x = odeint(dXdz, X0, self.parameter, rtol=rtol)[0:,0]
            self.envelope_y = odeint(dXdz, X0, self.parameter, rtol=rtol)[0:,2]

        elif solver == 'solve_ivp':
            def dXdz(z:np.arange, X:list,
                     dz:float=self.step,
                     gamma:interpolate.interp1d=self.gamma,
                     Bz:interpolate.interp1d=self.Bz,
                     Ez:interpolate.interp1d=self.Ez,
                     Gz:interpolate.interp1d=self.Gz,
                     beam_current:float=self.current,
                     emitt_n_x:float=self.normalized_emittans_x,
                     emitt_n_y:float=self.normalized_emittans_y,
                     c:float=self.clight,
                     mc:float=self.mass_rest_electron,
                     alfven_current:float=self.alfven_current) ->list:
                '''Located derivative for further integration.
                '''
                sigma_x = X[0]
                x = X[1]
                sigma_y = X[2]
                y = X[3]

                g = gamma(z)
                dgdz = Ez(z) / mc
                d2gdz2 = (Ez(z+dz) - Ez(z) / dz / mc)

                beta = np.sqrt(1 - 1 / (g * g))

                K_s = ( c * Bz(z) / (2*mc*1e6) / (g * beta))**2
                K_q = ( c * Gz(z) / (mc*1e6)) / (g * beta)
                K_x = K_s + K_q
                K_y = K_s - K_q

                emitt_x = emitt_n_x / (g * beta)
                emitt_y = emitt_n_y / (g * beta)

                P = 2 * beam_current / (alfven_current * beta * beta * beta * g * g * g)

                dxdz = 2 * P / (sigma_x + sigma_y) + emitt_x * emitt_x / (sigma_x * sigma_x * sigma_x) - K_x * sigma_x - \
                       dgdz * x / (beta * beta * g) - d2gdz2 * sigma_x/(2 * beta * beta * g)
                dsigma_xdz = x

                dydz = 2 * P/(sigma_x + sigma_y) + emitt_y * emitt_y / (sigma_y * sigma_y * sigma_y) - K_y * sigma_y - \
                       dgdz * y/(beta * beta * g) - d2gdz2 * sigma_y/(2 * beta * beta * g)
                dsigma_ydz = y

                return [dsigma_xdz, dxdz, dsigma_ydz, dydz]

            X0 = np.array([self.radius_x, self.angular_x, self.radius_y, self.angular_y])
            self.envelope_x = solve_ivp(dXdz, t_span=[self.parameter[0],self.parameter[-1]], y0=X0, t_eval=self.parameter, rtol=rtol).y[0,:]
            self.envelope_y = solve_ivp(dXdz, t_span=[self.parameter[0],self.parameter[-1]], y0=X0, t_eval=self.parameter, rtol=rtol).y[2,:]


Sim = Simulation
