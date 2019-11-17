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

        self.particles = beam.particles

        self.start = accelerator.start
        self.stop = accelerator.stop
        self.step = accelerator.step
        self.parameter = accelerator.parameter
        self.Bz = accelerator.Bz
        self.Ez = accelerator.Ez
        self.Gz = accelerator.Gz
        self.dBzdz = accelerator.dBzdz
        self.dEzdz = accelerator.dEzdz
        self.dGzdz = accelerator.dGzdz

        self.gamma_0 = self.energy/self.mass_rest_electron + 1
        self.gamma = self.gamma_0 + integrate.cumtrapz(-self.Ez(self.parameter)/self.mass_rest_electron, self.parameter)
        self.gamma = interpolate.interp1d(self.parameter[1:], self.gamma, fill_value=(self.gamma[0], self.gamma[-1]), bounds_error=False)

        self.particle_x = {}
        self.particle_angular_x = {}
        self.particle_y = {}
        self.particle_angular_y = {}
        self.particle_radius = {}

        self.particle_larmor_phase = {}
        self.particle_larmor_x = {}
        self.particle_larmor_angular_x = {}
        self.particle_larmor_y = {}
        self.particle_larmor_angular_y = {}

        self.envelope_x = self.beam_x = []
        self.beam_angular_x = []
        self.envelope_y = self.beam_y = []
        self.beam_angular_y = []

    def track(self, solver = 'odeint', rtol:float = 1e-8):
        '''Tracking!
        '''

        X0 = [self.radius_x, self.angular_x, self.radius_y, self.angular_y]

        if solver == 'odeint':

            def dXdz(X:list, z:np.arange,
                     dz:float=self.step,
                     gamma:interpolate.interp1d=self.gamma,
                     Bz:interpolate.interp1d=self.Bz,
                     Ez:interpolate.interp1d=self.Ez,
                     Gz:interpolate.interp1d=self.Gz,
                     dBzdz:interpolate.interp1d=self.dBzdz,
                     dEzdz:interpolate.interp1d=self.dEzdz,
                     dGzdz:interpolate.interp1d=self.dGzdz,
                     beam_current:float=self.current,
                     emitt_n_x:float=self.normalized_emittans_x,
                     emitt_n_y:float=self.normalized_emittans_y,
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

            def dXdz_particle(X:list, z:np.arange,
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

            sol = odeint(dXdz, X0, self.parameter, rtol=rtol)
            self.envelope_x = sol[0:,0]
            self.beam_angular_x = sol[0:,1]
            self.envelope_y = sol[0:,2]
            self.beam_angular_y = sol[0:,3]

            for particle in self.particles.values():
                X0_particle = [particle.x_position, particle.angular_x, particle.y_position, particle.angular_y, particle.phase]
                sol = odeint(dXdz_particle, X0_particle, self.parameter, rtol=rtol)
                self.particle_larmor_x[particle.name] =  sol[0:,0]
                self.particle_larmor_angular_x[particle.name] =  sol[0:,1]
                self.particle_larmor_y[particle.name] =  sol[0:,2]
                self.particle_larmor_angular_y[particle.name] =  sol[0:,3]
                self.particle_larmor_phase[particle.name] =  sol[0:,4]

                self.particle_x[particle.name] = sol[0:,0]*np.cos(sol[0:,4]) - sol[0:,2]*np.sin(sol[0:,4])
                self.particle_y[particle.name] = sol[0:,0]*np.sin(sol[0:,4]) + sol[0:,2]*np.cos(sol[0:,4])
                self.particle_angular_x[particle.name] = sol[0:,1]*np.cos(sol[0:,4]) - sol[0:,3]*np.sin(sol[0:,4])
                self.particle_angular_y[particle.name] = sol[0:,1]*np.sin(sol[0:,4]) + sol[0:,3]*np.cos(sol[0:,4])

                self.particle_radius[particle.name] = (sol[0:,0]**2 + sol[0:,2]**2)**0.5

        elif solver == 'solve_ivp':

            def dXdz(z:np.arange, X:list,
                     dz:float=self.step,
                     gamma:interpolate.interp1d=self.gamma,
                     Bz:interpolate.interp1d=self.Bz,
                     Ez:interpolate.interp1d=self.Ez,
                     Gz:interpolate.interp1d=self.Gz,
                     dBzdz:interpolate.interp1d=self.dBzdz,
                     dEzdz:interpolate.interp1d=self.dEzdz,
                     dGzdz:interpolate.interp1d=self.dGzdz,
                     beam_current:float=self.current,
                     emitt_n_x:float=self.normalized_emittans_x,
                     emitt_n_y:float=self.normalized_emittans_y,
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

            def dXdz_particle(z:np.arange, X:list,
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

            X0 = np.array(X0)
            sol = solve_ivp(dXdz, t_span=[self.parameter[0],self.parameter[-1]], y0=X0, t_eval=self.parameter, rtol=rtol).y
            self.envelope_x = sol[0,:]
            self.beam_angular_x = sol[1,:]
            self.envelope_y = sol[2,:]
            self.beam_angular_y = sol[3,:]

            for particle in self.particles.values():
                X0_particle = np.array([particle.x_position, particle.angular_x, particle.y_position, particle.angular_y, particle.phase])
                sol = solve_ivp(dXdz_particle, t_span=[self.parameter[0],self.parameter[-1]], y0=X0_particle, t_eval=self.parameter, rtol=rtol).y
                self.particle_larmor_x[particle.name] = sol[0,:]
                self.particle_larmor_angular_x[particle.name] =  sol[1,:]
                self.particle_larmor_y[particle.name] =  sol[2,:]
                self.particle_larmor_angular_y[particle.name] =  sol[3,:]
                self.particle_larmor_phase[particle.name] =  sol[4,:]

                self.particle_x[particle.name] = sol[0,:]*np.cos(sol[4,:]) - sol[2,:]*np.sin(sol[4,:])
                self.particle_y[particle.name] = sol[0,:]*np.sin(sol[4,:]) + sol[2,:]*np.cos(sol[4,:])
                self.particle_angular_x[particle.name] = sol[1,:]*np.cos(sol[4,:]) - sol[3,:]*np.sin(sol[4,:])
                self.particle_angular_y[particle.name] = sol[1,:]*np.sin(sol[4,:]) + sol[3,:]*np.cos(sol[4,:])


Sim = Simulation
