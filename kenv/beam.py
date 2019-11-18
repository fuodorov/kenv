# beam.py

'''Creating an electron beam.
'''

__all__ = ['Beam']

class Beam:
    '''Creating an electron beam.

    Creating an elliptical electron beam with parameters:
    current [A],
    energy [MeV],
    radius [m],
    angular [m],
    normalized emittans [m*rad],
    radius x [m],
    radius y [m],
    angular x [rad],
    angular y [rad],
    normalized emittans x [m*rad],
    normalized emittans y [m*rad]
    and
    x [m], y [m], phase [rad]
    '''

    def __init__(self,
                 current: float=.0e0,
                 energy: float=.0e0,
                 radius: float=.0e0,
                 radius_x: float=.0e0,
                 radius_y: float=.0e0,
                 angular: float=.0e0,
                 angular_x: float=.0e0,
                 angular_y: float=.0e0,
                 normalized_emittans: float=.0e0,
                 normalized_emittans_x: float=.0e0,
                 normalized_emittans_y: float=.0e0,
                 x: float=.0e0,
                 y: float=.0e0,
                 phase: float=.0e0):

        self.current = current
        self.energy = energy
        self.radius = radius
        self.radius_x = radius_x
        self.radius_y = radius_y
        self.angular = angular
        self.angular_x = angular_x
        self.angular_y = angular_y
        self.normalized_emittans = normalized_emittans
        self.normalized_emittans_x = normalized_emittans_x
        self.normalized_emittans_y = normalized_emittans_y
        if radius != .0e0:
            self.radius_x = radius
            self.radius_y = radius
        if angular !=.0e0:
            self.angular_x = angular
            self.angular_y = angular
        if normalized_emittans !=.0e0:
            self.normalized_emittans_x = normalized_emittans
            self.normalized_emittans_y = normalized_emittans
        self.particles = {}
        self.x = x
        self.y = y
        self.phase = phase

    def __str__(self):
            return 'Beam parameters:' + '\n' \
                    +'\tCurrent\t%0.1f A'%(self.current) + '\n' \
                    +'\tEnergy\t%0.1f MeV'%(self.energy) + '\n' \
                    +'\tx position\t%0.1f mm'%(self.x*1e3) + '\n' \
                    +'\ty position\t%0.1f mm'%(self.y*1e3) + '\n' \
                    +'\tRadius y\t%0.1f mm'%(self.radius_y*1e3) + '\n' \
                    +'\tRadius x\t%0.1f mm'%(self.radius_x*1e3) + '\n' \
                    +'\tAngular x\t%0.1f mrad'%(self.angular_x*1e3) + '\n' \
                    +'\tAngular y\t%0.1f mrad'%(self.angular_y*1e3) + '\n' \
                    +'\tLarmor phase\t%0.1f mrad'%(self.phase*1e3) + '\n' \
                    +'\tNormalized emittans x\t%0.1f mm*mrad'%\
                    (self.normalized_emittans_x*1e6) + '\n' \
                    +'\tNormalized emittans y\t%0.1f mm*mrad'%\
                    (self.normalized_emittans_y*1e6) + '\n' \
