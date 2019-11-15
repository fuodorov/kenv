# beam.py

'''Creating an electron beam.
'''

__all__ = ['Particle', 'Beam']


class Particle:
    '''Sets particle in the beam.

    Sets particle in the beam with parameters:
    x_position [m],
    y_position [m],
    phase [rad],
    name --- unique element name.
    '''
    def __init__(self,
                 x_position: float,
                 y_position: float,
                 angular_x: float,
                 angular_y: float,
                 phase: float,
                 name: str):
        self.x_position = x_position
        self.y_position = y_position
        self.angular_x = angular_x
        self.angular_y = angular_y
        self.phase = phase
        self.name = name

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
    normalized emittans y [m*rad].
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
                 normalized_emittans_y: float=.0e0):

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

    def add_particle(self,
                     name: str,
                     x_position: float,
                     y_position: float,
                     angular_x: float,
                     angular_y: float,
                     phase: float=.0e0) -> None:
        '''Create particle in the beam.

        Create particle in the beam with parameters:
        x_position [m],
        y_position [m],
        angular_x [rad],
        angular_y [rad],
        phase [rad],
        name --- unique element name.
        '''
        self.particles[name] = Particle(x_position, y_position,\
         angular_x, angular_y,\
          phase, name)

    def delete_particle(self,
                     name: str='all') -> None:
        '''Delete particle in the beam.

        Delete a particle in the beam with parameters:
        name --- qparticle's id

        *if there is no name that will be removed all
        '''
        if name == 'all':
            self.particles = {}
        else:
            self.particles.pop(name)

    def __str__(self):
            string=''
            string += '\tParticles:\n'
            for particle in self.particles.values():
                string +="\t[ %.1f mm, %.1f mm, %.1f mrad, %.1f mrad, '%s'],\n" % \
                (particle.x_position*1e3, particle.y_position*1e3,\
                 particle.angular_x*1e3, particle.angular_y*1e3,\
                  particle.name)
            return 'Beam parameters:' + '\n' \
                    +'\tCurrent\t%0.1f A'%(self.current) + '\n' \
                    +'\tEnergy\t%0.1f MeV'%(self.energy) + '\n' \
                    +'\tRadius x\t%0.1f mm'%(self.radius_x*1e3) + '\n' \
                    +'\tRadius y\t%0.1f mm'%(self.radius_y*1e3) + '\n' \
                    +'\tAngular x\t%0.1f mrad'%(self.angular_x*1e3) + '\n' \
                    +'\tAngular y\t%0.1f mrad'%(self.angular_y*1e3) + '\n' \
                    +'\tNormalized emittans x\t%0.1f mm*mrad'%\
                    (self.normalized_emittans_x*1e6) + '\n' \
                    +'\tNormalized emittans y\t%0.1f mm*mrad'%\
                    (self.normalized_emittans_y*1e6) + '\n' \
                    +string

    add_part = add_new_particle = add_particle
    del_part = del_particle = delete_particle
