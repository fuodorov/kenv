# beam.py

'''Creating an electron beam.
'''

__all__ = ['Beam']

class Beam:
    '''Creating an electron beam.
    
    Creating an elliptical electron beam with parameters:
    current [A],
    energy [MeV],
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
        self.radius_x = radius
        self.radius_y = radius
        self.angular_x = angular
        self.angular_y = angular
        self.normalized_emittans_x = normalized_emittans
        self.normalized_emittans_y = normalized_emittans
        if (radius_x or radius_y):
            self.radius_x = radius_x
            self.radius_y = radius_y   
        if (angular_x or angular_y):
            self.angular_x = angular_x
            self.angular_y = angular_y
        if (normalized_emittans_x or normalized_emittans_y):
            self.normalized_emittans_x = normalized_emittans_x
            self.normalized_emittans_y = normalized_emittans_y
            
    
    def __str__(self):
            return 'Beam parameters:' + '\n' \
                    +'\tCurrent\t%0.1f A'%(self.current) + '\n' \
                    +'\tEnergy\t%0.1f MeV'%(self.energy) + '\n' \
                    +'\tRadius x\t%0.1f mm'%(self.radius_x*1e3) + '\n' \
                    +'\tRadius y\t%0.1f mm'%(self.radius_y*1e3) + '\n' \
                    +'\tAngular x\t%0.1f mrad'%(self.angular_x*1e3) + '\n' \
                    +'\tAngular y\t%0.1f mrad'%(self.angular_y*1e3) + '\n' \
                    +'\tNormalized emittans x\t%0.1f mm*mrad'%(self.normalized_emittans_x*1e6) + '\n' \
                    +'\tNormalized emittans y\t%0.1f mm*mrad'%(self.normalized_emittans_y*1e6) + '\n'