# beam.py

'''Creating an electron beam.
'''

__all__ = ['Beam']

class Beam:
    '''Creating an electron beam.

    Creating a round electron beam with parameters:
    current [A],
    energy [MeV],
    radius [m],
    rp [rad],
    normalized emittance [m*rad],

    Creating an elliptical electron beam with parameters:
    radius x [m],
    radius y [m],
    radius_xp [rad],
    radius_yp [rad],
    normalized emittance x [m*rad],
    normalized emittance y [m*rad]

    and with shifted centroid:
    x [m], y [m], xp[rad], yp[rad], larmor_angle [rad]
    '''

    def __init__(self,
                 current: float=.0e0,
                 energy: float=.0e0,
                 radius: float=.0e0,
                 radius_x: float=.0e0,
                 radius_y: float=.0e0,
                 rp: float=.0e0,
                 radius_xp: float=.0e0,
                 radius_yp: float=.0e0,
                 normalized_emittance: float=.0e0,
                 normalized_emittance_x: float=.0e0,
                 normalized_emittance_y: float=.0e0,
                 x: float=.0e0,
                 y: float=.0e0,
                 xp: float=.0e0,
                 yp: float=.0e0,
                 larmor_angle: float=.0e0):

        self.current = current
        self.energy = energy
        self.radius = radius
        self.rp = rp
        self.radius_x = radius_x
        self.radius_y = radius_y
        self.radius_xp = radius_xp
        self.radius_yp = radius_yp
        self.normalized_emittance = normalized_emittance
        self.normalized_emittance_x = normalized_emittance_x
        self.normalized_emittance_y = normalized_emittance_y
        if radius != .0e0:
            self.radius_x = radius
            self.radius_y = radius
        if rp !=.0e0:
            self.radius_xp = rp
            self.radius_yp = rp
        if normalized_emittance !=.0e0:
            self.normalized_emittance_x = normalized_emittance
            self.normalized_emittance_y = normalized_emittance

        self.x = x
        self.y = y
        self.xp = xp
        self.yp = yp
        self.larmor_angle = larmor_angle

    def __str__(self):
            return 'Beam parameters:' + '\n' \
                    +'\tCurrent\t%0.1f A'%(self.current) + '\n' \
                    +'\tEnergy\t%0.1f MeV'%(self.energy) + '\n' \
                    +'\tRadius x\t%0.1f mm'%(self.radius_x*1e3) + '\n' \
                    +'\tRadius y\t%0.1f mm'%(self.radius_y*1e3) + '\n' \
                    +'\tRadius x prime\t%0.1f mrad'%(self.radius_xp*1e3) + '\n' \
                    +'\tRadius y prime\t%0.1f mrad'%(self.radius_yp*1e3) + '\n' \
                    +'\tHorizontal centroid position\t%0.1f mm'%(self.x*1e3) + '\n' \
                    +'\tVerical centroid position\t%0.1f mm'%(self.y*1e3) + '\n' \
                    +'\tHorizontal centroid angle\t%0.1f mm'%(self.xp*1e3) + '\n' \
                    +'\tVerical centroid angle\t%0.1f mm'%(self.yp*1e3) + '\n' \
                    +'\tLarmor angle\t%0.1f mrad'%(self.larmor_angle*1e3) + '\n' \
                    +'\tNormalized emittance x\t%0.1f mm*mrad'%\
                    (self.normalized_emittance_x*1e6) + '\n' \
                    +'\tNormalized emittance y\t%0.1f mm*mrad'%\
                    (self.normalized_emittance_y*1e6) + '\n' \
