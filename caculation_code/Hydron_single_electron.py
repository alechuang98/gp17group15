from vpython import *

'''
using larmor formula to caculate the deltaR in each deltaT
warning: this program use fucking long long long long long long time
'''

scene=canvas(width=800,height=800,background=color.white)

eV=1.6021766208e-19
k=8.9875517873681764e9
c=299792458
e0=8.854187817e-12

r={"s1": 5.2917721067e-11}

class electron(object):
    def __init__(self,_pos,_color,_qcenter):
        self.s=sphere(radius=0.25e-11,pos=_pos,color=_color)
        self.q=-eV
        self.m=9.10938356e-31
        self.p=vec(0,sqrt(-k*_qcenter*-eV/(self.m*mag(_pos))),0)*self.m

    def update(self,_qcenter,dt):
        self.p+=k*_qcenter*self.q/mag(self.s.pos)**2 *norm(self.s.pos)*dt
        self.s.pos+=self.p/self.m *dt

    def emmit(self,_qcenter,dt):
        dr = self.q**4 *dt / (12 * 2 * pi**2 * c**3 * e0**2 * self.m**2 * mag(self.s.pos)**2)
        self.s.pos = (mag(self.s.pos)-dr)*norm(self.s.pos)

class proton(object):
    def __init__(self,_color,_q):
        self.s=sphere(radius=1e-11,pos=vec(0,0,0),color=_color)
        self.q=_q*eV
        self.m=1.672621898e-27

p1=proton(color.orange,1)
e1=electron(vec(r["s1"],0,0),color.blue,p1.q)

dt=1*1e-19
t=0

while mag(e1.s.pos)>2.1*1e-15:
    if(mag(e1.s.pos)<1e-12):
        dt=1e-22
    if(mag(e1.s.pos)<1e-13):
        dt=1e-27
    if(mag(e1.s.pos)<1e-14):
        dt=1e-31
    e1.emmit(p1.q,dt)
    e1.update(p1.q,dt)
    t+=dt
print(t)
#t=1.5e-11
