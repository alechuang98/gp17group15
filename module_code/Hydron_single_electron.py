from vpython import *

'''
using larmor formula to caculate the deltaR in each deltaT
due to the long time, we let deltaR = deltaR * 5000 
'''

scene=canvas(width=800,height=800,background=color.black)

eV=1.6021766208e-19
k=8.9875517873681764e9
c=299792458
e0=8.854187817e-12

r={"s1": 5.2917721067e-11}

class electron(object):
	def __init__(self,_pos,_color,_qcenter):
		self.s=sphere(radius=0.25e-12*4,pos=_pos,color=_color,make_trail=True,interval=10)
		self.q=-eV
		self.m=9.10938356e-31
		self.p=vec(0,sqrt(-k*_qcenter*-eV/(self.m*mag(_pos))),0)*self.m

	def update(self,_qcenter,dt):
		self.p+=k*_qcenter*self.q/mag(self.s.pos)**2 *norm(self.s.pos)*dt
		self.s.pos+=self.p/self.m *dt

	def emmit(self,_qcenter,dt):
		dr = self.q**4 *dt / (12* pi**2 * c**3 * e0**2 * self.m**2 * mag(self.s.pos)**2)
		dr*=5000
		self.s.pos = (mag(self.s.pos)-dr)*norm(self.s.pos)

class proton(object):
	def __init__(self,_color,_q):
		self.s=sphere(radius=1e-12*4,pos=vec(0,0,0),color=_color)
		self.q=_q*eV
		self.m=1.672621898e-27

p1=proton(color.red,1)
e1=electron(vec(r["s1"],0,0),color.white,p1.q)

dt=0.5*1e-19
t=0

while mag(e1.s.pos)>2*1e-12:
	rate(1000)
	e1.emmit(p1.q,dt)
	e1.update(p1.q,dt)
	t+=dt
