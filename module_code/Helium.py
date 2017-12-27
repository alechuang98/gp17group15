from vpython import *
from random import random
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
	def __init__(self,_pos,_color,_qcenter,_angle):
		self.s=sphere(radius=0.25e-12*4,pos=_pos,color=_color,make_trail=True,interval=10)
		self.q=-eV
		self.m=9.10938356e-31
		a=_pos.x
		b=_pos.y
		self.p=vec(-b/sqrt(a**2+b**2)*sin(angle),a/sqrt(a**2+b**2)*sin(angle),cos(angle))*sqrt(-k*_qcenter*-eV/(self.m*mag(_pos)))*self.m

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

p1=proton(color.red,2)

e=[]
for i in range(2):
	angle=random()*2*pi
	e.append(electron(vec(r["s1"]/2*cos(angle),r["s1"]/2*sin(angle),0),color.white,p1.q,random()*2*pi))

dt=0.5*1e-19/3
t=0

while mag(e[0].s.pos)>2*1e-12:
	rate(1000)
	for i in range(2):
		e[i].emmit(p1.q,dt)
		e[i].update(p1.q,dt)
	t+=dt
