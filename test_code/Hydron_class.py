from vpython import *

scene=canvas(width=800,height=800,background=color.white)

eV=1.6021766208e-19
k=8.9875517873681764e9

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

class proton(object):
	def __init__(self,_color,_q):
		self.s=sphere(radius=1e-11,pos=vec(0,0,0),color=_color)
		self.q=_q*eV
		self.m=1.672621898e-27

p1=proton(color.orange,1)
e1=electron(vec(r["s1"],0,0),color.blue,p1.q)

dt=3*1e-19
t=0

while t<1e-15:
	rate(100)
	e1.update(p1.q,dt)
	t+=dt
