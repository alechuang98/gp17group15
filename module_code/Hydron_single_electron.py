from vpython import *

'''
using larmor formula to caculate the deltaE in each deltaT
due to the long time, we let deltaE = deltaE * 5000 
'''


scene=canvas(width=800,height=800,background=color.white)

eV=1.6021766208e-19
k=8.9875517873681764e9
c=299792458
e0=8.854187817e-12

r={"s1": 5.2917721067e-11}

class electron(object):
	def __init__(self,_pos,_color,_qcenter):
		self.s=sphere(radius=0.25e-11,pos=_pos,color=_color,make_trail=True)
		self.q=-eV
		self.m=9.10938356e-31
		self.p=vec(0,sqrt(-k*_qcenter*-eV/(self.m*mag(_pos))),0)*self.m

	def update(self,_qcenter,dt):
		self.p+=k*_qcenter*self.q/mag(self.s.pos)**2 *norm(self.s.pos)*dt
		self.s.pos+=self.p/self.m *dt

	def emmit(self,_qcenter,dt):
		dE=(self.q**2 * (k*_qcenter*self.q/mag(self.s.pos)**2 /self.m)**2)/(6*pi*e0*c**3) *dt
		dE*=5000
		self.p=sqrt(mag(self.p)**2-2*self.m*dE)*norm(self.p)
		#print(mag(self.p)**2,2*self.m*dE)
		#self.p+=-(self.q**6/(96 * pi**3 * e0**3 * c**3 * self.m) /(mag(self.s.pos)**4 *dot(self.p,self.s.pos))) *norm(self.s.pos) *dt
		#dr = exp(4) *dt / (12* pi**2 * c**3 * e0**2 * self.m**2 * mag(self.s.pos)**2)
		#print(dr)
		#self.s.pos = (mag(self.s.pos)-dr)*norm(self.s.pos)

class proton(object):
	def __init__(self,_color,_q):
		self.s=sphere(radius=1e-11,pos=vec(0,0,0),color=_color)
		self.q=_q*eV
		self.m=1.672621898e-27

p1=proton(color.orange,1)
e1=electron(vec(r["s1"],0,0),color.blue,p1.q)

dt=3*1e-19
t=0

while t<1e-11 and mag(e1.s.pos)>1e-11:
	rate(100)
	e1.emmit(p1.q,dt)
	e1.update(p1.q,dt)
	t+=dt
