from vpython import *
from random import random

'''
using larmor formula to caculate the deltaR in each deltaT
due to the long time, we let deltaR = deltaR * 5000 
'''

scene=canvas(width=800,height=800,background=vec(29/256,30/256,24/256))
g=graph(xtitle='t(s)', ytitle='R(m)',background=color.black)

eV=1.6021766208e-19
k=8.9875517873681764e9
c=299792458
e0=8.854187817e-12

r={"s1": 5.2917721067e-11/36}

class electron(object):
	def __init__(self,_pos,_color,_qcenter,_angle):
		self.s=sphere(radius=0.25e-12*2,pos=_pos,color=_color,make_trail=False,interval=10,retain=10)
		self.q=-eV
		self.m=9.10938356e-31
		a=_pos.x
		b=_pos.y
		self.p=vec(-b/sqrt(a**2+b**2)*sin(angle),a/sqrt(a**2+b**2)*sin(angle),cos(angle))*sqrt(-k*_qcenter*-eV/(self.m*mag(_pos)))*self.m

	def update(self,_qcenter,dt):
		self.p+=k*_qcenter*self.q/mag(self.s.pos)**2 *norm(self.s.pos)*dt
		self.s.pos+=self.p/self.m *dt

	def emmit(self,_qcenter,dt):
		dr = self.q**4 *dt / (12* pi**2 * c**3 * e0**2 * 2 * self.m**2 * mag(self.s.pos)**2)
		dr*=5000
		self.s.pos = (mag(self.s.pos)-dr)*norm(self.s.pos)

class proton(object):
	def __init__(self,_color,_q):
		self.s=sphere(radius=1e-12,pos=vec(0,0,0),color=_color)
		self.q=_q*eV
		self.m=1.672621898e-27

p1=proton(vec(249/256,38/256,114/256),36)

e=[]
for i in range(2):
	angle=random()*2*pi
	e.append(electron(vec(r["s1"]*cos(angle),r["s1"]*sin(angle),0),vec(174/256,129/256,255/256),p1.q,random()*2*pi))
for i in range(8):
	angle=random()*2*pi
	e.append(electron(vec(r["s1"]*4*cos(angle),r["s1"]*4*sin(angle),0),vec(102/256,217/256,239/256),p1.q,random()*2*pi))
for i in range(18):
	angle=random()*2*pi
	e.append(electron(vec(r["s1"]*9*cos(angle),r["s1"]*9*sin(angle),0),vec(166/256,226/256,42/256),p1.q,random()*2*pi))
for i in range(8):
	angle=random()*2*pi
	e.append(electron(vec(r["s1"]*16*cos(angle),r["s1"]*16*sin(angle),0),vec(253/256,151/256,31/256),p1.q,random()*2*pi))

g1=gcurve(graph=g,color=e[0].s.color)
g2=gcurve(graph=g,color=e[2].s.color)
g3=gcurve(graph=g,color=e[10].s.color)
g4=gcurve(graph=g,color=e[28].s.color)

dt=0.5*1e-22
t=0

while mag(e[35].s.pos)>1*1e-12:
	rate(1000)
	f=0
	for i in range(36):
		if mag(e[i].s.pos)<0.6*1e-12:
			if mag(e[i].s.pos)>1e-13:
				f=1
			e[i].s.pos=vec(0,0,0)
			continue
		e[i].emmit(p1.q,dt)
		e[i].update(p1.q,dt)

	g1.plot([t,mag(e[0].s.pos)])
	g2.plot([t,mag(e[2].s.pos)])
	g3.plot([t,mag(e[10].s.pos)])
	g4.plot([t,mag(e[28].s.pos)])

	if f==1 and dt<0.5e-20:
		dt*=5
	t+=dt
for i in range(8):
	e[i+28].s.pos=vec(0,0,0)