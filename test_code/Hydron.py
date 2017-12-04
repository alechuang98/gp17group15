from vpython import *

scene3=canvas(width=800,height=800,background=color.white)

proton={"M": 1.672621898e-27, "Q": 1.6021766208e-19, "R": 1e-11}
electron={"M": 9.10938356e-31, "Q": -1.6021766208e-19, "R": 0.25e-11}
r={"s1": 5.2917721067e-11}
k=8.9875517873681764e9

P=sphere(radius=proton["R"],pos=vec(0,0,0),color=color.orange)
E=sphere(radius=electron["R"],pos=vec(r["s1"],0,0),color=color.blue,make_trail=True)
V=vec(0,sqrt(-k*proton["Q"]*electron["Q"]/(electron["M"]*r["s1"])),0)
pE=V*electron["M"]
print(V,round(pE.y,31))

dt=3*1e-19
t=0
while t<1e-15:
	rate(100)
	pE+=k*proton["Q"]*electron["Q"]/mag(E.pos-P.pos)**2*norm(E.pos-P.pos)*dt
	V=pE/electron["M"]
	E.pos+=V*dt
	t+=dt
