# Problem set 1

Started on July $16^{th}$, 2020. These problems are a simple review of the basic concepts in special relativity.

* 3.1 (pg. 45)
> A free particle is moving in an inertial frame (x,y,z) in the x-y plane on a trajectory x = d, y = ut, where d and v are constants in time. Consider a rectangular frame (x',y',z') rotating with respect to the inertial frame with an angular velocity $\omega$ about a common z-axis (z'= z). What are the equations of motion obeyed by x'(t), y'(t), and z'(t) in the rotating frame? Sketch the trajectory of the particle in the x'-y' plane and show explicitly that it satisfies these equations of motion.

$$\frac{d}{dt}\frac{\partial L}{\partial q}-\frac{\partial L}{\partial \dot{q}}=0$$
$$x'=cos(\omega t)(x-d)+sin(\omega t)(y-ut)$$
$$y'=-sin(\omega t)(x-d)+cos(\omega t)(y-ut)$$

$$\dot{x}'=cos(\omega t)\dot{x}-sin(\omega t)\omega x+sin(\omega t)\omega d+cos(\omega t)\omega y+sin(\omega t)\dot{y}-cos(\omega t)\omega ut-sin(\omega t)u$$
$$\dot{x}'=cos(\omega t)[\dot{x}+\omega y-\omega ut]+sin(\omega t)[\dot{y}+\omega d-\omega x-u]$$
$$\dot{y}'=-sin(\omega t)\dot{x}-cos(\omega t)\omega x+cos(\omega t)\omega d+cos(\omega t)\dot{y}-sin(\omega t)\omega y+sin(\omega t)\omega ut-cos(\omega t)u$$
$$\dot{y}'=cos(\omega t)[\dot{y}+\omega d-\omega x-u]+sin(\omega t)[\omega ut-\dot{x}-\omega y]$$
<br>
$$L=\frac{1}{2}mv^2=\frac{1}{2}m\sqrt{\dot{x}^{2'}+\dot{y}^{2'}}$$
$$L=\frac{1}{2}m\sqrt{(cos(\omega t)[\dot{x}+\omega y-\omega ut]+sin(\omega t)[\dot{y}+\omega d-\omega x-u])^{2'}+(cos(\omega t)[\dot{y}+\omega d-\omega x-u]+sin(\omega t)[\omega ut-\dot{x}-\omega y])^{2'}}$$

* x'
$$\frac{d}{dt}\frac{\partial L}{\partial x}=\frac{d}{dt}\frac{\partial}{\partial x}[\frac{1}{2}m\sqrt{(cos(\omega t)[\dot{x}+\omega y-\omega ut]+sin(\omega t)[\dot{y}+\omega d-\omega x-u])^{2'}+(cos(\omega t)[\dot{y}+\omega d-\omega x-u]+sin(\omega t)[\omega ut-\dot{x}-\omega y])^{2'}}]$$
$$\frac{\partial L}{\partial \dot{x}}=\frac{\partial }{\partial \dot{x}}[\frac{1}{2}m\sqrt{(cos(\omega t)[\dot{x}+\omega y-\omega ut]+sin(\omega t)[\dot{y}+\omega d-\omega x-u])^{2'}+(cos(\omega t)[\dot{y}+\omega d-\omega x-u]+sin(\omega t)[\omega ut-\dot{x}-\omega y])^{2'}}]$$
* y'
$$\frac{d}{dt}\frac{\partial L}{\partial x}=\frac{d}{dt}\frac{\partial}{\partial y}[\frac{1}{2}m\sqrt{(cos(\omega t)[\dot{x}+\omega y-\omega ut]+sin(\omega t)[\dot{y}+\omega d-\omega x-u])^{2'}+(cos(\omega t)[\dot{y}+\omega d-\omega x-u]+sin(\omega t)[\omega ut-\dot{x}-\omega y])^{2'}}]$$
$$\frac{\partial L}{\partial \dot{x}}=\frac{\partial }{\partial \dot{y}}[\frac{1}{2}m\sqrt{(cos(\omega t)[\dot{x}+\omega y-\omega ut]+sin(\omega t)[\dot{y}+\omega d-\omega x-u])^{2'}+(cos(\omega t)[\dot{y}+\omega d-\omega x-u]+sin(\omega t)[\omega ut-\dot{x}-\omega y])^{2'}}]$$
* $z'(t)=z$


Stopped after creating L, may use code to solve.


```python
from sympy import *
import numpy as np

m = "mass"
w = "omega"
d = "x_constant"
u = "y_constant"

L = (1/2)*m*sqrt((-sin(w*t)*w*d+cos(w*t)w*u*t+sin(w*t)*u)**2+(-cos(w*t)*w*d-sin(w*t)w*u*t+cos(w*t)*u)**2)

dt_pL = [1,1]
```


      File "<ipython-input-1-11096b829e89>", line 9
        L = (1/2)*m*sqrt((-sin(w*t)*w*d+cos(w*t)w*u*t+sin(w*t)*u)**2+(-cos(w*t)*w*d-sin(w*t)w*u*t+cos(w*t)*u)**2)
                                                ^
    SyntaxError: invalid syntax



* 3.5 (pg. 46)
> Conside the functional $$S[x(t)] = \int_{0}^{T} [(\frac{dx(t)}{dt})^2 + x^2(t)]dt.$$ <br>
Find the curve x(t) satisfying the conditions $$x(0) = 0, x(T) = 1,$$ <br>
which makes S[x(t)] an extremum. What is the extremum value of S[x(t)]? Is it a maximum or minimum?

* 4.2 (pg. 73)
> A rocket ship of proper length L leaves Each vertically at speed $\frac{4}{5}c$. A light signal is sent vertically after it which arrives at the rocket's tail at t = 0 according to both rocket and Earth-based clocks. When does the signal reach the nose of the rocket according to <br>
(a) the rocket clocks <br>
(b) the Earth clocks

* 4.4 (pg. 73)
> A satellite orbits the Earth in the same direction it rotates in a circular orbit above the equator a distance 200 km from the surface. By how many seconds per day will a clock on such a satellite run slow compared to a clock on the Earth? (Compute just the special relativistic effects.)

* 4.5 (pg. 74)
> The radio source 3C345 is participating in the expansion of the universe, and its distance can be determined from the redshift arising from its recession velocity and assumptions about our universe. However, a rough idea of the distance can be obtained from Hubble's law relating distance d to observed recession velocity V: $$V = H_0d,$$ <br>
where $H_0\approx 72 (km/s)/Mpc$ is the Hubble constant. (Look at the endpapers for astronomical units such as megaparsec (Mpc).) V for 3C345 is about 0.6c. Use these facts together with the data in Box 4.3 on p. 61 to roughly estimate the velocity of the cloud C2 assuming (contrary to fact) that it is moving transverse to the line of sight.

* 4.8 (pg. 74)
> Calculate the hyperbolic angle betweeen the sides AC and AB of triangle ABC illustrated in Figure 4.8.

* 4.12 (pg. 75)
> (a) Show explicitly that the straight line path between any two points in flat three-dimensional space ($dS^2=dx^2+dy^2+dz^2$) is the shortest distance between them. <br>
(b) Is the straight line path between two spacelike separated points in flat spacetime the shortest distance between them?

* 4.18 (pg. 76)
> Show that for two timelike separated events, there is some inertial frame in which $\Delta t \neq 0,\Delta \vec{x} = 0$. Show that for two spacelike separated events there is an inertial frame where $\Delta t = 0,\Delta \vec{x} \neq 0$.
