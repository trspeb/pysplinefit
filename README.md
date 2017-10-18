PYSPLINEFIT
===========

Fitting a general spline to data in python/scipy.

# Example : Study of the handspinner speed with arduino and python/scipy

---

## Method

A proximity detector realized with a LED, a photoresistor and an arduino board is realized. The duration between pulses is computed and RPM is calculated.

Data are processed in order to define several parameters of a physical model.

## Physical model

We want to express the rotation speed $N$ of the handspinner.

Inertia $I$ is storing the rotational energy and a friction torque $T$ from air and bearing is removing energy.

$$T = I \frac{dN}{dt}$$

We can estimate that that friction torque $T$ will be a function of the speed $N$. This function shall be decreasing, as friction with air and brearing should be higher at higher speed. 

There is no acceleration due to friction when $N=0$, hence $T(0)=0$.

## Mathematical model

We have:

$\frac{dN}{dt} = f(N)$ with $f(0)=0$.

## Numerical test

We realize a first experimentation. We calculate the speed from intervals between detection pulses.

![experimental data](res/figure_1-1.svg)



