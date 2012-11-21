#!/usr/bin/python
# -*- coding: iso-8859-1 -*-

"Some things cannot be calculated. For everything else, you've got maths."

## @package pyspline
#
# @todo functions define functions in order to use the module
#        without classes for some users
# @todo addition work on mono+int, mono+mono, basis+mono
# @todo plots work on plot overloading
#
# @version 0.2
# @author P.E.Bontinck
# @licence GPL3

import scipy as sc
import pylab as pl


# utility functions

def convarg(arg):
    """
    convert any arguments in flat mat of float64
    1 ==> mat([[1.]])
    [1] ==> mat([1.]])
    [1,1] ==> mat([[1.,1.]])
    etc...
    """
    return(sc.mat(arg, dtype='f8').flatten())

def convarga(arg):
    """
    convert any arguments in float array of float64
    1 ==> array([1.])
    [1] ==> array([1.])
    [1.] ==> array(([1.])
    """
    return(sc.array(arg, dtype='f8').flatten())

# generic classes

class Mono(object):
    pass

class Basis(object):
    def __init__(self, t):
        """
        basis of mono functions
        t: array of mono functions
        """
        self.t = t
        self.lt = len(t)

    def __call__(self, dx):
        """
        basis evaluation with vector x
        """
        x = convarg(dx)
        y = sc.zeros((x.shape[1], self.lt))
        for i in range(self.lt):
            y[:,i] = self.t[i](x)
        return sc.mat(y)

    def __repr__(self):
        return "basis :{0}".format(self.t)


class genspline(object):
    def __init__(self, basis, coeffs):
        """
        generic linear combination from a basis of functions and coefficients for each function
        """
        self.b = basis
        self.c = convarg(coeffs)
        if self.b.lt!=self.c.shape[1]:
            raise NameError("invalid sizes between basis and coefficients.")

    def __call__(self, x):
        X = self.b(x)
        return(X*self.c)

class optspline(genspline):
    def __init__(self, basis, x, y):
        """
        spline fitting of points (x,y) on basis
        """
        self.b = basis
        X = self.b(x)
        self.c = (X.T*X).I*X.T*convarg(y).T


# specific classes

class hat(Mono):
    def __init__(self, left=-1.0, right=1.0, center=0.0):
        """
        simple normalized hat function
        """
        self.a = left - center
        self.b = right - center
        self.c = center
        if self.a>=.0:
            self.a = -sc.inf
        if self.b<=.0:
            self.b = sc.inf

    def __call__(self, dx):
        """
        hat function evaluation
        """
        x = convarg(dx) - self.c
        y = sc.zeros_like(x)
        i = sc.nonzero((x>self.a)&(x<=.0))
        y[i] = 1. - x[i]/self.a
        i = sc.nonzero((x>.0)&(x<self.b))
        y[i] = 1. - x[i]/self.b
        return(y)

    def __repr__(self):
        return "hat function ({0},{1},{2})".format(self.a+self.c,
                                                   self.c, self.b+self.c)

class hatbasis(Basis):
    def __init__(self, bp):
        """
        hat functions orthogonal basis based on breakpoints bp
        @todo growing not constant array is of interest ?
        """
        bpa = convarga(bp)
        x = sc.hstack((bpa.min(), bpa, bpa.max()))
        self.t = []
        for i in range(len(bp)):
            self.t.append(hat(x[i], x[i+2], x[i+1]))
        self.lt = len(self.t)

class absf(Mono):
    def __init__(self, c):
        self.c = c
    def __call__(self, x):
        return(abs(x-self.c))
    def __repr__(self):
        return("|x-{0}|".format(self.c))

class absbasis(Basis):
    def __init__(self, bp):
        bpa = convarga(bp)
        self.t = []
        for i in bp:
            self.t.append(absf(i))
        self.lt = len(self.t)

class mpoly(Mono):
    def __init__(self, order):
        self.n = order

    def __call__(self, x):
        return(sc.power(x, self.n))

    def __repr__(self):
        return "x^{0}".format(self.n)

class polybasis(Basis):
    def __init__(self, n):
        """
        polynomial basis of order n
        """
        self.t = []
        for i in range(n+1):
            self.t.append(mpoly(i))
        self.lt = len(self.t)

class absf3(Mono):
    def __init__(self, c):
        self.c = c
    def __call__(self, x):
        return(sc.power(abs(x-self.c), 3))
    def __repr__(self):
        return("|x-{0}|^3".format(self.c))

class cubsplinebasis(Basis):
    def __init__(self, bp):
        self.t = [mpoly(0),mpoly(1)]
        bpa = convarga(bp)
        for i in bpa:
            self.t.append(absf3(i))
        self.lt = len(self.t)

if __name__ == "__main__":
    """ module self testing
    """
    x = sc.array((0, 0.1, 1, 1.5, 2, 2.3, 2.6, 2.9, 3))
    xp = sc.linspace(-.5, 3.5, 200)
    y = sc.sin(x/1.5)
    b0 = hatbasis((0, 1, 2, 3))
    y0 = optspline(b0, x, y)
    print("hat optimisation OK")
    b1 = polybasis(2)
    y1 = optspline(b1, x, y)
    print("polynomial optimisation OK")
    b2 = absbasis((0, 1, 2, 3))
    y2 = optspline(b2, x, y)
    print("abs functions optimisation OK")
    b3 = cubsplinebasis((0, 1, 2, 3))
    y3 = optspline(b3, x, y)
    print("cubic spline optimisation OK")

    pl.ioff()
    pl.subplot(2,1,1)
    pl.plot(x, y, 'bo-')
    pl.plot(xp, y0(xp), 'g-')
    print("hat evaluation OK")
    pl.plot(xp, y1(xp), 'r-')
    print("polynomial evaluation OK")
    pl.plot(xp, y2(xp), 'm-')
    print("abs evaluation OK")
    pl.plot(xp, y3(xp), 'c-')
    print("cubic evaluation OK")
    pl.grid()
    pl.subplot(2,1,2)
    pl.plot(x, y-y0(x), 'g')
    pl.plot(x, y-y1(x), 'r')
    pl.plot(x, y-y2(x), 'm')
    pl.plot(x, y-y3(x), 'c')
    pl.show()
