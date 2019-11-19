from numpy import sqrt, square, divide, inf
from scipy.optimize import curve_fit

class dROverR:
    def __init__(self, i=None, c1_2=None, c3=None, c4=None, c5=None):
        self.i=i
        self.c1_2=c1_2
        self.c3=c3
        self.c4=c4
        self.c5=c5

    def __call__(self, q, R, f=None):
        if f is None:
            f = self.func_dR_R1
        dR = True
        return dR

    def new(self, i=None, c1_2=None, c3=None, c4=None, c5=None):
        self.i=i
        self.c1_2=c1_2
        self.c3=c3
        self.c4=c4
        self.c5=c5
        return True

    def s1_2(self,q):
        return self.c1_2*q # q

    def t(self,q):
        return self.c3 + self.c4*square(q) # q

    def nb(self, q, R):
        return ( self.i*self.c5*
                self.t(q)*
                square(self.s1_2(q)) )

    def ns(self, q, R):
        return ( self.i*R*
                self.t(q)*
                square(self.s1_2(q)) )

    def func_dR_R1(self, q, R):
        out = divide( sqrt(R+2*c5),
                     R*sqrt(i*self.t(q)*
                                 square( self.s1_2(q) )) )
        return out

    def func_dR_R2(self, q, R):
        ns = self.ns(q,R)
        nb = self.nb(q,R)
        return divide( sqrt( ns + 2*nb), ns )

    def dR_R_func(self, q, R): # difrent
        ns = self.ns(q,R)
        return divide( sqrt(ns+2*self.c5), ns)



class minimiser:
    def __init__(self, classToMini, initial):
        print(initial)
        self.func = classToMini(*initial)

    def __call__(self):
        pass









