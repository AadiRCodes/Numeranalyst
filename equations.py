import numpy as np
from Func import Func
def bisection_method(f: Func, a: float, b: float, tol=10**(-9), max_iter = 100):
    """Calculates a value c in [a, b] such that f(c) = 0 for a function f: R -> R, 
    given f(a)*f(b)<=0"""
    assert f.eval(a)*f.eval(b) <= 0, "Method only works if f(a) and f(b) have opposite signs!"
    if f.eval(a) == 0:
        return a
    if f.eval(b) == 0:
        return b
    i = 0
    left, right = a, b
    left_val = f.eval(left)
    while i < max_iter:
        mdpt = left+(right-left)/2 # Avoids overflow errors
        md_val = f.eval(mdpt)
        if md_val == 0 or (right-left)/2 < tol:
            return mdpt
        if left_val * md_val > 0:
            left = mdpt
        else:
            right = mdpt
        i+=1
    if i == max_iter:
        print("Method failed, surpassed iteration limit")
        pass

f = Func(lambda x: x**2-3)
print(bisection_method(f, 0, 2, tol=10**(-12)))


def newton_method(f: Func, init: float, tol=10**(-10), max_iter = 100):
    """Uses Newton's Method to find a solution to f(x)=0. Note that
    an initial approximation init is required."""
    i = 0
    root = init
    while i < max_iter:
        deriv = f.derivative(init)
        val = f.eval(init)
        root = init-(val/deriv)
        if (abs(root-init) < tol):
            return root
        init = root
        i+=1
    if i == max_iter:
        print(f"Method did not converge in {max_iter} iterations")
        return
    
print(newton_method(f, 1.5, tol=10**(-12)))
g = Func(lambda x: np.exp(-x)-x)
print(newton_method(g, 1.5, tol=10**(-12)))