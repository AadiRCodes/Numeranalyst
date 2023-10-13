import numpy as np
import matplotlib.pyplot as plt
def clamped_spline(x_vals: np.ndarray, func_vals: np.ndarray, clamp1, clamp2):
    assert len(x_vals) == len(func_vals), "x-values and function values do not match!"
    n = len(x_vals)-1
    alpha = [0]*(n+1)
    steps = [0]*n
    for i in range(n):
        steps[i] = x_vals[i+1]-x_vals[i]
    alpha[0] = 3*(func_vals[1]-func_vals[0])/steps[0] - 3*clamp1
    alpha[n] = 3*clamp2 - 3*(func_vals[n]-func_vals[n-1])/steps[n-1]
    for i in range(1, n):
        alpha[i] = (3/steps[i])*(func_vals[i+1]-func_vals[i])-(3/steps[i-1])*(func_vals[i]-func_vals[i-1])
    l, mu, z = [0]*(n+1), [0]*n, [0]*(n+1)
    l[0], mu[0] = 2*steps[0], 0.5
    z[0] = alpha[0]/l[0]
    for i in range(1, n-1):
        l[i] = 2*(x_vals[i+1]-x_vals[i-1])-steps[i-1]*mu[i-1]
        mu[i] = steps[i]/l[i]
        z[i] = (alpha[i]-steps[i-1]*z[i-1])/l[i]
    l[n] = steps[n-1]*(2-mu[n-1])
    z[n] = (alpha[n]-steps[n-1]*z[n-1])/l[n]
    b_vals, c_vals, d_vals = [0]*n, [0]*(n+1), [0]*n
    c_vals[n] = z[n]
    for j in reversed(range(n)):
        c_vals[j] = z[j]-mu[j]*c_vals[j+1]
        b_vals[j] = (func_vals[j+1]-func_vals[j])/steps[j]-(steps[j]/3)*(c_vals[j+1]+2*c_vals[j])
        d_vals[j] = (c_vals[j+1]-c_vals[j])/(3*steps[j])

    polys = zip(func_vals, b_vals, c_vals[:n], d_vals)
    return np.array([list(t) for t in polys])

def spline_cubic(coeffs: np.ndarray, point: float):
    def spline(x):
        return coeffs[0]+coeffs[1]*(x-point)+coeffs[2]*(x-point)**2+coeffs[3]*(x-point)**3
    
    return spline

def plot_spline(x_vals, func_vals, coeffs):
    plt.plot(x_vals, func_vals, 'ro')
    n = len(x_vals)
    for i in range(n-1):
        S = spline_cubic(coeffs[i], x_vals[i])
        x_int = np.arange(x_vals[i], x_vals[i+1], 0.001)
        plt.plot(x_int, S(x_int), "b")
    plt.show()

