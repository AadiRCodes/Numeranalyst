import numpy as np
import cmath
pi = np.pi

def fft(x: np.ndarray) -> np.ndarray:
    assert len(x) & (len(x)-1) == 0, "len(x) is not power of 2"
    n = len(x)
    if (n == 1):
        return np.array([x[0]])
    x_even = x[::2]
    x_odd = x[1::2]
    even_fft, odd_fft = fft(x_even), fft(x_odd)
    vals = []
    for k in range(n):
        root = np.cos((2*pi*k)/n)+np.sin((2*pi*k)/n)*1j
        vals.append(even_fft[(2*k % n)//2]+root*odd_fft[(2*k % n)//2])
    return np.array(vals)
    
def hyperceil(n: int) -> int:
    """Find smallest power of 2 greater than n."""
    lg = np.ceil(np.log2(n))
    return pow(2, lg)
  

def poly_mul(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    # Padding
    degp, degq = len(p)-1, len(q)-1
    pow_2 = max(hyperceil(degp+1), hyperceil(degq+1))
    padded_p, padded_q = p+[0]*int(2*pow_2-len(p)), q+[0]*int(2*pow_2-len(q))
    print(padded_p)
    vals_p, vals_q = fft(padded_p), fft(padded_q)
    vals_prod = vals_p*vals_q
    conj = np.conjugate(vals_prod)
    prod = np.conjugate(fft(conj))
    return np.array((1/(2*pow_2))*np.array(prod[0:degp+degq+1]), dtype=np.double)

    
print(poly_mul([1, 1, 1, 1, 1, 1, 1], [11, 12, 13]))