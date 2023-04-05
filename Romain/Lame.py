#!/usr/bin/python
'''
 - All functions depend only on $\Psi^p_n$, which are expressible in terms of Cartesian coordinates (Table III)
  - Need to know class and parity

'''
from math import sqrt, copysign, pi
import sys
import numpy as np
#try:
#  from pylab import plot, show
#except:
#  print 'Runnning without matplotlib'

# Romberg quadratures for numeric integration.
#
# Written by Scott M. Ransom <ransom@cfa.harvard.edu>
# last revision: 14 Nov 98
#
# Cosmetic changes by Konrad Hinsen <hinsen@cnrs-orleans.fr>
# last revision: 1999-7-21
#
# Adapted to scipy by Travis Oliphant <oliphant.travis@ieee.org>
# last revision: Dec 2001

def _difftrap(function, interval, numtraps):
    """
    Perform part of the trapezoidal rule to integrate a function.
    Assume that we had called difftrap with all lower powers-of-2
    starting with 1.  Calling difftrap only returns the summation
    of the new ordinates.  It does _not_ multiply by the width
    of the trapezoids.  This must be performed by the caller.
        'function' is the function to evaluate (must accept vector arguments).
        'interval' is a sequence with lower and upper limits
                   of integration.
        'numtraps' is the number of trapezoids to use (must be a
                   power-of-2).
    """
    from numpy import arange
    if numtraps <= 0:
        raise ValueError("numtraps must be > 0 in difftrap().")
    elif numtraps == 1:
        return 0.5*(function(interval[0])+function(interval[1]))
    else:
        numtosum = numtraps/2
        h = float(interval[1]-interval[0])/numtosum
        lox = interval[0] + 0.5 * h;
        points = lox + h * arange(0, numtosum)
        s = sum(function(points),0)
        return s

def _romberg_diff(b, c, k):
    """
    Compute the differences for the Romberg quadrature corrections.
    See Forman Acton's "Real Computing Made Real," p 143.
    """
    tmp = 4.0**k
    return (tmp * c - b)/(tmp - 1.0)

def _printresmat(function, interval, resmat):
    # Print the Romberg result matrix.
    i = j = 0
    print 'Romberg integration of', `function`,
    print 'from', interval
    print ''
    print '%6s %9s %9s' % ('Steps', 'StepSize', 'Results')
    for i in range(len(resmat)):
        print '%6d %9f' % (2**i, (interval[1]-interval[0])/(2.**i)),
        for j in range(i+1):
            print '%9f' % (resmat[i][j]),
        print ''
    print ''
    print 'The final result is', resmat[i][j],
    print 'after', 2**(len(resmat)-1)+1, 'function evaluations.'

def vectorize1(func, args=(), vec_func=False):
    """Vectorize the call to a function.

    This is an internal utility function used by `romberg` and
    `quadrature` to create a vectorized version of a function.

    If `vec_func` is True, the function `func` is assumed to take vector
    arguments.

    Parameters
    ----------
    func : callable
        User defined function.
    args : tuple
        Extra arguments for the function.
    vec_func : bool
        True if the function func takes vector arguments.

    Returns
    -------
    vfunc : callable
        A function that will take a vector argument and return the
        result.

    """
    if vec_func:
        def vfunc(x):
            return func(x, *args)
    else:
        def vfunc(x):
            from numpy import asarray, empty, isscalar
            if isscalar(x):
                return func(x, *args)
            x = asarray(x)
            # call with first point to get output type
            y0 = func(x[0], *args)
            n = len(x)
            if hasattr(y0, 'dtype'):
                output = empty((n,), dtype=y0.dtype)
            else:
                output = empty((n,), dtype=type(y0))
            output[0] = y0
            for i in xrange(1, n):
                output[i] = func(x[i], *args)
            return output
    return vfunc

def romberg(function, a, b, args=(), tol=1.48e-8, rtol=1.48e-8, show=False, divmax=10, vec_func=False):
  """
  Romberg integration of a callable function or method.

  Returns the integral of `function` (a function of one variable)
  over the interval (`a`, `b`).

  If `show` is 1, the triangular array of the intermediate results
  will be printed.  If `vec_func` is True (default is False), then `function` is
  assumed to support vector arguments.

  Parameters
  ----------
  function : callable
    Function to be integrated.
  a : float
    Lower limit of integration.
  b : float
    Upper limit of integration.

  Returns
  --------
  results  : float
    Result of the integration.

  Other Parameters
  ----------------
  args : tuple, optional
    Extra arguments to pass to function. Each element of `args` will
    be passed as a single argument to `func`. Default is to pass no
    extra arguments.
  tol, rtol : float, optional
    The desired absolute and relative tolerances. Defaults are 1.48e-8.
  show : bool, optional
    Whether to print the results. Default is False.
  divmax : int, optional
    Maximum order of extrapolation. Default is 10.
  vec_func : bool, optional
    Whether `func` handles arrays as arguments (i.e whether it is a
    "vector" function). Default is False.

  References
  ----------
  .. [1] 'Romberg's method' http://en.wikipedia.org/wiki/Romberg%27s_method

  Examples
  --------
  Integrate a gaussian from 0 to 1 and compare to the error function.

    >>> from scipy.special import erf
    >>> gaussian = lambda x: 1/np.sqrt(np.pi) * np.exp(-x**2)
    >>> result = romberg(gaussian, 0, 1, show=True)
    Romberg integration of <function vfunc at 0x101eceaa0> from [0, 1]

    ::

       Steps  StepSize  Results
           1  1.000000  0.385872
           2  0.500000  0.412631  0.421551
           4  0.250000  0.419184  0.421368  0.421356
           8  0.125000  0.420810  0.421352  0.421350  0.421350
          16  0.062500  0.421215  0.421350  0.421350  0.421350  0.421350
          32  0.031250  0.421317  0.421350  0.421350  0.421350  0.421350  0.421350

    The final result is 0.421350396475 after 33 function evaluations.

    >>> print 2*result,erf(1)
    0.84270079295 0.84270079295

  """
  from numpy import isinf
  if isinf(a) or isinf(b):
    raise ValueError("Romberg integration only available for finite limits.")
  vfunc = vectorize1(function, args, vec_func=vec_func)
  n = 1
  interval = [a,b]
  intrange = b-a
  ordsum = _difftrap(vfunc, interval, n)
  result = intrange * ordsum
  resmat = [[result]]
  err = np.inf
  for i in xrange(1, divmax+1):
    n = n * 2
    ordsum = ordsum + _difftrap(vfunc, interval, n)
    resmat.append([])
    resmat[i].append(intrange * ordsum / n)
    for k in range(i):
      resmat[i].append(_romberg_diff(resmat[i-1][k], resmat[i][k], k+1))
    result = resmat[i][i]
    lastresult = resmat[i-1][i-1]

    err = abs(result - lastresult)
    if err < tol or err < rtol*abs(result):
      break
    else:
      #print("Accuracy Warning: divmax (%d) exceeded. Latest difference = %e" % (divmax, err))
      pass

  if show:
    _printresmat(vfunc, interval, resmat)
  return result

def integrateTrapezoid(f, a, b, steps):
  '''Trapezoid Rule'''
  h = (b-a)/steps
  s = sum(f(a+i*h) for i in range(1, steps))
  return (h/2.0)*(f(a)+f(b)+2.0*s2)

def integrateMidpoint(f, a, b, steps):
  '''Midpoint Rule'''
  h  = (b-a)/steps
  a1 = a+h/2
  s  = sum(f(a1+i*h) for i in range(0, steps))
  return h*s

def integrate(f, a, b, steps):
  '''Simpson's Rule'''
  h  = (b-a)/steps
  a1 = a+h/2
  s1 = sum(f(a1+i*h) for i in range(0, steps))
  s2 = sum(f(a +i*h) for i in range(1, steps))
  return (h/6.0)*(f(a)+f(b)+4.0*s1+2.0*s2)

class EllipsoidalSystem(object):
  def __init__(self, a, b, c = None):
    '''The ellipsoidal coordinate system is defined by three axes a >= b >= c'''
    if not isinstance(a, int) and not isinstance(a, float):
      a, b, c = self.estimateShape(a, b)
    assert(a >= b)
    assert(b >= c)
    self.a = a
    self.b = b
    self.c = c
    # The h and k lengths, 0 < h < k, are defined by Romain (4) or Hobson (???)
    self.h2 = a*a - b*b
    self.h = sqrt(self.h2)
    self.k2 = a*a - c*c
    self.k = sqrt(self.k2)
    return

  def estimateShape(self, points, radii):
    '''function [a,b,c,pqrdata,V] = calcInertialEllipsoid(pqrdata)'''
    N, dim = points.shape
    masses = radii**3
    M      = sum(masses)
    r0     = np.dot(masses, points) / M
    r      = np.array(points)
    for i in range(len(r)):
      r[i] = points[i] - r0

    I = np.zeros((dim, dim))
    I[0,0] = sum(masses * (r[:,1]**2 + r[:,2]**2 + (2.0/5.0) * radii**2))
    I[1,1] = sum(masses * (r[:,0]**2 + r[:,2]**2 + (2.0/5.0) * radii**2))
    I[2,2] = sum(masses * (r[:,0]**2 + r[:,1]**2 + (2.0/5.0) * radii**2))
    I[0,1] = -sum(masses * r[:,0] * r[:,1])
    I[1,0] = I[0,1]
    I[0,2] = -sum(masses * r[:,0] * r[:,2])
    I[2,0] = I[0,2]
    I[1,2] = -sum(masses * r[:,1] * r[:,2])
    I[2,1] = I[1,2]

    ev, V = np.linalg.eig(I)
    ev.sort()
    Ixx = ev[0]
    Iyy = ev[1]
    Izz = ev[2]

    g = 5.0/(2.0*M)
    a = sqrt(g * (-Ixx + Iyy + Izz))
    b = sqrt(g * ( Ixx - Iyy + Izz))
    c = sqrt(g * ( Ixx + Iyy - Izz))
    return a, b, c

  def cartesianCoords(self, l, m, n):
    '''\
Takes in a triplet of ellipsoidal coordinates (l, m, n)
Returns a triplet of Cartesian coordinates (x, y, z)
 - The coordinates satisfy k^2 \le l < \infty, h^2 \le m \le k^2, 0 \le n \le h^2'''
    h2 = self.h2
    k2 = self.k2
    # Romain (7)
    x = sqrt((l*l*m*m*n*n)/(h2*k2))
    y = sqrt(((l*l - h2)*(m*m - h2)*(h2 - n*n))/(h2*(k2 - h2)))
    z = sqrt(((l*l - k2)*(k2 - m*m)*(k2 - n*n))/(k2*(k2 - h2)))
    # These determine what quadrant we are in
    #   This is based on the Cartesian representation given in Romain Table III
    if l*m*n < 0:
      x = -x
    if l*n   < 0:
      y = -y
    if l*m   < 0:
      z = -z
    return x, y, z

  def ellipsoidalCoords(self, x, y, z):
    '''\
Takes in a triplet of Cartesian coordinates (x, y, z)
Returns a triplet of ellipsoidal coordinates (l, m, n)'''
    from math import pi, cos, acos
    h2 = self.h2
    k2 = self.k2
    # Use Romain (6) for approximate values
    a1  = -(x*x + y*y + z*z + h2 + k2)
    a2  = x*x*(h2 + k2) + y*y*k2 + z*z*h2 + h2*k2
    a3  = -x*x*h2*k2
    Q   = (a1*a1 - 3.*a2)/9.
    R   = (9.*a1*a2 - 27.*a3 - 2.*a1*a1*a1)/54.
    cth = R/(sqrt(Q*Q*Q))
    th  = acos(cth)
    #   Wikipedia http://en.wikipedia.org/wiki/List_of_trigonometric_identities
    #     cos x = 4 cos^3 x/3 - 3 cos x/3
    #     4 x^3 - 3x - v = 0  ==> {a = 4, b = 0, c = -3, d = -v}
    #     Discriminant = 18abcd - 4b^3d + b^2c^2 - 4ac^3 - 27a^2d^2
    #                  = 0 - 0 + 0 + 432 - 432*v^2 = 432*(1 - v^2) > 0 so 3 real roots
    #     r1 = -1/3a (cbrt((27a^2d + sqrt(729a^4d^2 + 108a^3c^3))/2) + cbrt((27a^2d - sqrt(729a^4d^2 + 108a^3c^3))/2))
    #        = -1/12 (cbrt((432d + 432 sqrt(d^2 - 1))/2) + cbrt((432d - 432 sqrt(d^2 - 1))/2))
    #        = -cbrt(216)/12 (cbrt((d + sqrt(v^2 - 1))) + cbrt((d - sqrt(v^2 - 1))))
    #        = cbrt(216)/(12) (cbrt(cth - i sth) + cbrt(cth + i sth))
    lambda2 = [cos(th/3), cos((th + 4*pi)/3), cos((th + 2*pi)/3)]
    l, m, n = map(sqrt, map(lambda x: 2*sqrt(Q)*x - a1/3, lambda2))
    # Get $\min(\lambda^2_i, |\lambda^2 - h^2|, |\lambda^2_i - k^2|)$
    # Rewrite ellipsoid cubic so that smallest value is root
    # Solve the equation
    # Transform back
    # Check sign from Romain Table III
    #   l m n = h k x
    #   sqrt((l2 - h2)*(m2 - h2)*(h2 - n2)) = h sqrt(k2 - h2) y
    #   sqrt((l2 - k2)*(k2 - m2)*(k2 - n2)) = k sqrt(k2 - h2) z
    #   l m n sqrt((l2 - k2)*(l2 - h2)*(k2 - m2)*(m2 - h2)*(k2 - n2)*(h2 - n2)) = h2 k2 (k2 - h2) x y z
    # We can treat only the variable signs using addition mod 2,
    #   l+m+n    = x
    #   lh+mh+nh = y
    #   lk+mk+nk = z
    # and the N comes from adding al the equations. Now
    #   l+m+n    = x
    #   l+h+m+h+n+h = y = h+l+m+n = h+x
    #   l+k+m+k+n+k = z = k+l+m+n = k+x
    # So
    #   h = sign(x y)  and k = sign(x z)
    # and now let m = h and n = k
    #   l + m + n = x
    #   l + 0 + n = y
    #   l + m + 0 = z
    # so
    #   l = sign(x y z)
    #   m = sign(x y)
    #   n = sign(x z)
    # Check
    #   l + m + n = sign(x y z x y x z) = sign(x)
    # Which means
    #   lh = l+m
    #   mh = 0
    #   nh = m+n
    #   lk = l+n
    #   mk = m+n
    #   nk = 0
    if x*y*z < 0:
      l = -l
    if x*y < 0:
      m = -m
    if x*z < 0:
      n = -n
    return l, m, n

  def getLameType(self, n, p):
    '''Gives the type of Lame function, comes from Romain Table I'''
    r = n/2
    if p < r+1:
      return 'K', p
    elif p < (n-r) + (r+1):
      return 'L', p - (r+1)
    elif p < (n-r) + (n-r) + (r+1):
      return 'M', p - (n-r) - (r+1)
    elif p < 2*n+1:
      return 'N', p - (n-r) - (n-r) - (r+1)
    raise RuntimeError('Invalid entry n,p: %d,%d' % (n, p))

  def getLameSize(self, n, t):
    '''This comes from Romain Table I'''
    if   t == 'K': return r+1
    elif t == 'L': return n-r
    elif t == 'M': return n-r
    elif t == 'N': return r
    raise RuntimeError('Invalid Lame type '+str(t))

  def getLameCoefficientMatrix(self, t, n):
    # OPT: We can memoize this later
    alpha = self.h2
    beta  = self.k2 - self.h2
    gamma = alpha - beta
    r     = n/2
    if t == 'K':
      i = np.array(range(0, r+1))
      g = (-(2*i + 2)*(2*i + 1)*beta)[:-1]
      if n%2: # n is odd
        d = ((2*r + 1)*(2*r + 2) - 4*i*i)*alpha + (2*i + 1)*(2*i + 1)*beta
        f = (-alpha*(2*(r - i) + 2)*(2*(r + i) + 1))[1:]
      else: # n is even
        d = 2*r*(2*r + 1)*alpha - 4*i*i*gamma
        f = (-alpha*(2*(r - i) + 2)*(2*(r + i) - 1))[1:]
    elif t == 'L':
      i = np.array(range(0, n-r))
      g = (-(2*i + 2)*(2*i + 3)*beta)[:-1]
      if n%2: # n is odd
        d = (2*r + 1)*(2*r + 2)*alpha - (2*i + 1)*(2*i + 1)*gamma
        f = (-alpha*(2*r - 2*i + 2)*(2*r + 2*i + 1))[1:]
      else: # n is even
        d = (2*r*(2*r + 1) - (2*i + 1)*(2*i + 1))*alpha + (2*i + 2)*(2*i + 2)*beta
        f = (-alpha*(2*r - 2*i)*(2*r + 2*i + 1))[1:]
    elif t == 'M':
      i = np.array(range(0, n-r))
      g = (-(2*i + 2)*(2*i + 1)*beta)[:-1]
      if n%2: # n is odd
        d = ((2*r + 1)*(2*r + 2) - (2*i + 1)*(2*i + 1))*alpha + 4*i*i*beta
        f = (-alpha*(2*r - 2*i + 2)*(2*r + 2*i + 1))[1:]
      else: # n is even
        d = 2*r*(2*r + 1)*alpha - (2*i + 1)*(2*i + 1)*gamma
        f = (-alpha*(2*r - 2*i)*(2*r + 2*i + 1))[1:]
    elif t == 'N':
      i = np.array(range(0, r))
      g = (-(2*i + 2)*(2*i + 3)*beta)[:-1]
      if n%2: # n is odd
        d = (2*r + 1)*(2*r + 2)*alpha - (2*i + 2)*(2*i + 2)*gamma
        f = (-alpha*(2*r - 2*i)*(2*r + 2*i + 3))[1:]
      else: # n is even
        d = 2*r*(2*r + 1)*alpha - (2*i + 2)*(2*i + 2)*alpha + (2*i + 1)*(2*i + 1)*beta
        f = (-alpha*(2*r - 2*i)*(2*r + 2*i + 1))[1:]
    M = np.diag(g, 1) + np.diag(d) + np.diag(f, -1)
    p, V = np.linalg.eig(M)
    return np.transpose(V)

  def computeLameCoefficients(self, n, p):
    '''This comes Romain Annex 3 '''
    t, tp = self.getLameType(n, p)
    B     = self.getLameCoefficientMatrix(t, n)
    return B[tp]

  def evalLame(self, n, p, l, signm = 1, signn = 1):
    '''Evaluate E^p_n(l)'''
    signh = copysign(1, signm*l)
    signk = copysign(1, signn*l)
    t, tp = self.getLameType(n, p)
    B     = self.getLameCoefficientMatrix(t, n)
    r     = n/2
    l2    = float(l*l)
    # This comes from Romain Table II
    if   t == 'K':
      m   = r
      psi = pow(l, n-2*r)
    elif t == 'L':
      m   = n-r-1
      psi = pow(l, 1-n+2*r)*signh*sqrt(abs(l2 - self.h2))
    elif t == 'M':
      m   = n-r-1
      psi = pow(l, 1-n+2*r)*signk*sqrt(abs(l2 - self.k2))
    elif t == 'N':
      m   = r-1
      psi = pow(l, n-2*r)*signh*signk*sqrt(abs((l2 - self.k2)*(l2 - self.h2)))
    else: raise RuntimeError('Invalid Lame type '+str(t))
    # Romain (30) \sum^m_{j=0} b_j (1 - \frac{\lambda^2}{h^2}) = \sum^m_{j=0} b_j \Lambda^j
    Lambda_Romain = (1.0 - (l2/self.h2)) # Romain bottom of p.252
    b = B[tp]
    # Normalize b so that the highest power of lambda has coefficient unity
    # Romain p242 and Dassios (D9,D10 and their gamma in B16-B20)
    b = b/(b[-1]/pow(-self.h2,len(b)-1))
    P = b[m]
    for j in range(m-1, -1, -1):
      P = P*Lambda_Romain + b[j]
    return psi*P

  def evalLameDerivative(self, n, p, l, signm = 1, signn = 1):
    signh = copysign(1, signm*l)
    signk = copysign(1, signn*l)
    '''Evaluate E^p_n(l)'''
    # OPT: We can return the correct function, built on the fly
    t, tp = self.getLameType(n, p)
    B     = self.getLameCoefficientMatrix(t, n)
    r     = n/2
    l2    = float(l*l)
    # This comes from Romain Table II
    if   t == 'K':
      m      = r
      psi    = pow(l, n-2*r)
      psider = (n-2*r)*pow(l, n-2*r-1)
    elif t == 'L':
      m      = n-r-1
      psi    = pow(l, 1-n+2*r)*signh*sqrt(abs(l2 - self.h2))
      psider = (1-n+2*r)*pow(l, -n+2*r)*signh*sqrt(abs(l2 - self.h2)) + pow(l, 2-n+2*r)*signh/sqrt(abs(l2 - self.h2))
    elif t == 'M':
      m      = n-r-1
      psi    = pow(l, 1-n+2*r)*signk*sqrt(abs(l2 - self.k2))
      psider = (1-n+2*r)*pow(l, -n+2*r)*signk*sqrt(abs(l2 - self.k2)) + pow(l, 2-n+2*r)*signk/sqrt(abs(l2 - self.k2))
    elif t == 'N':
      m      = r-1
      psi    = pow(l, n-2*r)*signh*signk*sqrt(abs((l2 - self.k2)*(l2 - self.h2)))
      psider = ((n-2*r)*pow(l, n-2*r-1)*signh*signk*sqrt(abs((l2 - self.k2)*(l2 - self.h2))) +
                pow(l, n-2*r+1)*signh*signk*sqrt(abs((l2 - self.k2)/(l2 - self.h2))) +
                pow(l, n-2*r+1)*signh*signk*sqrt(abs((l2 - self.h2)/(l2 - self.k2))))
    else: raise RuntimeError('Invalid Lame type '+str(t))
    # Romain (30) \sum^m_{j=0} b_j (1 - \frac{\lambda^2}{h^2}) = \sum^m_{j=0} b_j \Lambda^j
    Lambda_Romain = (1.0 - (l2/self.h2)) # Romain bottom of p.252
    b    = B[tp]
    # Normalize b so that the highest power of lambda has coefficient unity
    # Romain p242 and Dassios (D9,D10 and their gamma in B16-B20)
    b    = b/(b[-1]/pow(-self.h2,len(b)-1))
    P    = b[m]
    for j in range(m-1, -1, -1):
      P    = P*Lambda_Romain + b[j]
    # P' = \sum^{m-1}_{k=0} c_k \Lambda^k where c_k = b_{k+1} (-2 (k+1) \frac{\lambda}{h^2})
    Pder = b[m]*(-2*m*float(l)/self.h2)
    for k in range(m-2, -1, -1):
      Pder = Pder*Lambda_Romain + b[k+1]*(-2*(k+1)*l/self.h2)
    return psider*P + psi*Pder

  def calcI(self, n, p, l, signm = 1, signn = 1):
    '''I^p_n(\lambda) = \int^{1\lambda}_0 \frac{t^{2n} dt}{(E^p_n(1/t))^2 \sqrt{1 - k^2 t^2}\sqrt{1 - h^2 t^2}}'''
    l  = float(l)
    k2 = self.k2
    h2 = self.h2
    N  = 10000
    def integrand(t):
      t2 = t*t
      s  = 1e6 if t == 0.0 else 1.0/t
      E  = self.evalLame(n, p, s, signm, signn)
      return 1.0 / (E*E*sqrt(1.0 - k2*t2)*sqrt(1.0 - h2*t2))
    return integrate(integrand, 0.0, 1.0/l, N)

  def calcIDerivative(self, n, p, l, signm = 1, signn = 1):
    '''\dot I^p_n(\lambda) = \frac{-1}{(E^p_n(\lambda))^2 \sqrt{\lambda^2 - k^2}\sqrt{\lambda^2 - h^2}}'''
    l = float(l)
    l2 = l*l
    k2 = self.k2
    h2 = self.h2
    E  = self.evalLame(n, p, l, signm, signn)
    return -1.0 / (E*E*sqrt(l2 - k2)*sqrt(l2 - h2))

  def calcEigenvalue(self, n, p, l, signm = 1, signn = 1):
    '''Calculate eigenvalues of the surface operator'''
    a, b, c = self.a, self.b, self.c
    I    = self.calcI(n, p, l, signm, signn)
    E    = self.evalLame(n, p, l, signm, signn)
    Eder = self.evalLameDerivative(n, p, l, signm, signn)
    ev   = (2*a*b*c*(Eder/a)*I*E - 1)/2
    return ev

  def calcNormalization(self, n, p):
    '''mathcalI1  0.6225, 3.9913, 0.6992, 1.9596
       Imn       -1.3919, 0.3931, 1.5634, 0.6870
       alpha     -0.0973, 1.2390'''
    h  = float(self.h)
    h2 = float(self.h2)
    k  = float(self.k)
    k2 = float(self.k2)
    N  = 10000
    # The four integrals making up Romain (51)
    I1b = romberg(lambda l: self.evalLame(n, p, l)**2        / sqrt((l**2 - h2)*(k2 - l**2)), h+1e-5, k-1e-5, tol=1.48e-8, rtol=1.48e-8, divmax=10)
    I1 = integrateMidpoint(lambda l: self.evalLame(n, p, l)**2        / sqrt((l**2 - h2)*(k2 - l**2)), h, k, N)
    I2 = integrateMidpoint(lambda l: self.evalLame(n, p, l)**2 * l**2 / sqrt((l**2 - h2)*(k2 - l**2)), h, k, N)
    I3 = integrateMidpoint(lambda l: self.evalLame(n, p, l)**2        / sqrt((h2 - l**2)*(k2 - l**2)), 0, h, N)
    I4 = integrateMidpoint(lambda l: self.evalLame(n, p, l)**2 * l**2 / sqrt((h2 - l**2)*(k2 - l**2)), 0, h, N)
    # The basic elliptic integrals Romain (53)
    #   He has an error in the equation, the prefactor should be h/2
    I20 = 0.5*h*integrateMidpoint(lambda Lambda: 1.0    / (sqrt(1 - Lambda)*sqrt(h2*Lambda + (k2 - h2))*sqrt(-Lambda)), 0.0, 1.0 - k2/h2, N)
    I21 = 0.5*h*integrateMidpoint(lambda Lambda: Lambda / (sqrt(1 - Lambda)*sqrt(h2*Lambda + (k2 - h2))*sqrt(-Lambda)), 0.0, 1.0 - k2/h2, N)
    I30 = 0.5*h*integrateMidpoint(lambda Lambda: 1.0    / (sqrt(1 - Lambda)*sqrt(h2*Lambda + (k2 - h2))*sqrt( Lambda)), 0.0, 1.0, N)
    I31 = 0.5*h*integrateMidpoint(lambda Lambda: Lambda / (sqrt(1 - Lambda)*sqrt(h2*Lambda + (k2 - h2))*sqrt( Lambda)), 0.0, 1.0, N)
    # Solve system for coefficients
    matrix = np.matrix([[I20, I21], [I30, I31]])
    alpha, beta = np.linalg.solve(matrix, [I1, I3])
    A,     B    = np.linalg.solve(matrix, [I2, I4])
    # Romain (54)
    #   The factor 8 is from Dassios (see discussion in Deng), since we integrate over the whole ellipsoid rather than just one octant
    gamma = 8 * (pi/2.0)*(alpha*B - beta*A)
    return gamma

class BEM(object):
  def __init__(self, filename, sourceFilename = None):
    data = np.loadtxt(filename, unpack = False)
    self.dim       = (data.shape[1]-1)/2
    self.N         = data.shape[0]
    self.centroids = np.array(data[:,0:self.dim])
    assert(self.centroids.shape == (self.N, self.dim))
    self.normals   = np.array(data[:,self.dim:self.dim*2])
    assert(self.normals.shape == (self.N, self.dim))
    self.areas     = np.array(data[:,self.dim*2])
    assert(self.areas.shape == (self.N,))
    if sourceFilename:
      self.readSources(sourceFilename)
      self.e       = EllipsoidalSystem(self.sourcePoints, self.sourceRadii)
    return

  def readSources(self, filename):
    sourcePoints  = []
    sourceCharges = []
    sourceRadii   = []
    for line in file(filename):
      if len(line) > 6 and (line.startswith('ATOM') or line.startswith('HETATM')):
        # Field name, Atom number, Atom name, Residue name, Chain ID, Residue number, Source point, Source charge, Source radius
        fieldName,atomNum,atomName,resName,chainID,resNum,x,y,z,charge,radius = line.split()
        sourcePoints.append(map(float, [x,y,z]))
        sourceCharges.append(float(charge))
        sourceRadii.append(float(radius))
    self.sourcePoints  = np.array(sourcePoints)
    self.sourceCharges = np.array(sourceCharges)
    self.sourceRadii   = np.array(sourceRadii)
    self.Nq = len(self.sourceRadii)
    return

  def calcOperators(self):
    '''Generate the single-layer (V) and double-layer (K) operators'''
    V = np.zeros((self.N, self.N))
    K = np.zeros((self.N, self.N))
    for i in range(self.N):
      for j in range(self.N):
        if i == j: continue
        vr = self.centroids[i] - self.centroids[j]
        r  = np.linalg.norm(vr)
        V[i,j] = self.areas[j] / (4*pi*r)
        K[i,j] = self.areas[j] * np.dot(vr, self.normals[j]) / (4*pi*r*r*r)
    return np.matrix(V), np.matrix(K)

  def calcSourceOperators(self):
    '''Generate the map from sources to normal electric field at the surface (C), and map from surface charge to potential (B)'''
    B = np.zeros((self.N, self.Nq))
    C = np.zeros((self.N, self.Nq))
    for i in range(self.N):
      for j in range(self.Nq):
	 vr   = self.centroids[i] - self.sourcePoints[j]
	 r    = np.linalg.norm(vr)
         G    = 0
         dGdn = 0
	 if r >= 1e-10:
           G    = 1.0/(4*pi*r)
           dGdn = -np.dot(vr, self.normals[i])/(4*pi*r**3)
	 C[i,j] = G
	 B[i,j] = dGdn
    return np.matrix(B), np.matrix(C)

import unittest

class TestLame(unittest.TestCase):
  def setUp(self):
    pass

  def getLameFormula(self, e, n):
    formula = ''
    for p in range(2*n+1):
      t, tp = e.getLameType(n, p)
      formula += t+'^%d_%d ' % (tp, n)
    return formula.strip()

  def outputLameTypes(self, nMax):
    e = EllipsoidalSystem(3, 2, 1)
    for n in range(nMax):
      print 'n = '+str(n)+':',getLameFormula(e, n)
    return

  def testLameType(self):
    formulas = {1: 'K^0_1 L^0_1 M^0_1',
                2: 'K^0_2 K^1_2 L^0_2 M^0_2 N^0_2',
                5: 'K^0_5 K^1_5 K^2_5 L^0_5 L^1_5 L^2_5 M^0_5 M^1_5 M^2_5 N^0_5 N^1_5'}
    e = EllipsoidalSystem(3, 2, 1)
    for n, f in formulas.iteritems():
      self.assertEqual(f, self.getLameFormula(e, n))
    return

  def testCoordinateTransform(self):
    e   = EllipsoidalSystem(3, 2, 1)
    eps = 1.0e-10
    # Loop over octants
    for sx in [1, -1]:
      for sy in [1, -1]:
        for sz in [1, -1]:
          # Loop over brick in cartesian space
          for x in np.arange(4, 5+eps, 1.5):
            for y in np.arange(4, 5+eps, 1.5):
              for z in np.arange(4, 5+eps, 1.5):
                l, m, n = e.ellipsoidalCoords(sx*x, sy*y, sz*z)
                try:
                  nx, ny, nz = e.cartesianCoords(l, m, n)
                  relError = [abs(sx*x - nx)/x, abs(sy*y - ny)/y, abs(sz*z - nz)/z]
                  self.assertAlmostEqual(sx, nx/x, 10)
                  self.assertAlmostEqual(sy, ny/y, 10)
                  self.assertAlmostEqual(sz, nz/z, 10)
                except ValueError, ex:
                  print 'Invalid ellipsoidal coordinates',l,m,n
    return

  def testI(self):
    e = EllipsoidalSystem(3, 2, 1)
    a = e.a
    b = e.b
    c = e.c
    h = e.h
    k = e.k
    l = a
    rootSum            = sqrt((a**4 - b**2*c**2) + (b**4 - a**2*c**2) + (c**4 - a**2*b**2))
    DassiosLambda      = (1/3.)*(a**2+b**2+c**2) + (1/3.)*rootSum
    DassiosLambdaPrime = (1/3.)*(a**2+b**2+c**2) - (1/3.)*rootSum
    # Dassios D1 error
    lhs = 3 * (DassiosLambda + DassiosLambdaPrime)
    rhs = 2 * (a**2 + b**2 + c**2)
    self.assertAlmostEqual(lhs, rhs)
    # Dassios D2 error
    lhs = 3 * DassiosLambda * DassiosLambdaPrime
    rhs = a**2*b**2 + a**2*c**2 + b**2*c**2;
    self.assertAlmostEqual(lhs, rhs)
    # Dassios D3 error
    lhs = ((-1)*(b**2-c**2)*(DassiosLambda-a**2)) + ((-1)**2*k**2*(DassiosLambda-b**2)) + ((-1)**3*h**2*(DassiosLambda-c**2))
    rhs =  ((-1)*(b**2-c**2)*(DassiosLambdaPrime-a**2)) + ((-1)**2*k**2*(DassiosLambdaPrime-b**2)) + ((-1)**3*h**2*(DassiosLambdaPrime-c**2))
    self.assertAlmostEqual(lhs, 0.0)
    self.assertAlmostEqual(rhs, 0.0)
    # Dassios D4 error
    lhs = ((-1)**1*a**2*(b**2-c**2)*(DassiosLambda-a**2)) + ((-1)**2*b**2*k**2*(DassiosLambda-b**2)) + ((-1)**3*c**2*h**2*(DassiosLambda-c**2))
    rhs = ((-1)**1*a**2*(b**2-c**2)*(DassiosLambdaPrime-a**2)) + ((-1)**2*b**2*k**2*(DassiosLambdaPrime-b**2)) + ((-1)**3*c**2*h**2*(DassiosLambdaPrime-c**2))
    analytical = h**2*k**2*(b**2-c**2);
    self.assertAlmostEqual(lhs, analytical)
    self.assertAlmostEqual(rhs, analytical)
    # Dassios D5 error
    alpha = np.array([a, b, c])
    lhs = sum(alpha*alpha / (alpha*alpha - DassiosLambda))
    rhs = sum(alpha*alpha / (alpha*alpha - DassiosLambdaPrime))
    self.assertAlmostEqual(lhs, 3)
    self.assertAlmostEqual(rhs, 3)
    # Dassios D6 error
    anal1 = (-1)**(1+1)*h**2*k**2*(b**2-c**2);
    num1  = 3*(b**2-c**2)*(DassiosLambda - a**2)*(DassiosLambdaPrime-a**2);
    anal2 = (-1)**(2+1)*h**2*k**2*(b**2-c**2);
    num2  = 3*k**2*(DassiosLambda - b**2)*(DassiosLambdaPrime-b**2);
    anal3 = (-1)**(3+1)*h**2*k**2*(b**2-c**2);
    num3  = 3*h**2*(DassiosLambda - c**2)*(DassiosLambdaPrime-c**2);
    self.assertAlmostEqual(num1, anal1)
    self.assertAlmostEqual(num2, anal2)
    self.assertAlmostEqual(num3, anal3)
    # Calculate I^p_n
    nmax = 2
    I    = []
    for n in range(nmax+1):
      for p in range(2*n+1):
        I.append(e.calcI(n, p, l))
    # I = [0.50849469942291659, 0.02603346073135392, 0.044488203423648773, 0.095941114824094301, 0.0029857137346618505, 0.014964156240992845, 0.00369094853845898, 0.0087384567615925567, 0.017150970466815185]
    # Dassios D7: their h_3 = our h and their h_2 = our k 
    analyticalSumTest1 = 1.0/(l*sqrt(l**2-h**2)*sqrt(l**2-k**2))
    numericalSumTest1  = sum(I[1:4])
    self.assertAlmostEqual(numericalSumTest1, analyticalSumTest1, 3)
    # Dassios D8: their alpha1 = our a, their alpha2 = our b, their alpha3 = our c
    analyticalSumTest2 = I[0] - (l**2-a**2)/(l*sqrt(l**2-h**2)*sqrt(l**2-k**2));
    numericalSumTest2 = a**2*I[1] + b**2*I[2] + c**2*I[3]
    self.assertAlmostEqual(numericalSumTest2, analyticalSumTest2, 3)
    # Dassios D9: their Lambda = our DassiosLambda, their LambdaPrime = our DassiosLambdaPrime
    analyticalSumTest3 = 1/(2*(DassiosLambda-a**2+l**2)*l*sqrt(l**2-h**2)*sqrt(l**2-k**2)) - (1.0/2.0)*(I[1]/(DassiosLambda-a**2) + I[2]/(DassiosLambda-b**2) + I[3]/(DassiosLambda-c**2))
    numericalSumTest3 = I[4]
    self.assertAlmostEqual(numericalSumTest3, analyticalSumTest3, 3)
    # Dassios D10: their Lambda = our DassiosLambda, their LambdaPrime = our DassiosLambdaPrime
    analyticalSumTest4 = 1/(2*(DassiosLambdaPrime-a**2+l**2)*l*sqrt(l**2-h**2)*sqrt(l**2-k**2)) - (1.0/2.0)*(I[1]/(DassiosLambdaPrime-a**2) + I[2]/(DassiosLambdaPrime-b**2) + I[3]/(DassiosLambdaPrime-c**2))
    numericalSumTest4 = I[5]
    self.assertAlmostEqual(numericalSumTest4, analyticalSumTest4, 3)
    # Dassios D11: their h1**2 = our (b**2-c**2)
    analyticalSumTest5 = (1.0/h**2) * (I[2]-I[1])
    numericalSumTest5  = I[6]
    self.assertAlmostEqual(numericalSumTest5, analyticalSumTest5)
    # Dassios D12: their h1**2 = our (b**2-c**2)
    analyticalSumTest6 = (1.0/k**2) * (I[3]-I[1])
    numericalSumTest6  = I[7]
    self.assertAlmostEqual(numericalSumTest6, analyticalSumTest6)
    # Dassios D13: their h1**2 = our (b**2-c**2)
    analyticalSumTest7 = (1.0/(b**2-c**2)) * (I[3]-I[2])
    numericalSumTest7  = I[8]
    self.assertAlmostEqual(numericalSumTest7, analyticalSumTest7)
    return

  def testEigenvalues(self):
    '''This tests the eigenvalues of the surface operator'''
    e = EllipsoidalSystem(3, 2, 1)
    a,  b,  c  = e.a, e.b, e.c
    a2, b2, c2 = a*a, b*b, c*c
    l = a
    n = 1
    # Get analytic solutions
    analytic = []
    N = 50000
    L = 1.0e3
    # Integrals are from Ritter, need to transform integrals to finite interval
    integral = integrate(lambda s: 1.0 / ((s + a2)*sqrt((s + a2)*(s + b2)*(s + c2))), 0, L, N)
    analytic.append((a*b*c * integral - 1.0)/2)
    integral = integrate(lambda s: 1.0 / ((s + b2)*sqrt((s + a2)*(s + b2)*(s + c2))), 0, L, N)
    analytic.append((a*b*c * integral - 1.0)/2)
    integral = integrate(lambda s: 1.0 / ((s + c2)*sqrt((s + a2)*(s + b2)*(s + c2))), 0, L, N)
    analytic.append((a*b*c * integral - 1.0)/2)
    for p in range(2*n+1):
      ev = e.calcEigenvalue(n, p, l)
      self.assertAlmostEqual(ev, analytic[p], 3)
    return

  def testNormalization(self):
    import random
    #c, b, a = sorted([random.uniform(1, 7), random.uniform(1, 5), random.uniform(1, 3)])
    a = 3
    b = 2
    c = 1
    print a, b, c
    e = EllipsoidalSystem(a, b, c)
    # Dassios (B14)
    firstTerm    = (a**2 + b**2 + c**2)/3.0
    secondTerm   = sqrt((a**4 - b**2*c**2) + (b**4 - a**2*c**2) + (c**4 - a**2*b**2))/3.0
    LambdaD      = firstTerm + secondTerm
    LambdaDprime = firstTerm - secondTerm
    # Dassios (B16)-(B20)
    hx = sqrt(b**2 - c**2)
    hy = e.k
    hz = e.h
    analytic = [[4*pi],
                [4*pi/3 * hy**2*hz**2, 4*pi/3 * hx**2*hz**2, 4*pi/3 * hx**2*hy**2],
                [-8*pi/5 * (LambdaD - LambdaDprime)*(LambdaD      - a**2)*(LambdaD      - b**2)*(LambdaD      - c**2),
                  8*pi/5 * (LambdaD - LambdaDprime)*(LambdaDprime - a**2)*(LambdaDprime - b**2)*(LambdaDprime - c**2),
                  4*pi/15 *hx**2*hy**2*hz**4,
                  4*pi/15 *hx**2*hy**4*hz**2,
                  4*pi/15 *hx**4*hy**2*hz**2]]
    for n in range(0, 3):
      for p in range(2*n+1):
        self.assertAlmostEqual(e.calcNormalization(n, p), analytic[n][p], 3)
    return

  def testChargeInEllipsoid(self):
    '''Checking accuracy of approximation of surface operator by expansion in ellipsoidal harmonics:
    - Put point charges in a dielectric ellipsoid
    - Compute B_{np} expansion coefficients for the reaction potential
    - Calculate both point BEM and approximate energy and eigenspaces
    '''
    nmax = 3
    eps1 = 4.
    eps2 = 80.
    e = EllipsoidalSystem(15, 12, 10)
    charges = []
    for r, q in [((3, 4, 5), 1)]:
      xi, mu, nu = e.ellipsoidalCoords(r[0], r[1], r[2])
      print xi, mu, nu
      signm = copysign(1, mu)
      signn = copysign(1, nu)
      charges.append(((xi, mu, nu, signm, signn), q))
    print 'charges', charges
    # Compute E_{np}
    E = np.zeros((nmax, 2*nmax+1))
    for n in range(0, nmax):
      for p in range(0, 2*n+1):
        print 'n',n,'p',p
        for (xi, mu, nu, signm, signn), q in charges:
          factor = q*(4*pi/(2*n+1))
          print '  factor',factor
          lame   = e.evalLame(n, p, xi, signm, signn)*e.evalLame(n, p, mu, signm, signn)*e.evalLame(n, p, nu, signm, signn)
          print '  lame',lame
          gamma  = e.calcNormalization(n, p)
          print '  gamma',gamma
          E[n,p] += factor*lame/gamma
          print '  E',E[n,p]
    # Compute B_{np}
    B = np.zeros((nmax, 2*nmax+1))
    for n in range(0, nmax):
      for p in range(0, 2*n+1):
        print 'n',n,'p',p
        factor   = (eps1 - eps2)/(eps1*eps2)
        print '  factor',factor
        lameI    = e.calcI(n, p, e.a)
        ratio    = (2*n+1)*lameI
        print '  ratio',ratio
        # \tilde F = \dot F / F = (\dot E I + E \dot I) / (I E)
        # \tilde E = \dot E / E
        # \tilde E / \tilde F = (\dot E I) / (\dot E I + E \dot I)
        lameE    = e.evalLame(n, p, e.a)
        lameEDer = e.evalLameDerivative(n, p, e.a)
        lameIDer = e.calcIDerivative(n, p, e.a)
        print '  lameE',lameE,'lameEDer',lameEDer,'lameIDer',lameIDer
        ratio2   = (lameEDer * lameI)/(lameEDer * lameI + lameE * lameIDer)
        print '  ratio2',ratio2
        B[n,p] = factor*ratio*E[n,p]/(1 - (eps1/eps2)*ratio2)
        print '  B',B[n,p]
    return

class TestBEM(unittest.TestCase):
  def setUp(self):
    pass

  def calcSphereEnergy(self, filename, sourceFilename):
    bem    = BEM(filename, sourceFilename)
    epsIn  = 4.0
    epsOut = 80.0
    Area   = np.diag(bem.areas)
  
    V,  K  = bem.calcOperators()
    Bp, Cp = bem.calcSourceOperators()
    B = -(4*pi/epsIn) * Area * Bp
    C =  (4*pi/epsIn) * Cp.T * Area
    A =  (4*pi/epsIn) * K.T * Area + (2*pi*((epsIn+epsOut)/(epsIn*(epsIn-epsOut)))) * Area
    calculatedEnergy = (0.5 * 332.112) * C * np.linalg.solve(A, B)
    # Get analytic info
    a = float(file(filename).readlines()[0][1:].split()[1])
    analyticEnergy = 0.5 * 332.112 * (1/epsOut - 1/epsIn) * (1/a)
    return calculatedEnergy, analyticEnergy

  def testSphere(self):
    for npanels in [100, 250, 400]:
    #for npanels in [100]:
      calcE, analE = self.calcSphereEnergy('../geometry/ellipsoid_points/sphere_mesh_rad3_%d.dat' % npanels, '../geometry/ellipsoid_points/monopole.pqr')
      # Absolute error (kcal/mol)
      self.assertTrue(abs(calcE[0,0] - analE)/analE < 0.06)
    return

if __name__ == '__main__':
  if 1:
    unittest.main()
  else:
    e = EllipsoidalSystem(15, 12, 10)
    xi, mu, nu = e.ellipsoidalCoords(1.0, 0.0, 0.0)
    e.evalLame(1, 2, xi)
    e.evalLame(1, 2, mu)
    e.evalLame(1, 2, nu)
    sys.exit(0)
    bem = BEM('../geometry/ellipsoid_points/mesh_test1_res8.dat')
    print bem.e.a, bem.e.b, bem.e.c
    print bem.calcOperators()
    sys.exit(0)
    e = EllipsoidalSystem(3, 2.95, 2.9)
    k = sqrt(3*3 - 1*1)
    h = sqrt(3*3 - 2*2)
    for n in range(2, 5):
      for p in range(2*n+1):
        print 'n',n,'p',p,e.calcEigenvalue(n, p, e.a)
    
