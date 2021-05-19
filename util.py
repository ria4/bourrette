#!/usr/bin/python

"""
Find a, b, c so that f(x1) = y1, f(x2) = y2, f(x3) = y3
where f(x) = exp(a*x**2 + b*x + c).
"""

import math

"""
# fibers length conversion
x1 = 1;     y1 = 1/1000.0
x2 = 50;    y2 = 1/10.0
x3 = 100;   y3 = 2.0
"""

# fibers width conversion
x1 = 1;     y1 = 1/1000.0
x2 = 50;    y2 = 1/50.0
x3 = 100;   y3 = 1/2.0

z1 = math.log(y1)
z2 = math.log(y2)
z3 = math.log(y3)

a = (x1*(z3-z2)+x2*(z1-z3)+x3*(z2-z1))/((x1-x2)*(x1-x3)*(x2-x3))
b = (z2-z1)/(x2-x1)-a*(x1+x2)
c = z1-a*x1**2-b*x1

def f(x):
    return math.exp(a*x**2 + b*x + c)

print(a)
print(b)
print(c)
print(abs(y1-f(x1)), abs(y2-f(x2)), abs(y3-f(x3)))
print(-b/2.0/a)


"""
a = ((y2-y1)*(x3-x2)-(y3-y2)*(x2-x1)) / ((x2**3-x1**3)*(x3-x2)-(x3**3-x2**3)*(x2-x1))
c = (y2-y1-a*(x2**3-x1**3)) / (x2-x1)
d = y1 - a*x1**3 - c*x1
"""

