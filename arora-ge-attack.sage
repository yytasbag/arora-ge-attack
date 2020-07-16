# from sage.all import *

q = 127
F = GF(q)
n = 16
B = 1 # Errors are sampled from [-B, B]
s = vector(F, [F.random_element() for _ in range(n)])

# this is an implementation of lwe
def lwe_oracle():
    e = randint(-B, B)
    a = vector(F, [F.random_element() for _ in range(n)])
    return (a, a.dot_product(s) + e)


R = PolynomialRing(F, names=tuple([f"x{i}" for i in range(n)]))
R.inject_variables() # defining variables
xs = [eval(f"x{i}") for i in range(n)] # Horrible, horrible hack :(

def gen_poly(a, c):
    assert len(a) == n
    f = 1
    for e in range(-B, B + 1):
        g = 0
        for i in range(n):
            g += a[i] * xs[i]
        g = g + e - c
        f *= g
    return f # Notice that s is a root of f

# test that s is a root
# a, c = lwe_oracle()
# assert gen_poly(a, c)(*s) == 0 

# generating monomials
mons = []
for b in range(1, 2*B + 1 + 1): # dont include 1
    temp = list(WeightedIntegerVectors(b, [1] * 16))
    for dist in temp:
        vs = list(map(lambda t: t[1]^dist[t[0]], enumerate(xs[::-1])))
        mons += [mul(vs)] 
no_of_monomials = len(mons)

A = []
V = []

# this is where the linearization occur
for i in range(no_of_monomials):
    a, c = lwe_oracle()
    f = gen_poly(a, c)
    A.append([f.monomial_coefficient(m) for m in mons])
    V.append(-f.constant_coefficient())
    

A = Matrix(F, A)
V = vector(F, V)
T = A.solve_right(V)

print(s)
print(T[:n])
assert s == T[:n]