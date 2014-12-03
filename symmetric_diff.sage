import sys
from sage.all import *
from collections import defaultdict


R.<z1, z2, dz1, dz2> = PolynomialRing(QQ, ['z1', 'z2', 'dz1', 'dz2'])

#INPUT: List of variables, var_list, and degree of symmetric differentials, degree.
#OUTPUT: List of symmetric differentials in variables dictated by var_list that
# 1.) Are invariant under Z2 action 
# 2.) Have degrees (i1, i2, m1, m2) where i1 + i2 < m1 + m2.

def invariant_monomials_less_than(var_list, degree):
    differential_list = monomials(var_list[2: 4], [degree+1 for var in range(2)])
    differential_list = [diff for diff in differential_list 
                              if diff.total_degree() == degree]
      
    z1, z2, dz1, dz2 = var_list
    inv_monomials = []
    for differential in differential_list:
        for i in range(degree):
            if i % 2 == degree % 2:
                for j in range(i + 1):
                    if differential.degree(dz2) - (i-j) <= degree:
                        inv_monomials += [z1^(i-j)*z2^j*differential]
    return inv_monomials

#INPUT: Polynomial ring, ring, and list of monomials.
#OUTPUT: Image of monomials after pullback by 
# f(u, v, du, dv) = (u, uv, du, udv + vdu)

def fstar(ring, monomial_list):
    f = ring.hom([u, u*v, du, u*dv + v*du])
    return [f(monomial) for monomial in monomial_list]

#INPUT: Polynomial ring, ring, and list of polynomials
#OUTPUT: Image of polynomials after modding out by equivalence by:
# g(u, v, du, dv) = (u^2, v, 2du, dv)
def mod_gstar(ring, polynomial_list):

    basis = defaultdict(lambda: False)
    count = 1
    n = len(polynomial_list)
    coeff_matrix = copy(matrix(ZZ, n, n))
    for index, polynomial in enumerate(polynomial_list):
        for monomial in map(mul, zip(polynomial.coefficients(), polynomial.monomials())):
            if monomial.degree(u) < monomial.degree(du):
                if not basis[monomial.lm()]:
                    basis[monomial.lm()] = count
                    count += 1
                coeff_matrix[index, basis[monomial.lm()] - 1] = monomial.lc()
    print coeff_matrix
    return coeff_matrix.rank()



S.<u, v, du, dv> = PolynomialRing(QQ, ['u', 'v', 'du', 'dv'])    

def generate_sums(i, m):
    sum1, sum2, sum3 = 0, 0, 0
    for c in range(i + 1):
        sum1 += min(c + 1, (m - i) / 2)
    for c in range(i + 1, m - 1 + 1):
        sum2 += min(i + 1, (m - i) /2)
    for c in range(m, m + i + 1):
        sum3 += min(m - c + i + 1, (m - i)/2)
    return sum1 + sum2 + sum3

def test_sums(lower_bound, upper_bound):
    for m in range(lower_bound, upper_bound):
        total_sum = 0
        i = m - 2
        while i >= 0:
            total_sum += generate_sums(i, m)
            i -= 2   

        invariants = invariant_monomials_less_than([z1, z2, dz1, dz2], m)
        fstar_inv = fstar(S, invariants)
        inv_mod_gstar = mod_gstar(S, fstar_inv)

        print total_sum, inv_mod_gstar
        assert(inv_mod_gstar == total_sum)

test_sums(22, 50)




