from sympy import Matrix
from sympy.abc import t

def all_factors(n):
    i = 2
    factor_set = {1}
    while i ** 2 <= n:
        if n % i == 0:
            factor_set = factor_set.union({factor * i for factor in factor_set})
            n = n // i
        else:
            i +=1
    factor_set = factor_set.union({factor * n for factor in factor_set})
    return factor_set

def factor_pairs(n):
    factors = all_factors(n)
    for factor in factors:
        if factor ** 2 < n:
            yield (factor, n // factor)

for b in range(5):
    for a,c in factor_pairs(5):
        S_matrix = Matrix([[a,b],[b+1,c]])
        int_form = t * S_matrix - S_matrix.transpose()

class SeifertMatrix(Matrix):
    def int_form(self, eval=None):
        if eval is None:
            eval = t
        form = t * self - self.transpose()
        return form

    def __repr__(self):
        return "Seifert" + super().__repr__()
    def __str__(self):
        return "Seifert" + super().__str__()

    def has_non_standard_intersection(self, tests):
        non_standard = False
        for test in tests:
            if not are_congruent(A.int_form(test), b.int_form(test)):
                non_congruent = True
                break
        return non_congruent

standard = SeifertMatrix([[0, 1], [0, 0]])

def congruent_to_standard(A, t):




