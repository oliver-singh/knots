import sympy


def _homology_generators(braid):
    hom_generators = [0]*(len(braid)-1)
    for i in range(len(braid)-1):
        hom_generators[i] = 0
        j = i+1
        while hom_generators[i] == 0 and j < len(braid):
            if abs(braid[i]) == abs(braid[j]):
                hom_generators[i] = j
            j = j + 1
    return hom_generators


def seifert_matrix(braid):
    hom_generators = _homology_generators(braid)
    s_matrix = sympy.zeros(len(hom_generators), len(hom_generators))

    for i in range(len(hom_generators)):
        if hom_generators[i] == 0:
            continue

        if braid[i] > 0 and braid[hom_generators[i]] > 0:
            s_matrix[i, i] = -1
        elif braid[i] < 0 and braid[hom_generators[i]] < 0:
            s_matrix[i, i] = 1

        for j in range(i+1, len(hom_generators)):
            if hom_generators[i] == j:
                if braid[j] > 0:
                    s_matrix[j, i] = 1
                else:
                    s_matrix[i, j] = -1

            if i < j < hom_generators[i] < hom_generators[j]:
                separating_strands = abs(braid[i]) - abs(braid[j])
                if separating_strands == 1:
                    s_matrix[j, i] = -1
                elif separating_strands == -1:
                    s_matrix[i, j] = 1

    for i in range(len(hom_generators) - 1, 0, -1):
        if hom_generators[i] == 0:
            s_matrix.row_del(i)
            s_matrix.col_del(i)
    return s_matrix

def _get_reduced_burau_matrices(n):
    """
    Generates reduced Burau matrices for braid group B_n
    :param n: number of strands of braid group
    :return: dictionary "reps" of reduced Burau matrices with reps[i] the Burau matrices for the ith braid generator
    """
    t = sympy.Symbol('t')
    assert n >= 2
    reps = {}
    if n == 2:
        reps[1], reps[-1] = sympy.eye(1), sympy.eye(1),
        reps[1][0, 0], reps[-1][0, 0]  = -t, -1/t
    else:
        reps[1], reps[-1] = sympy.eye(n - 1), sympy.eye(n - 1)
        reps[1][0, 0], reps[-1][0, 0] = -t, -1/t
        reps[1][0, 1], reps[-1][0, 1] = 1, 1/t
        for i in range(2, n - 1):
            reps[i], reps[-i] = sympy.eye(n - 1), sympy.eye(n - 1)
            reps[i][i - 1, i - 2:i + 1] = [[t, -t, 1]]
            reps[-i][i - 1, i - 2:i + 1] = [[1, -1/t, 1/t]]
        reps[n - 1], reps[-n + 1] = sympy.eye(n - 1), sympy.eye(n - 1)
        reps[n - 1][n - 2, n - 3:n - 1] = [[t, -t]]
        reps[-n + 1][n - 2, n - 3:n - 1] = [[1, -1/t]]
    return reps

def burau_rep(braid):
    t = sympy.Symbol('t')
    reps = _get_reduced_burau_matrices(braid.n_strands)
    matrix = sympy.eye(braid.n_strands - 1)
    for gen in braid:
        matrix = matrix * reps[gen]
    return matrix

def burau_to_alexander(matrix):
    a, b = matrix.shape
    n = a + 1
    t = sympy.Symbol('t')
    matrix = sympy.eye(a) - matrix
    alex_poly = matrix.det()
    alex_poly = alex_poly * (1 - t) * t**(2*n)
    alex_poly = sympy.pexquo(alex_poly, ( 1 - t**n))
    alex_poly = sympy.PurePoly(alex_poly.expand())
    while alex_poly.nth(0) == 0:
        alex_poly = sympy.exquo(alex_poly, sympy.PurePoly(t))
    if alex_poly.nth(0) < 0:
        alex_poly = - alex_poly
    return alex_poly

def seifert_to_alexander(seifert_matrix, method = "seifert"):
    if method == "seifert":
        t = sympy.Symbol('t')
        n, m = seifert_matrix.shape
        assert n == m, "non square matrix received"
        matrix = sympy.SparseMatrix(t * seifert_matrix - seifert_matrix.transpose())
        alex_poly = t**(-n//2) * matrix.det(method = "berkowitz")
        alex_poly = sympy.PurePoly(sympy.expand(alex_poly))
    if method == "bureau":
        pass
    return alex_poly


def signature(seifert_matrix):
    n, m = seifert_matrix.shape
    assert n == m, "non square matrix received"
    matrix = seifert_matrix + seifert_matrix.transpose()

    t = sympy.Symbol('t')
    characteristic_poly = matrix.charpoly(t)
    #remove 0 eigenvalues
    while characteristic_poly.nth(0) == 0:
        characteristic_poly = sympy.exquo(characteristic_poly, sympy.PurePoly(t))

    sig = sympy.degree(characteristic_poly)
    negative_root_intervals = characteristic_poly.intervals(sup=0)
    for interval, multiplicity in negative_root_intervals:
        assert interval[1] <= 0, ("Error in signature computation, root of characteristic equation was not isolated" +
                                  " correctly below zero. Could be a sympy version problem (working in 1.7.1).")
        sig += -2 * multiplicity
    return sig

