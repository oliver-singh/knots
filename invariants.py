import sympy
from sympy import PurePoly, eye, pexquo, zeros, SparseMatrix
from sympy.abc import t

import sys
import os

sys.path.append(os.path.dirname(__file__))

import sparsedet

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
    s_matrix = sparsedet.SparseMatrix.zeros(len(hom_generators), len(hom_generators))

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

            elif i < j < hom_generators[i] < hom_generators[j]:
                separating_strands = abs(braid[i]) - abs(braid[j])
                if separating_strands == 1:
                    s_matrix[j, i] = -1
                elif separating_strands == -1:
                    s_matrix[i, j] = 1

    for i in range(len(hom_generators) - 1, 0, -1):
        if hom_generators[i] == 0:
            s_matrix = s_matrix.row_del(i)
            s_matrix = s_matrix.col_del(i)
    return s_matrix

def _get_reduced_burau_matrices(n):
    """
    Generates reduced Burau matrices for braid group B_n
    :param n: number of strands of braid group
    :return: dictionary "reps" of reduced Burau matrices with reps[i] the Burau matrices for the ith braid generator
    """
    assert n >= 2
    reps = {}
    if n == 2:
        reps[1], reps[-1] = eye(1), eye(1),
        reps[1][0, 0], reps[-1][0, 0]  = -t, -1/t
    else:
        reps[1], reps[-1] = eye(n - 1), eye(n - 1)
        reps[1][0, 0], reps[-1][0, 0] = -t, -1/t
        reps[1][0, 1], reps[-1][0, 1] = 1, 1/t
        for i in range(2, n - 1):
            reps[i], reps[-i] = eye(n - 1), eye(n - 1)
            reps[i][i - 1, i - 2:i + 1] = [[t, -t, 1]]
            reps[-i][i - 1, i - 2:i + 1] = [[1, -1/t, 1/t]]
        reps[n - 1], reps[-n + 1] = eye(n - 1), eye(n - 1)
        reps[n - 1][n - 2, n - 3:n - 1] = [[t, -t]]
        reps[-n + 1][n - 2, n - 3:n - 1] = [[1, -1/t]]
    return reps

def burau_rep(braid):
    reps = _get_reduced_burau_matrices(braid.n_strands)
    matrix = eye(braid.n_strands - 1)
    for gen in braid:
        matrix = matrix * reps[gen]
    return matrix

def _normalise_laurent(polynomial, symbol):
    """
    multiplies laurent polynomial by +/- symbol^n so it has no negative powers of symbol, and non zero positive constant
    :param polynomial: sympy Poly or PurePoly
    :param symbol: sympy symbol
    :return: normalised polynomial, sympy Poly or PurePoly
    """
    polynomial = polynomial.as_poly(t, 1/t)
    negative_degree = polynomial.degree(1/symbol)
    polynomial = (polynomial.as_expr() * symbol**negative_degree).as_poly(symbol)
    if negative_degree == 0:
        min_degree = min([mono[0] for mono in polynomial.monoms()])
        polynomial = (polynomial.as_expr() * symbol**(-min_degree)).as_poly(symbol)

    if polynomial.nth(0) < 0:
        polynomial = - polynomial
    return polynomial


def burau_to_alexander(matrix):
    a, b = matrix.shape
    n = a + 1
    matrix = eye(a) - matrix
    matrix = matrix.expand()
    alex_poly = matrix.det()
    alex_poly = alex_poly * (1 - t)
    alex_poly = _normalise_laurent(alex_poly, t)
    alex_poly = pexquo(alex_poly.as_expr(), ( 1 - t**n)).expand()
    alex_poly = PurePoly(alex_poly, t)
    if alex_poly.nth(0) < 0:
        alex_poly = - alex_poly
    return alex_poly


def decompose_poly_matrix(poly_matrix):
    poly_matrix = sympy.Matrix(poly_matrix)
    rows = poly_matrix.rows
    cols = poly_matrix.cols
    matrix_dict = {}

    for row in range(rows):
        for col in range(cols):
            poly = sympy.PurePoly(poly_matrix[row,col],t,1/t)
            for degrees, coeff in zip(poly.monoms(), poly.coeffs()):
                assert degrees[0] == 0 or degrees[1] == 0, "mixed expression!!"
                if degrees[0] > 0:
                    degree = degrees[0]
                elif degrees[1] > 0:
                    degree = - degrees[1]
                else:
                    degree = 0
                if not degree in matrix_dict.keys():
                    matrix_dict[degree] = zeros(rows, cols)

                matrix_dict[degree][row, col] = coeff
    return matrix_dict

def min_max_deg_degrees(poly_matrix):
    matrix_dict = decompose_poly_matrix(poly_matrix)
    degrees = sorted(matrix_dict.keys())
    for degree in degrees:
        det = matrix_dict[degree].det()
        if det != 0:
            min_degree = degree
            break
    for degree in degrees[::-1]:
        det = matrix_dict[degree].det()
        if det != 0:
            max_degree = degree
            break
    return min_degree, max_degree


def seifert_to_alexander(seifert_matrix):
    n, m = seifert_matrix.shape
    assert n == m, "non square matrix received"
    matrix = SparseMatrix(t * seifert_matrix - seifert_matrix.transpose())
    alex_poly = t**(-n//2) * matrix.det(method = "berkowitz")
    alex_poly = PurePoly(sympy.expand(alex_poly))
    return alex_poly


def signature(seifert_matrix):
    n, m = seifert_matrix.shape
    assert n == m, "non square matrix received"
    matrix = seifert_matrix + seifert_matrix.transpose()
    characteristic_poly = matrix.charpoly(t)
    #remove 0 eigenvalues
    while characteristic_poly.nth(0) == 0:
        characteristic_poly = pexquo(characteristic_poly, PurePoly(t))

    sig = characteristic_poly.degree(t)
    negative_root_intervals = characteristic_poly.intervals(sup=0)
    for interval, multiplicity in negative_root_intervals:
        assert interval[1] <= 0, ("Error in signature computation, root of characteristic equation was not isolated" +
                                  " correctly below zero. Could be a sympy version problem (working in 1.7.1).")
        sig += -2 * multiplicity
    return sig

