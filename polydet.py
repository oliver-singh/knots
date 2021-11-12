from sympy import Matrix, Add
from itertools import permutations

class PolyMatrixDet:
    def __init__(self, matrix=None, factor=1):
        assert isinstance(matrix, Matrix), "matrix must be a sympy matrix object!"
        self.mat = matrix
        self.factor = Add(factor)

    @property
    def rows(self):
        return self.mat.rows

    @property
    def cols(self):
        return self.mat.rows

    def __repr__(self):
        return "PolyMatrixDet(\nmatrix={},\nfactor={})".format(self.mat.__repr__(), self.factor)

    def row_swap(self, i, j):
        for k in range(0, self.mat.cols):
            self.mat[i, k], self.mat[j, k] = self.mat[j, k], self.mat[i, k]
        self.factor *= -1
        return self

    def col_swap(self, i, j):
        for k in range(0, self.mat.rows):
            self.mat[k, i], self.mat[k, j] = self.mat[k, j], self.mat[k, i]
        self.factor *= -1
        return self

    def row_mult(self,i, multiple):
        for k in range(0, self.mat.cols):
            self.mat[i, k] *= multiple
        self.factor = (self.factor / multiple).expand()
        return self

    def col_mult(self,i, multiple):
        for k in range(0, self.mat.rows):
            self.mat[k, i] *= multiple
        self.factor = (self.factor / multiple).expand()
        return self

    def row_add(self, i, j, multiple=1):
        for k in range(0, self.mat.cols):
            self.mat[i, k] = (self.mat[i, k] + multiple * self.mat[j, k]).expand()
        return self

    def off_diag_row(self, i):
        row = list(self.mat[i, :i]) + list(self.mat[i, i + 1:])
        row = [item.expand() for item in row]
        return row

    def off_diag_col(self, i):
        col = list(self.mat[:i, i]) + list(self.mat[i + 1:, i])
        col = [item.expand() for item in col]
        return col


    def clear_rowcol(self, i):
        remaining_row = self.off_diag_row(i)
        remaining_col = self.off_diag_col(i)
        row_zero = all([item ==0 for item in remaining_row])
        col_zero = all([item == 0 for item in remaining_col])
        assert row_zero or col_zero, "not all zeros!"

        self.factor *= self.mat[i, i]
        self.mat.col_del(i)
        self.mat.row_del(i)
        return self


def calculate_poly_det(matrix):
    matrix = Matrix(matrix)
    assert matrix.rows == matrix.cols
    det = PolyMatrixDet(matrix=matrix)
    determinant = 1

    while matrix.rows > 0:
        print("iteration")
        print(any([item == 1 for item in det.mat[:,0]]))
        top_val = det.mat[0,0].expand()
        if top_val == 0:
            rest_of_row = det.off_diag_row(0)
            if all([item == 0 for item in rest_of_row]):
                determinant = 0
                break
            else:
                j=0
                while rest_of_row[j] == 0:
                    j += 1
                det.row_swap(0, j+1)
                continue

        for j in range(1,det.rows):
            left_val = det.mat[j,0]
            if left_val.expand() != 0:
                det.row_add(j, 0, multiple=-left_val/top_val)
        det.clear_rowcol(0)
    determinant *= det.factor
    print("done")
    return determinant



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



def _minor_sign(coords):
    positive = True
    for coord in coords:
        row_even = (coord[0] % 2) == 0
        col_even = (coord[1] % 2) == 0
        if row_even != col_even:
            positive = not positive
    if positive:
        sign = 1
    else:
        sign = -1
    return sign

def iter_minor_coords(n, m):
    for col_coords in permutations(range(m)):
        yield [[row, col] for row, col in zip(range(n), col_coords)]

def calc_det(matrix):
    def _iter_minor_coords(coords_so_far, n,m):
        if len(coords_so_far == 11):
            yield coords_so_far













