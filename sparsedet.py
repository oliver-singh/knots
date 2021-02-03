from copy import deepcopy


def _index_value_dict_to_row_col_dict(index_value_dict):
    row_col_dict = {}
    for index, value in index_value_dict.items():
        if index[0] not in row_col_dict:
            row_col_dict[index[0]] = {index[1]: value}
        else:
            row_col_dict[index[0]][index[1]] = value
    return row_col_dict

class SparseMatrix:
    def __init__(self, row_col_dict, rows=None, cols=None):
        """
        index value dict looks like
        {(i,j):value, (i',j'):value'}
        """
        self._row_col_dict = row_col_dict
        self.rows = rows
        self.cols = cols

    @classmethod
    def from_indexed_values(cls, index_value_dict, rows=None, cols=None):
        row_col_dict = _index_value_dict_to_row_col_dict(index_value_dict)
        matrix = cls(row_col_dict, rows=rows, cols=cols)
        return matrix

    @classmethod
    def from_list_of_lists(cls, dense_mat, rows=None, cols=None):
        row_col_dict = {}
        for row_ix, row in enumerate(dense_mat):
            if not all(val == 0 for val in row):
                row_col_dict[row_ix] = {}
            for col_ix, value in enumerate(row):
                if value != 0:
                    row_col_dict[row_ix][col_ix] = value
        matrix = cls(row_col_dict, rows=rows, cols=cols)
        return matrix

    def index_value_dict(self):
        index_val_dict = {}
        for row_ix, row in self._row_col_dict.items():
            for col_ix, value in row.items():
                index_val_dict[(row_ix,col_ix)] = value
        return index_val_dict

    def max_non_zero_row(self):
        if self._row_col_dict == {}:
            max_row = -1
        else:
            max_row = max(self._row_col_dict.keys())
        return max_row

    def max_non_zero_col(self):
        col_union = [col for row_inx, row in self._row_col_dict.items() for col in row.keys()]
        if col_union == []:
            max_col = -1
        else:
            max_col = max(col_union)
        return max_col

    @property
    def rows(self):
        return self._rows

    @rows.setter
    def rows(self, n):
        if n is None:
            self._rows = self.max_non_zero_row() + 1
        else:
            assert isinstance(n, int)
            assert self.max_non_zero_row() < n
            self._rows = n

    @property
    def cols(self):
        return self._cols

    @cols.setter
    def cols(self, m):
        if m is None:
            self._cols = self.max_non_zero_col() + 1
        else:
            assert isinstance(m, int)
            assert self.max_non_zero_col() < m
            self._cols = m

    def __setitem__(self, key, value):
        assert len(key) == 2, "was expecting two indices but received {}.".format(key)
        row, col = key
        assert row < self.rows and col < self.cols, "index {},{} out of range for {}x{} matrix.".format(row, col, self.rows, self.cols)
        if value == 0:
            if row in self._row_col_dict and col in self._row_col_dict[row]:
                del self._row_col_dict[row][col]
                if self._row_col_dict[row] == {}:
                    del self._row_col_dict[row]
        else:
            if row not in self._row_col_dict:
                self._row_col_dict[row] = {col: value}
            else:
                self._row_col_dict[row][col] = value

    def __getitem__(self, key):
        row, col = key
        if isinstance(row, int) and isinstance(col, int):
            assert row < self.rows and col < self.cols, "index {},{} out of range for {}x{} matrix.".format(row, col, self.rows, self.cols)
            if row in self._row_col_dict and col in self._row_col_dict[row]:
                value = self._row_col_dict[key[0]][key[1]]
            else:
                value = 0
        elif isinstance(row, slice) or isinstance(col, slice):
            if isinstance(col, slice):
                if (col.start is None) or col.start < 0:
                    col_start = 0
                else:
                    col_start = col.start
                if col.stop is None or col.stop > self.cols:
                    col_stop = self.cols
                else:
                    col_stop = col.stop
                if col.step is not None:
                    raise NotImplementedError
            else:
                assert isinstance(col, int), "can only slice with slices or integers"
                assert 0 <= col < self.cols
                col_start = col
                col_stop = col + 1
            if isinstance(row, int):
                value = SparseMatrix.zeros(1, col_stop - col_start)
                for col_ix in self._row_col_dict[row]:
                    if col_start <= col_ix < col_stop:
                        value[0, col_ix - col_start] = self[row, col_ix]
            else:
                assert isinstance(row, slice), "can only slice with integers or slices"
                if (row.start is None) or row.start < 0:
                    row_start = 0
                else:
                    row_start = row.start
                if row.stop is None or row.stop > self.rows:
                    row_stop = self.rows
                else:
                    row_stop = row.stop
                if row.step is not None:
                    raise NotImplementedError

                value = SparseMatrix.zeros(row_stop - row_start, col_stop - col_start)
                for row_ix, _row in self._row_col_dict.items():
                    if row_start <= row_ix < row_stop:
                        for col_ix, val in _row.items():
                            if col_start <= col_ix < col_stop:
                                value[row_ix - row_start, col_ix - col_start] = val
        else:
            raise TypeError("can only slice with int or slices")
        return value

    @classmethod
    def zeros(cls, *args):
        arg_len = len(args)
        if arg_len == 0:
            matrix = cls({}, rows = 0, cols = 0)
        elif arg_len == 1:
            matrix = cls({}, rows=args[0], cols=args[0])
        elif arg_len == 2:
            matrix = cls({}, rows=args[0], cols=args[1])
        else:
            raise ValueError("too many values to unpack (expected 0, 1 or 2)")
        return matrix

    @classmethod
    def eye(cls, n):
        assert isinstance(n, int), "argument of eye must be an integer"
        matrix = cls.zeros(n)
        for i in range(n):
            matrix[i, i] = 1
        return matrix

    def __repr__(self):
        return ("SparseMatrix.from_list_of_lists(\n"
                + str(self)
                + ", rows={}, cols={})".format(self.rows, self.cols))

    def __add__(self, other):
        result = SparseMatrix(deepcopy(self._row_col_dict), self.rows, self.cols)

        for row in other._row_col_dict:
            for col in other._row_col_dict[row]:
                result[row, col] = result[row, col] + other[row, col]
        return result

    def to_dense(self):
        dense_mat = []
        for row in range(self.rows):
            dense_mat += [[0]*self.cols]

        for row_ix, row in self._row_col_dict.items():
            for col_ix, value in row.items():
                dense_mat[row_ix][col_ix] = value
        return dense_mat

    def __str__(self):
        return "[{}]".format(",\n ".join([str(row) for row in self.to_dense()]))

    def __mul__(self, other):
        if isinstance(other, SparseMatrix):
            result = SparseMatrix.zeros(self.rows, other.cols)
            assert self.cols == other.rows

            for row_l in self._row_col_dict:
                for row_r in other._row_col_dict:
                    for col_r in other._row_col_dict[row_r]:
                        result[row_l, col_r] += self[row_l, row_r]*other[row_r, col_r]
        else:
            result = SparseMatrix(deepcopy(self._row_col_dict), self.rows, self.cols)
            for row_ix in self._row_col_dict:
                result.row_mult(row_ix, other)
        return result

    def __rmul__(self, other):
        return self * other

    def swap_rows(self, i, j):
        self._row_col_dict[i], self._row_col_dict[j] = self._row_col_dict[j], self._row_col_dict[i]
        return self

    def del_row(self, i):
        new_dict = {}
        assert 0 <= i < self.rows, "row index {} out of range".format(i)
        for row_ix, row in self._row_col_dict.items():
            print(row_ix)
            if row_ix < i:
                new_dict[row_ix] = row
            elif row_ix > i:
                new_dict[row_ix-1] = row
        return SparseMatrix(new_dict, rows=self.rows - 1, cols=self.cols)

    def row_mult(self, i, multiple):
        if i in self._row_col_dict:
            row_i = self._row_col_dict[i].copy()
            for col in row_i:
                self[i, col] *= multiple
        return self

    def row_add_op(self, i, j, multiple):
        assert i < self.rows and j < self.rows, "index out of range"
        assert i != j, "cant perform row col between same rows"
        if j in self._row_col_dict:
            if not i in self._row_col_dict:
                self._row_col_dict[i] = {col_ix: value*multiple for col_ix, value in self._row_col_dict[j].items()}
            else:
                for col_ix, value in self._row_col_dict[j].items():
                    self[i, col_ix] += value*multiple
        return self

    def transpose(self):
        result = SparseMatrix.zeros(self.cols, self.rows)
        for row_ix, row in self._row_col_dict.items():
            for col_ix, value in row.items():
                result[col_ix, row_ix] = value
        return result




def test(*args, value):
    print(args)
    print(value)

a= {1:{1}}
b = a.copy