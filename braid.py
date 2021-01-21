
class Braid:
    def __init__(self, braid_list=None, n_strands=None):
        if braid_list is None:
            braid_list = []
        assert isinstance(braid_list, list), ("braid must be formatted as a list, " +
                                              "unless using Braid.from_string")
        assert all(isinstance(generator, int) for generator in braid_list)
        assert all(generator != 0 for generator in braid_list)
        self._braid_list = braid_list

        if n_strands is None:
            n_strands = self.max_generator() + 1
        self.n_strands = n_strands

    def max_generator(self):
        if self._braid_list == []:
            max_gen = -1
        else:
            max_gen = max([abs(gen) for gen in self._braid_list])
        return max_gen

    @classmethod
    def from_string(cls, braid_string):
        knot = cls()
        knot.braid_string = braid_string
        return knot

    @classmethod
    def empty(cls, n_strands):
        return cls(n_strands=n_strands)

    def __str__(self):
        return ";".join([str(generator) for generator in self._braid_list])

    @property
    def n_strands(self):
        return self._n_strands

    @n_strands.setter
    def n_strands(self, n):
        assert n >= 0, "number of strands cannot be negative"
        assert n > self.max_generator(), (f"{n} strands not enough for braid " +
                                          f"on {self.max_generator()} generators")
        self._n_strands = n

    def __repr__(self):
        return "Braid({}, n_strands={})".format(self._braid_list, self.n_strands)

    def __len__(self):
        return len(self._braid_list)

    def __getitem__(self, key):
        return self._braid_list[key]

    def __setitem__(self, key, value):
        if isinstance(value, int):
            maxval = value
            assert value != 0, "0 is not a valid generator"
        else:
            assert all(isinstance(gen, int) for gen in value), "braid generators must be int"
            maxval = max([abs(gen) for gen in value])
            assert all(gen != 0 for gen in value), "0 is not a valid generator"

        assert maxval < self.n_strands, (f"cannot use generator {maxval} on a braid with " +
                                         f"{self.n_strands} strands. You may increase the number of " +
                                         "strands using Braid.n_strands = {}".format(maxval + 1))
        self._braid_list[key] = value

    def _bump_gen(self, diff):
        braid_list = []

        for gen in self._braid_list:
            if gen < 0:
                braid_list += [gen - diff]
            else:
                braid_list += [gen + diff]

        return Braid(braid_list, self.n_strands + diff)

    def __add__(self, other):
        assert isinstance(other, Braid)
        return self * other._bump_gen(self.n_strands)

    def __mul__(self, other):
        if isinstance(other, Braid):
            n_strands = max(self.n_strands, other.n_strands)
            braid_list = self._braid_list + other._braid_list
            product = Braid(braid_list, n_strands=n_strands)
        elif isinstance(other, int):
            product = Braid()
            for i in range(other):
                product += self
        else:
            raise TypeError("unsupported operand type(s) for *: 'Braid' and '{}'".format(type(other)))
        return product

    def __rmul__(self, other):
        assert isinstance(other, int)
        product = Braid()
        for i in range(other):
            product += self
        return product

    def __pow__(self, power):
        assert isinstance(power, int)
        braid_list = self._braid_list
        if power == 0:
            braid_list = []
        elif power < 0:
            power = abs(power)
            braid_list = [-gen for gen in braid_list]
            braid_list = braid_list[::-1]
        return Braid(braid_list*power, n_strands=self.n_strands)

    def inverse(self):
        return self**-1

    def simplify(self):
        braid_list = [0]
        for gen in self._braid_list:
            if gen != -braid_list[-1]:
                braid_list += [gen]
            else:
                braid_list.pop(-1)
        braid_list.pop(0)
        while braid_list != [] and braid_list[0] == -braid_list[-1]:
            braid_list.pop(0)
            braid_list.pop(-1)

        return Braid(braid_list, n_strands=self.n_strands)


