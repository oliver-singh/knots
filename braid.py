
import sys
import os

sys.path.append(os.path.dirname(__file__))

import invariants

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
    def from_string(cls, braid_string, n_strands=None):
        braid_string = braid_string.replace(";", ",")
        braid_string = braid_string.replace("{", "")
        braid_string = braid_string.replace("}", "")
        braid_string = braid_string.replace(" ", "")
        braid_list = braid_string.split(",")
        braid_list = [int(gen) for gen in braid_list]
        return cls(braid_list, n_strands)

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

    def __neg__(self):
        braid_list = [-gen for gen in self._braid_list]
        return Braid(braid_list, n_strands=self.n_strands)

    def __mul__(self, other):
        if isinstance(other, Braid):
            n_strands = max(self.n_strands, other.n_strands)
            braid_list = self._braid_list + other._braid_list
            product = Braid(braid_list, n_strands=n_strands)
        elif isinstance(other, int):
            product = Braid()
            if other < 0:
                other = abs(other)
                self = - self
            for i in range(other):
                product += self
        else:
            raise TypeError("unsupported operand type(s) for *: 'Braid' and '{}'".format(type(other)))
        return product

    def __rmul__(self, other):
        assert isinstance(other, int)
        return self * other

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
        braid_list = self._braid_list

        def check_unobstructed(position, gen):
            abs_gens = [abs(gen) for gen in braid_list]
            if position > last_position[gen]:
                intermediate_braids = abs_gens[last_position[gen]+1:position]
            else:
                intermediate_braids= abs_gens[last_position[gen]+1:] + abs_gens[:position]
            ob_above = gen + 1 in intermediate_braids
            ob_below = gen - 1 in intermediate_braids
            ob = gen in intermediate_braids
            obstructed = ob_above or ob_below or ob
            return obstructed

        last_position = {}
        last_sign = {}
        i = 0
        while i < 2 * len(braid_list):
            current_position = i % len(braid_list)
            gen = abs(braid_list[current_position])
            sign = braid_list[current_position] // gen
            if gen in last_position.keys():
                if current_position == last_position[gen]:
                    i += 1
                    continue
                if sign != last_sign[gen]:
                    if not check_unobstructed(current_position, gen):
                        first = min(current_position, last_position[gen])
                        last = max(current_position, last_position[gen])
                        braid_list.pop(last)
                        braid_list.pop(first)

                        i=0
                        last_position = {}
                        last_sign = {}
                        continue
            last_position[gen] = current_position
            last_sign[gen] = sign
            i += 1
        return Braid(braid_list, n_strands=self.n_strands)

    def super_simplify(self):

        def _generator_bump(strand, generator):
            abs_gen = abs(generator)
            sign = generator // abs_gen
            over_under = None
            if abs_gen == strand:
                strand += 1
                if sign == 1:
                    over_under = "under"
                else:
                    over_under = "over"
            elif abs_gen == strand - 1:
                strand += -1
                if sign == -1:
                    over_under = "under"
                else:
                    over_under = "over"
            return strand, over_under

        def _find_disk():
            over_under = []
            obs = False
            for k, lower, upper in zip(range(i + 1, j-1), lower_s[:-1], upper_s[:-1]):
                sign = braid_list[k]//abs(braid_list[k])
                abs_gen = abs(braid_list[k])
                print("********************************")
                print("i:",i)
                print("braid:",braid_list)
                print("over_under",over_under)
                print("upper:",upper_s)
                print("lower",lower_s)

                if abs_gen == lower - 1:
                    if sign == 1:
                        over_under = ["under"] + over_under
                    else:
                        over_under = ["over"] + over_under
                elif abs_gen == upper:
                    if sign == 1:
                        over_under = over_under + ["over"]
                    else:
                        over_under = over_under + ["under"]
                elif abs_gen == lower:
                    if over_under[0] == "over" and sign == -1:
                        obs = True
                        break
                    if over_under[0] == "under" and sign == +1:
                        obs = True
                        break
                    over_under = over_under[1:]
                elif abs_gen == upper - 1:
                    if over_under[0] == "over" and sign == 1:
                        obs = True
                        break
                    if over_under[0] == "under" and sign == -1:
                        obs = True
                        break
                    over_under = over_under[:-1]
                elif lower + 1 < abs_gen < upper:
                    if over_under[abs_gen - lower -1] == over_under[abs_gen - lower]:
                        obs = True
                        break
                    elif over_under[abs_gen - lower -1] == "under" and sign == -1:
                        obs = True
                        break
                    elif over_under[abs_gen - lower - 1] == "upper" and sign == 1:
                        obs = True
                        break
                    else:
                        a, b = over_under[abs_gen - lower], over_under[abs_gen - lower - 1]
                        over_under[abs_gen - lower], over_under[abs_gen - lower - 1] = b, a
            return obs

        i = -1
        braid_list = self._braid_list
        while i < len(braid_list) - 1:
            i = (i + 1)
            position = i % len(braid_list)
            base_gen = abs(braid_list[position])
            lower_s = [base_gen]
            upper_s = [base_gen + 1]
            base_sign = braid_list[position]//base_gen

            base_strands_all_over = True
            base_strands_all_under = True

            j=position
            while (position - 1 - j) % len(braid_list) != 0 and lower_s[-1] < upper_s[-1]:

                j = (j + 1) % len(braid_list)
                lower_end = lower_s[-1]
                upper_end = upper_s[-1]


                if base_strands_all_over or base_strands_all_under:
                    if lower_end + 1 == upper_end:
                        base_strands_all_over = True
                        base_strands_all_under = True
                else:
                    break
                current_gen = braid_list[j]

                lower_end, lower_over_under = _generator_bump(lower_end, current_gen)
                upper_end, upper_over_under = _generator_bump(upper_end, current_gen)

                if lower_end < upper_end:
                    if lower_over_under == "over" or upper_over_under == "over":
                        base_strands_all_under = False
                    if (lower_over_under == "under") or (upper_over_under == "under"):
                        base_strands_all_over = False

                lower_s += [lower_end]
                upper_s += [upper_end]
            if braid_list[j]//abs(braid_list[j]) == base_sign:
                continue
            if lower_s[-1] > upper_s[-1]:
                obstructed = True
                if base_strands_all_over or base_strands_all_under:
                    """
                    print("*************************cancellation************************")
                    print("position", position)
                    print("braid:", braid_list)
                    print("upper strand: ", upper_s)
                    print("lower strand:", lower_s)
                    print("all under:", base_strands_all_under)
                    print("all over:", base_strands_all_over)
                    """
                    first, last = sorted([position, j])
                    braid_list.pop(last)
                    braid_list.pop(first)
                    i = -1

        return Braid(braid_list, n_strands=self.n_strands)

    # invariants
    def seifert_matrix(self): # maybe cache this one later - store as property of braid, and add a clear_cache method
        # if not hasattr(self, "_seifert_matrix") or self._seifert_matrix is None:
        #    self._seifert_matrix = invariants.seifert_matrix(self)
        # return self._seifert_matrix

        matrix = invariants.seifert_matrix(self)
        return matrix

    def signature(self):
        matrix = self.seifert_matrix()
        sig = invariants.signature(matrix)
        return sig

    def burau_rep(self):
        matrix = invariants.burau_rep(self)
        return matrix

    def alexander_poly(self, method="burau"):
        if method == "burau":
            matrix = self.burau_rep()
            poly = invariants.burau_to_alexander(matrix)
        elif method == "seifert":
            matrix = self.seifert_matrix()
            poly = invariants.seifert_to_alexander(matrix)
        return poly


def empty(n):
    return Braid(n_strands=n)


def braid_range(*args, positive=True):
    """
    Returns a braid using list given by in built python range function. Excludes 0.
    :param args: Args to pass to range
    :param positive: If false negates the whole braid (negative braids may also be obtained by range(-4,-1)
    :return: a braid with a generator for each
    """
    braid_list = list(range(*args))
    braid_list = [gen for gen in braid_list if gen != 0]
    if not positive:
        braid_list = [-gen for gen in braid_list]
    return Braid(braid_list)


def full_twist(n):
    positive = True
    if n < 0:
        positive = False
        n = abs(n)
    twist = braid_range(n) ** n

    if not positive:
        twist = twist**(-1)
    return twist

class Sbraid(Braid):
    def __init__(self):
        super().__init__()