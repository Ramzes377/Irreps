from Partitions import *
from sympy import sqrt, Rational


class Tablue:

    def __init__(self, shape, tablue):
        assert isinstance(shape, Partition) and isinstance(tablue, (list, tuple))
        self.shape = shape
        self.n = self.shape.n
        self.tablue = tablue
        self.conjugate = {}

    def get_tablue_as_multi_level_list(self):
        if not hasattr(self, 'tablue_as_mll'):
            self.tablue_as_mll = []
            shift = 0
            for lenght in self.shape:
                self.tablue_as_mll.append(list(self.tablue[shift: lenght + shift]))
                shift += lenght
        return self.tablue_as_mll

    def conjugated(self, x, y = None):
        if not self.conjugate.get(x):
            y = x + 1 if not y else y
            tablue = list(self.tablue[:])
            x_index, y_index = tablue.index(x), tablue.index(y)
            tablue[x_index], tablue[y_index] = y, x
            t = Tablue(self.shape, tablue)
            self.conjugate[x] = t
        return self.conjugate.get(x)

    def axial_distance(self, x, y):
        assert 1 <= x <= self.n and 1 <= y <= self.n
        tablue = self.get_tablue_as_multi_level_list( )
        x_row = y_row = None
        for index, row in enumerate(tablue):
            if x_row and y_row:
                break
            if not x_row and x in row:
                x_row, x_column = index, row.index(x)
            if not y_row and y in row:
                y_row, y_column = index, row.index(y)
        return x_column + y_row - x_row - y_column

    def get_block_values(self, x, y = None):
        y = x + 1 if y == None else y
        axl_dist = self.axial_distance(x, y)
        diagonal = Rational(1, axl_dist)
        non_diagonal = sqrt(1 - Rational(1, axl_dist * axl_dist))
        return diagonal, non_diagonal

    def is_standart(self):
        t = self.shape.template()
        for i in range(self.n):
            case1, case2 = t.get(i), t.get(str(i))
            if (case1 and self[i] > self[case1]) or (case2 and self[i] > self[case2]):
                return False
        return True

    def __getitem__(self, item):
        return self.tablue[item]

    def __str__(self):
        # from texttable import Texttable
        # result = ""
        # n = len(self.tablue_as_mll[0])
        # for row in self.tablue_as_mll:
        #     table = Texttable(0)
        #     table.add_row(row)
        #     k = table.draw()
        #     for i in range(len(k) - 1, 0, -1):
        #         if k[i] == "|":
        #             #print(i)
        #             k = k[:i + 1]
        #             break
        #     result += k + "\n"
        # return result
        return "|\n".join(
            "| ".join(str(number) for number in string) for string in self.get_tablue_as_multi_level_list( )) + "\n"

    def __eq__(self, other):
        return self.shape == other.shape and tuple(self.tablue) == tuple(other.tablue)

    def get_grid(self):
        mll = self.get_tablue_as_multi_level_list()
        if not hasattr(self, 'grid'):
            self.grid = {}
            for i in range(len(mll)):
                level = mll[i]
                for j in range(len(level)):
                    self.grid[level[j]] = (i, j)
        return self.grid

    def restrict(self, restriction_length = 2):
        def reflect_and_concat(grid, remainder):
            new_grid = {}
            for x, pos in grid.items():
                new_grid[x] = (pos[1], pos[0])
            for x, pos in remainder:
                new_grid[x] = pos
            return new_grid

        grid = self.get_grid()
        remainder = []

        for i in range(self.n - restriction_length + 1, self.n + 1):
            remainder.append((i, grid.pop(i)))

        reflected_grid = reflect_and_concat(grid, remainder)
        seq = [key for key, value in sorted(sorted(reflected_grid.items( ), key = lambda x: x[0]), key = lambda x: x[1])]
        new_tablue = Tablue(self.shape, seq)
        return new_tablue

class Standart_tableauxes:
    def __init__(self, shape):
        assert isinstance(shape, (Partition, list, tuple))
        self._shape = shape if isinstance(shape, Partition) else Partition(shape)
        self._n = self._shape.n
        self._shape.template( )
        self.conjugate_transpositions = {}
        self.get_standart()

    def get_standart(self):
        '''Find standart Young tableauxes for given partition\n'''
        tablue = Tablue(self._shape, list(range(1, self._n + 1)))
        elems = [tablue.tablue_as_mll.index(row) for index in range(1, self._n + 1) for row in tablue.get_tablue_as_multi_level_list() if index in row]
        temp = map(lambda p: self.t_from_seq(p), self._helper_unique_permutations(elems))
        self.tableauxes = filter(lambda t: t.is_standart( ), temp)
        for i in range(1, self._n):
            self.tableauxes = sorted(self.tableauxes, key = lambda f: f[i], reverse = True)
        return self.tableauxes

    def t_from_seq(self, seq):
        yamanuchi = dict(sorted(zip(range(1, len(seq) + 1), seq), key = lambda x: x[1]))
        tablue = Tablue(self._shape, tuple(yamanuchi.keys( )))
        tablue.yamanuchi = yamanuchi
        return tablue

    def find_conjugate_transpositions(self):
        self.conjugate_transpositions[1] = {}
        for n in range(2, self._n):
            temp = {}
            for i, tablue in enumerate(self.tableauxes):
                try:
                    t = self.tableauxes[:i].index(tablue.conjugated(n))
                    temp[t] = i
                except ValueError:
                    continue
            self.conjugate_transpositions[n] = temp
        return self.conjugate_transpositions

    def find_axial_distances(self):
        result = {}
        for x in range(1, self._n):
            result[x] = {t: self.tableauxes[t].get_block_values(x) for t in range(len(self.tableauxes))}
        return result

    def __iter__(self):
        assert self.tableauxes
        for tablue in self.tableauxes:
            yield tablue

    def __str__(self):
        return f"\n{len(self.tableauxes)} standart Young's frames:\n" + "\n".join(map(str, self))


    def restrict(self):
        n = len(self.tableauxes) // 2
        result = {}
        for idx in range(1, self._n - 2):
            result[idx] = {x: y for x, y in self.conjugate_transpositions[idx].items( ) if x < n - 1}

        for idx in range(self._n - 2, self._n):
            temp = {}
            for x, y in self.conjugate_transpositions[idx].items( ):
                if x < n - 1:
                    if y < n:
                        temp[x] = y
                    elif y >= n:
                        tablue = self.tableauxes[x]
                        restricted_and_flipped = tablue.restrict(2)
                        try:
                            index = self.tableauxes.index(restricted_and_flipped)
                        except ValueError:
                            restricted_and_flipped = tablue.restrict(1)
                            index = self.tableauxes.index(restricted_and_flipped)
                        if index < n:
                            temp[x] = None
                            temp[str(x)] = index
            result[idx] = temp
        self.conjugate_transpositions = result
        return self.conjugate_transpositions

    def _helper_unique_permutations(self, elements):
        def unique_permutations(elements, t = 1):
            if len(elements) == 1:
                yield (elements[0],)
            unique_elements = set(elements)
            for first_element in unique_elements:
                if first_element < t:
                    remaining_elements = list(elements)
                    remaining_elements.remove(first_element)
                    for sub_permutation in unique_permutations(remaining_elements, t + 1):
                        yield (first_element,) + sub_permutation
        return unique_permutations(elements)


if __name__ == "__main__":
    from time import time
    shape = (3, 1, 1)
    print(f"Example\nFind standart tables for partition {shape}.\n")
    t = time()
    y = Standart_tableauxes(shape)
    print(y)
    print(time() - t)
    #print(y.find_conjugate_transpositions( ))


''' 




'''