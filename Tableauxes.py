from Partitions import *
from sympy import sqrt, Rational


# class unique_permutations:
#
#     class unique_element:
#         def __init__(self, value, occurrences):
#             self.value = value
#             self.occurrences = occurrences
#
#     def __init__(self, elements):
#         eset = set(elements)
#         self.list = [self.unique_element(i, elements.count(i)) for i in eset]
#         self.u = len(elements)
#
#     def __iter__(self):
#         def unique_perm(listunique, result_list, d):
#             if d < 0:
#                 yield tuple(result_list)
#             for i in reversed(listunique):
#                 if i.occurrences > 0:
#                     result_list[d] = i.value
#                     i.occurrences -= 1
#                     yield from unique_perm(listunique, result_list, d - 1)
#                     i.occurrences += 1
#         return unique_perm(self.list, [None] * self.u, self.u - 1)

def action_with_flatten(function, multi_level_list):
    return [function(x) for item in multi_level_list for x in item]


class Yamanuchi:
    def __init__(self, sequence = None, tablue = None, **kwargs):
        if kwargs.get('sequence', sequence):
            self._init_from_sequence(kwargs.get("sequence", sequence))
        elif kwargs.get('tablue', tablue):
            self._init_from_tablue(kwargs.get('tablue', tablue))
        else:
            raise AttributeError("Wrong arguments input")

    def _init_from_sequence(self, sequence):
        if isinstance(sequence, (list, tuple)):
            self.sequence = sequence
            self.yamanuchi = dict(sorted(zip(range(1, len(self.sequence) + 1), self.sequence), key = lambda x: x[1]))
        elif isinstance(sequence, dict):
            self.yamanuchi = dict(sorted(sequence.items( ), key = lambda x: x[1]))
            self.sequence = self.yamanuchi.values( )
        else:
            raise AttributeError('sequence must be list or dict type!')

    def _init_from_tablue(self, tablue):
        assert isinstance(tablue, Tablue), "tablue must be Tablue type!"
        self.tablue = tablue
        self.shape = tablue.shape
        from functools import reduce
        self.tablue_as_mll = tablue.get_tablue_as_multi_level_list( )
        n = reduce(lambda prev, cur: prev + cur, map(len, self.tablue_as_mll))
        self.sequence = [self.tablue_as_mll.index(row) for index in range(1, n + 1)\
                         for row in self.tablue_as_mll if index in row]
        self.yamanuchi = dict(sorted(zip(range(1, len(self.sequence) + 1), self.sequence), key = lambda x: x[1]))

    def get_yamanuchi_sequence(self):
        if not hasattr(self, 'yamanuchi_sequence'):
            self.yamanuchi_sequence = tuple(self.yamanuchi.values( ))
        return self.yamanuchi_sequence

    def get_shape(self):
        if not hasattr(self, 'shape'):
            t = self.get_yamanuchi_sequence( )
            self.shape = Partition(tuple(t.count(i) for i in set(t)))
        return self.shape

    def get_tablue_values(self):
        if not hasattr(self, 'tablue_values'):
            self.tablue_values = tuple(self.yamanuchi.keys( ))
        return self.tablue_values

    def get_tablue(self):
        if not hasattr(self, 'tablue'):
            shape = self.get_shape( )
            tablue = self.get_tablue_values( )
            self.tablue = Tablue(shape, tablue)
        return self.tablue

    def flip_and_concatenate(self, remainder):
        if not hasattr(self, 'tablue_as_mll'):
            self.tablue_as_mll = self.get_tablue( ).get_tablue_as_multi_level_list()

        temp = [[] for x in self.get_shape( )]
        copy = self.tablue_as_mll[:]
        for j in range(len(copy[0])):
            for i in range(len(copy)):
                try:
                    temp[i].append(copy[j].pop(0))
                except IndexError:
                    continue
        result = []
        for i in range(len(temp)):
            result.extend(temp[i] + copy[i])

        y = Yamanuchi(tablue = Tablue(self.get_shape( ), result)) + remainder
        return y.get_tablue( )

    def restrict(self, restiction_lenght = 2):
        restriction = self[:len(self) - restiction_lenght]
        remainder = self[len(self) - restiction_lenght:]
        y = Yamanuchi(sequence = restriction)
        shape = y.get_shape( )
        if restiction_lenght > 1 and not shape.is_self_conjugated( ):
            return self.restrict(1)
        elif restiction_lenght == 1:
            return None
        flipped = y.flip_and_concatenate(remainder)
        return flipped

    def __str__(self):
        return str(self.yamanuchi)

    def __len__(self):
        return len(self.yamanuchi)

    def __getitem__(self, key):
        items = sorted(self.yamanuchi.items( ), key = lambda x: x[0])
        return dict(items[key])

    def __add__(self, other):
        return Yamanuchi(sequence = {**self.yamanuchi, **other})


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
        # fy = lambda sequence: tuple(x[0] for x in sorted(tuple(enumerate(sequence, start = 1)), key = lambda x: x[1]))
        elems = Yamanuchi(tablue = Tablue(self._shape, list(range(1, self._n + 1)))).sequence
        # temp = map( lambda t: Tablue(self._shape, fy(t)), unique_permutations(elems) )
        temp = map(lambda p: Yamanuchi(sequence = p).get_tablue( ), self._helper_unique_permutations(elems))
        self.tableauxes = filter(lambda t: t.is_standart( ), temp)
        for i in range(1, self._n):
            self.tableauxes = sorted(self.tableauxes, key = lambda f: f[i], reverse = True)
        return self.tableauxes

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

    def calling_to_restrict(self):
        n = len(self.tableauxes) // 2
        result = {}
        for idx in range(1, self._n - 2):
            result[idx] = {x: y for x, y in self.conjugate_transpositions[idx].items( ) if x < n - 1}

        for idx in range(self._n - 2, self._n):
            temp = {}
            already_calc = []
            for x, y in self.conjugate_transpositions[idx].items( ):
                if x < n - 1 and x not in already_calc:
                    if  y > n - 1:
                        tablue = self.tableauxes[x]
                        yamanuchi = Yamanuchi(tablue = tablue)
                        flipped = yamanuchi.restrict()
                        try:
                            index = self.tableauxes.index(flipped)
                        except AttributeError:
                            index = None
                        if index and index < n:
                            temp[str(x)] = index
                            temp[x] = None
                            already_calc.append(index)
                    else:
                        temp[x] = y
            result[idx] = temp
        self.conjugate_transpositions = result
        return self.conjugate_transpositions

    def _helper_unique_permutations(self, elements):
        def unique_permutations(elements, t = 0):
            if len(elements) == 1:
                yield (elements[0],)
            unique_elements = set(elements)
            for first_element in unique_elements:
                if first_element > t:
                    continue
                remaining_elements = list(elements)
                remaining_elements.remove(first_element)
                for sub_permutation in unique_permutations(remaining_elements, t + 1):
                    yield (first_element,) + sub_permutation
        return unique_permutations(elements)


if __name__ == "__main__":
    from time import time
    shape = (3, 1, 1)
    print(f"Example\nFind standart tables for partition {shape}.\n")
    y = Standart_tableauxes(shape)
    print(y)
    #print(y.find_conjugate_transpositions( ))


''' 


def find_conjugate_transpositions2(self):
    t = {}
    t[1] = {}
    valid = lambda t1,t2, n: t1.index(n) == t2.index(n+1) and sum([abs(x-y) for x,y in zip(t1,t2)]) == 2
    for n in range(2, self._n):
        t[n] = {a: b for a, t1 in enumerate(self) for b, t2 in enumerate(self) if a < b and valid(t1.tablue, t2.tablue, n)}
    return t


def from_yamanuchi(sequence):
    tablue = [[] for x in range(len(set(sequence)))]
    for i, ys in enumerate(sequence, start = 1):
        tablue[ys].append(i)
    return [j for i in tablue for j in i]


old checking
        def valid(sequence):
            form = []; shift = 0
            for lenght in self._shape:
                temp = []
                for index in range(lenght):
                    temp.append(sequence[ index + shift ])
                form.append(temp)
                shift += lenght

            for i in range(len(form)):
                for j in range(len(form[i])):
                    try:
                        if form[i][j] > form[i+1][j] or form[i][j] > form[i][j+1]: 
                            return False
                    except IndexError:
                        try:
                            if form[i][j] > form[i][j+1]: 
                                return False
                        except IndexError:
                            try:
                                if form[i][j] > form[i+1][j]: 
                                    return False
                            except: 
                                continue
            #something.append(form)
            #print(something)
            return True


            def from_yamanuchi(sequence):
                tablue = [[] for x in range(len(shape))]
                for i, ys in enumerate(yamanuchi_sequence, start = 1):
                    tablue[ys].append(i)
                for i in tablue:
                    for j in i:
                        result.append(j)
                return result


            # def from_yamanuchi(sequence):
            #     result = []
            #     i = start = 0
            #     while len(result) < len(sequence):
            #         try:
            #             result.append(sequence.index(i, start) + 1)
            #             start = result[-1]
            #         except ValueError:
            #             i += 1
            #             start = 0
            #     return result

'''