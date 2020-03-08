from Tableauxes import *
import numpy as np
from sympy import Rational, expand, Matrix, pretty
from sympy.combinatorics import Permutation
from itertools import permutations
from numba import jit

flatted = lambda x: [j for i in x for j in i]
cycle_format = lambda permutation: "\nPermutation" + "".join(map(lambda x: str(tuple(x)), permutation.cyclic_form))
np.set_printoptions(precision = 3, suppress = True, linewidth = 100)

def get_matrix(matrix, show_type = False):
    if not show_type:
        if matrix.dtype == 'complex128':
            return '\n' + '\n'.join('\t'.join("{:.2g}".format(x) for x in y) for y in matrix)
        return '\n' + '\n'.join('\t'.join('%0.3f' % x for x in y) for y in matrix)
    else:
        # from texttable import Texttable
        # table = Texttable(0)
        # table.set_deco(Texttable.HEADER)
        matrix = expand(Matrix(matrix))#
        # n = matrix.shape[0]
        # for row in range(n):
        #     print("asdsad", pprint(matrix[row]))
        #     s = str(matrix[row])[1:-1]  # " ".join(map(lambda x: str(round(x, 4)), matrix[row]))#
        #     s = s.replace("sqrt(", "√").replace(")", "").split(" ")
        #     table.add_row(list(filter(lambda x: x, s)))  #
        #     print(s)
        # return '\n' + table.draw( )
        return '\n' + pretty(matrix)

class FiniteGroup:
    elements = {}

    def __iter__(self):
        for elem_arr_form, perm_el in self.elements.items():
            yield elem_arr_form, perm_el

    def __repr__(self):
        return f"{type(self).__name__}({self.n})"

    def __getitem__(self, element_array_form):
        return self.elements[element_array_form]
    #
    # def __setitem__(self, element_array_form, matrix_repr):
    #     self.elements[element_array_form].representation = matrix_repr

    def _display_calculated_elemets_representations(self):
        t = 0
        not_calc = []
        self.output = repr(self)
        for elem_arr_form, perm_el in self:
            if not hasattr(perm_el, 'representation'):
                not_calc.append(perm_el)
                t += 1
            else:
                self.output += '\n' + cycle_format(perm_el) + get_matrix(perm_el.representation, self.show_exact_values)
        self.output += f"\n\n{t} elements not calculated!"
        self.output += f"\nNot calculated elements:\n{not_calc}\n" if t > 0 else '\n'
        self.output += '_'*50
        print(self.output)
        return self.output

class SymmetricGroup(FiniteGroup):
    def __init__(self, N: int):
        assert isinstance(N, int) and N > 0
        self.n = N
        self.elements = {}
        for x in permutations(range(1, self.n + 1)):
            el = Permutation((0,) + x, size = self.n + 1)
            temp = tuple(el._array_form)
            self.elements[temp] = el

    def find_simple_transps_reprs(self):
        n = len(self.standart_tableauxes)
        conjugate = self.tableauxes.find_conjugate_transpositions( )
        axial_distances = self.tableauxes.find_axial_distances( )
        for x in range(1, self.n):  # going thru transpositions x, x+1
            m = np.zeros((n, n), dtype = Rational if self.show_exact_values else np.float)  # create matrix a filled zeroes
            conj = conjugate.get(x)
            axl_dist = axial_distances.get(x)
            for index in range(n):
                block = conj.get(index)
                diag, non_diag = axl_dist.get(index)
                if m[index][index] == 0:
                    m[index][index] = -diag
                if block:
                    m[index][block] = m[block][index] = non_diag
                    m[block][block] = diag
            self.simple_transpositions[Permutation(x, x + 1, size = self.n + 1)] = m
            try:
                self.elements[tuple(Permutation(x, x + 1, size = self.n + 1)._array_form)].representation = m
            except KeyError:
                pass

    def find_remain_simple_transps(self):
        def calculate_transposition(x, y):
            from numpy.linalg import multi_dot
            p = Permutation(x, y, size = self.n + 1)
            if p in self.simple_transpositions:
                return self.simple_transpositions[p]
            if x == 1:
                if Permutation(1, y - 1, size = self.n + 1) in self.simple_transpositions:
                    t1 = self.simple_transpositions[Permutation(1, y - 1, size = self.n + 1)]
                    return multi_dot([t1, self.simple_transpositions[Permutation(y - 1, y, size = self.n + 1)], t1])
                t2 = calculate_transposition(1, y - 1)
                return multi_dot([t2, self.simple_transpositions[Permutation(y - 1, y, size = self.n + 1)], t2])
            t3 = calculate_transposition(1, x)
            return multi_dot([t3, calculate_transposition(1, y), t3])

        for i in range(1, self.n):
            for j in range(i + 1, self.n + 1):
                if not (Permutation(i, j, size = self.n + 1) in self.simple_transpositions):
                    m = calculate_transposition(i, j)
                    self.simple_transpositions[Permutation(i, j, size = self.n + 1)] = m
                    if tuple(Permutation(i, j, size = self.n + 1)._array_form) in self.elements:
                        self.elements[tuple(Permutation(i, j, size = self.n + 1)._array_form)].representation= m
                        #self.group[tuple(Permutation(i, j, size = self._n + 1)._array_form)] = m

class AlternatingGroup(SymmetricGroup):
    def __init__(self, N: int):
        assert isinstance(N, int) and N > 0
        self.n = N
        self.elements = {}
        for x in permutations(range(1, self.n + 1)):
            el = Permutation((0,) + x, size = self.n + 1)
            if el.is_even:
                self.elements[tuple(el._array_form)] = el

    def find_simple_transps_reprs(self):
        if not self._shape.is_self_conjugated():
            super(AlternatingGroup, self).find_simple_transps_reprs()
            return

        def reflect_matrix(m):
            result = np.zeros((2 * n, 2 * n), datatype)
            for i in range(n):
                for j in range(n):
                    result[i][j] = m[i][j]
                    result[2 * n - i - 1][2 * n - j - 1] = -m[i][j] if i == j else m[i][j]
            return result

        self.simple_transpositions = {}
        n = len(self.standart_tableauxes) //2
        self.tableauxes.find_conjugate_transpositions( )
        conjugate = self.tableauxes.restrict()
        axial_distances = self.tableauxes.find_axial_distances( )
        mult, is_complex = self._shape.multiplier()
        datatype = Rational if self.show_exact_values else np.complex if is_complex else np.float
        #print(*conjugate.items( ), sep = "\n")
        for x in range(1, self.n):  # going thru transpositions x, x+1
            m = np.zeros((n, n), datatype) #create matrix filled with zeroes
            conj = conjugate[x]
            axl_dist = axial_distances.get(x)
            state = False
            for index in range(n):
                diag, non_diag = axl_dist.get(index)
                if m[index, index] == 0:
                    m[index, index] = -diag
                if index in conj:
                    if conj[index]:
                        block = conj.get(index)
                        m[block, block] = diag
                        m[index, block] = m[block, index] = non_diag
                    else:
                        block = conj.get(str(index))
                        state = True
                        m[block, block] = diag
                        m[index, block] = non_diag * mult
                        m[block, index] = non_diag * mult ** (-1)
            if x == 1:
                p12 = m
                self.simple_transpositions[Permutation(1, 2, size = self.n + 1)] = reflect_matrix(m)
            else:
                perm = Permutation(x, x + 1, size = self.n + 1)
                perm_array_form = tuple((Permutation(1, 2, size = self.n + 1) * perm)._array_form)
                if not state:
                    self.simple_transpositions[perm] = reflect_matrix(m)
                    self.elements[perm_array_form].representation = reflect_matrix(p12.dot(m))
                else:
                    m = reflect_matrix(m)
                    self.simple_transpositions[perm] = m#.dot(self.simple_transpositions[Permutation(1, 2, size = self.n + 1)])
                    self.elements[perm_array_form].representation = m

class Representation:

    def __init__(self, group_instance, get_exact_output = False, shape = None):
        self.group = group_instance
        self.show_exact_values = self.group.show_exact_values = get_exact_output
        if not shape:
            self._n = self.group.n
            _Partitions = Partitions(self._n, isinstance(self.group, AlternatingGroup))
            print(_Partitions)
            self._shape = self.group._shape = _Partitions[int(input("\nInput partition index: "))]
        else:
            self._n = self.group.n
            self._shape = self.group._shape = Partition(shape, isinstance(self.group, AlternatingGroup)) if isinstance(shape, (tuple, list)) and sum(shape) == self.group.n else shape

        self.group.tableauxes = Standart_tableauxes(self._shape)
        self.group.standart_tableauxes = self.group.tableauxes.get_standart( )
        print(self.group.tableauxes)

        self.find_simple_transps_reprs( )

    def find_simple_transps_reprs(self):
        self.group.simple_transpositions = {}
        self.group.find_simple_transps_reprs()
        self.group.find_remain_simple_transps()

    def find_remain_elems(self, parallel_type = None):
        if not parallel_type:
            tuple(self.compute_element(elem) for elem in self.group.elements)
        else:
            from multiprocessing import Manager
            self.d = Manager( ).dict( )
            if parallel_type == 'process':
                from multiprocessing import Pool as ProcessPool
                pool = ProcessPool(processes = 5)
            elif parallel_type == 'thread':
                from multiprocessing.dummy import Pool as ThreadPool
                pool = ThreadPool(processes = 3)
            else:
                raise AttributeError(f'{parallel_type} in not supported!')
            _pool = pool.map(self.compute_element_parallel, self.group.elements)
            for el, repres in self.d.items( ):
                self.group.elements[el].representation = repres
        return 'Function ends work'

    def _display_calculated_elemets_representations(self):
        return self.group._display_calculated_elemets_representations()


    def compute_element(self, element_array_form):
        permutation_element = self.group[element_array_form]
        permutation_element.transp_factor = [Permutation([x], size = permutation_element.size) for x in permutation_element.transpositions( )[::-1]]

        if len(permutation_element.transp_factor) > 1:
            permutation_element.representation = np.linalg.multi_dot([self.group.simple_transpositions[x] for x in permutation_element.transp_factor])
        elif len(permutation_element.transp_factor) == 0:
            permutation_element.representation = np.eye(len(self.group.standart_tableauxes))


    def compute_element_parallel(self, element_array_form):
        permutation_element = self.group[element_array_form]
        permutation_element.transp_factor = [Permutation([x], size = permutation_element.size) for x in permutation_element.transpositions( )[::-1]]

        if len(permutation_element.transp_factor) > 1:
            self.d[tuple(permutation_element._array_form)] = np.linalg.multi_dot([self.group.simple_transpositions[x] for x in permutation_element.transp_factor])
        elif len(permutation_element.transp_factor) == 0:
            self.d[tuple(permutation_element._array_form)] = np.eye(len(self.group.standart_tableauxes))


    def __getitem__(self, element_array_form):

        if not hasattr(self.group, 'simple_transpositions'):
            self.group.find_simple_transps_reprs( )

        temp_perm = Permutation(element_array_form._array_form, size = self._n + 1) if isinstance(element_array_form,
                                                                                                  Permutation) \
            else Permutation(element_array_form, size = self._n + 1)
        if tuple(temp_perm._array_form) != element_array_form:
            element_array_form = tuple(temp_perm._array_form)

        assert element_array_form in self.group.elements, f"Element {self.group[element_array_form]} not in {repr(self.group)}"

        permutation_element = self.group[element_array_form]
        if hasattr(permutation_element, 'representation'):
            return permutation_element.representation

        if not hasattr(permutation_element, 'transp_factor'):
            permutation_element.transp_factor = [Permutation([x], size = permutation_element.size) for x in
                                                 permutation_element.transpositions( )[::-1]]

        if len(permutation_element.transp_factor) > 1:
            permutation_element.representation = np.linalg.multi_dot(
                [self.group.simple_transpositions[x] for x in permutation_element.transp_factor])
        elif len(permutation_element.transp_factor) == 0:
            permutation_element.representation = np.eye(len(self.group.standart_tableauxes))

        try:
            return permutation_element.representation
        except AttributeError:
            pass


if __name__ == "__main__":
    from time import time
    S = SymmetricGroup(5)
    A = AlternatingGroup(5)
    repres = Representation(S, shape = (3, 1, 1), get_exact_output = True)

    t = time( )
    print(repres.find_remain_elems(parallel_type = None))
    #print(repres.find_remain_elems_parallel(parallel_type = 'process'))
    #print(repres.find_remain_elems_parallel(parallel_type = 'thread'))
    t = time( ) - t
    repres._display_calculated_elemets_representations()
    print('Time:', t)

    # # m1 = SGR.get_matrix_represents_permutation([5, 4, 2, 1, 3, 6])

