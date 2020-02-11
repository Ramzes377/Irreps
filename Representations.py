from Tableauxes import *
import numpy as np
from sympy import Rational, expand, Matrix, pretty
from sympy.combinatorics import Permutation
from itertools import permutations


flatted = lambda x: [j for i in x for j in i]
cycle_format = lambda permutation: "\nPermutation" + "".join(map(lambda x: str(tuple(x)), permutation.cyclic_form))



def get_matrix(matrix, show_type = False):
    if not show_type:
        return '\n' + '\n'.join('\t'.join('%0.3f' % x for x in y) for y in matrix)
    else:
        matrix = expand(Matrix(matrix))#
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

    def _display_calculated_elemets_representations(self, exact_output = False):
        t = 0
        not_calc = []
        self.output = repr(self)
        for elem_arr_form, perm_el in self:
            if not hasattr(perm_el, 'representation'):
                not_calc.append(perm_el)
                t += 1
            else:
                self.output += '\n' + cycle_format(perm_el) + get_matrix(perm_el.representation, exact_output)
        self.output += f"\n\n{t} elements not calculated!"
        self.output += f"\nNot calculated elements:\n{not_calc}\n" if t > 0 else '\n'
        self.output += '_'*50
        print(self.output)

class SymmetricGroup(FiniteGroup):
    def __init__(self, N: int):
        assert isinstance(N, int) and N > 0
        self.n = N
        self.elements = {}
        for x in permutations(range(1, self.n + 1)):
            el = Permutation((0,) + x, size = self.n + 1)
            self.elements[tuple(el._array_form)] = el

    def find_simple_transps_reprs(self):
        n = len(self.standart_tableauxes)
        conjugate = self.tableauxes.find_conjugate_transpositions( )
        axial_distances = self.tableauxes.find_axial_distances( )
        for x in range(1, self.n):  # going thru transpositions x, x+1
            m = np.zeros((n, n), Rational if self.show_exact_values else np.float)  # create matrix a filled zeroes
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
        if not self._shape.is_self_conj:
            super(AlternatingGroup, self).find_simple_transps_reprs()
            return
        self.simple_transpositions = {}
        n = len(self.standart_tableauxes)  # //2
        self.tableauxes.find_conjugate_transpositions( )
        #conjugate = self.tableauxes.restrict_conjugate_to_alternating( ) if self._n > 4 else self.tableauxes.conjugate_transpositions
        conjugate = self.tableauxes.calling_to_restrict()
        axial_distances = self.tableauxes.find_axial_distances( )
        #print(*conjugate.items( ), sep = "\n")
        for x in range(1, self.n):  # going thru transpositions x, x+1
            # a = zeros(n, n) #create matrix filled with zeroes
            m = np.zeros((n, n), Rational if self.show_exact_values else np.float)
            conj = conjugate[x]
            axl_dist = axial_distances.get(x)
            state = False
            for index in range(n // 2 if self.n > 4 else n):
                diag, non_diag = axl_dist.get(index)
                if m[index, index] == 0:
                    m[index, index] = -diag
                    m[n-index-1, n-index-1] = diag
                if index in conj:
                    if conj[index]:
                        block = conj.get(index)

                        m[block, block] = diag
                        m[index, block] = m[block, index] = non_diag

                        m[n - block - 1, n - block - 1] = -diag
                        m[n - index - 1, n - block - 1] = m[n - block - 1, n - index - 1] = non_diag
                    else:
                        block = conj.get(str(index))
                        state = True

                        m[block, block] = diag
                        m[index, block] = non_diag * self._shape.mult
                        m[block, index] = non_diag * self._shape.mult ** (-1)

                        m[n - block - 1, n - block - 1] = -diag
                        m[n - block - 1, n - index - 1] = non_diag * self._shape.mult
                        m[n - index - 1, n - block - 1] = non_diag * self._shape.mult ** (-1)

            if x == 1:
                self.simple_transpositions[Permutation(1, 2, size = self.n + 1)] = m#result
            else:
                if not state:
                    self.simple_transpositions[Permutation(x, x + 1, size = self.n + 1)] = m
                    self.elements[tuple((Permutation(1, 2, size = self.n + 1) * Permutation(x, x + 1, size = self.n + 1))._array_form)].representation = self.simple_transpositions[Permutation(1, 2, size = self.n + 1)].dot(m)
                else:
                    self.simple_transpositions[Permutation(x, x + 1, size = self.n + 1)] = m.dot(self.simple_transpositions[Permutation(1, 2, size = self.n + 1)])
                    self.elements[tuple((Permutation(1, 2, size = self.n + 1) * Permutation(x, x + 1, size = self.n + 1))._array_form)].representation = m


class Representation:
    show_exact_values = False

    def __init__(self, group_instance, get_exact_output = False, shape = None):
        self.group = group_instance
        self.show_exact_values = self.group.show_exact_values = get_exact_output
        if not shape:
            self._n = self.group.n
            _Partitions = Partitions(self._n, isinstance(self.group, AlternatingGroup))
            print(_Partitions)
            self._shape = self.group._shape = _Partitions[int(input("\nInput partition index: "))]
            self.group.tableauxes = Standart_tableauxes(self._shape)
            self.group.standart_tableauxes = self.group.tableauxes.get_standart( )
            print(self.group.tableauxes)
        elif shape:
            self._shape = Partition(shape) if isinstance(shape, (tuple, list)) else shape
        self.find_simple_transps_reprs( )

    def find_simple_transps_reprs(self):
        self.group.simple_transpositions = {}
        self.group.find_simple_transps_reprs()
        self.group.find_remain_simple_transps()

    def find_remain_elems(self):
        [self.__getitem__(elem) for elem in self.group.elements]
        return 'Function ends work'

    def find_remain_elems_parallel(self):
        from functools import partial
        from multiprocessing import Manager, Pool as ProcessPool

        manager = Manager( )
        self.d = manager.dict()
        pool = ProcessPool(processes = 3)
        _pool = pool.map(self.find_elem, self.group.elements)
        for el, repres in self.d.items():
            self.group.elements[el].representation = repres
        return 'Function ends work'


    def _display_calculated_elemets_representations(self):
        return self.group._display_calculated_elemets_representations(self.show_exact_values)

    def __getitem__(self, element_array_form):

        if not hasattr(self.group, 'simple_transpositions'):
            self.group.find_simple_transps_reprs()

        temp_perm = Permutation(element_array_form._array_form, size = self._n + 1) if isinstance(element_array_form, Permutation)\
            else Permutation(element_array_form, size = self._n + 1)
        if tuple(temp_perm._array_form) != element_array_form:
            element_array_form = tuple(temp_perm._array_form)

        assert element_array_form in self.group.elements, f"Element {self.group[element_array_form]} not in {repr(self.group)}"

        permutation_element = self.group[element_array_form]
        if hasattr(permutation_element, 'representation'):
            return permutation_element.representation

        if not hasattr(permutation_element, 'transp_factor'):
            permutation_element.transp_factor = [Permutation([x], size = permutation_element.size) for x in permutation_element.transpositions( )[::-1]]

        if len(permutation_element.transp_factor) > 1:
            permutation_element.representation = np.linalg.multi_dot([self.group.simple_transpositions[x] for x in permutation_element.transp_factor])
        elif len(permutation_element.transp_factor) == 0:
            permutation_element.representation = np.eye(len(self.group.standart_tableauxes))

        try:
            return permutation_element.representation
        except AttributeError:
            pass

    def find_elem(self, element_array_form):

        if not hasattr(self.group, 'simple_transpositions'):
            self.group.find_simple_transps_reprs( )

        temp_perm = Permutation(element_array_form._array_form, size = self._n + 1) if isinstance(
            element_array_form, Permutation) \
            else Permutation(element_array_form, size = self._n + 1)
        if tuple(temp_perm._array_form) != element_array_form:
            element_array_form = tuple(temp_perm._array_form)

        assert element_array_form in self.group.elements, f"Element {self.group[element_array_form]} not in {repr(self.group)}"

        permutation_element = self.group[element_array_form]
        if hasattr(permutation_element, 'representation'):
            self.d[tuple(permutation_element._array_form)] = permutation_element.representation
            #return permutation_element.representation

        if not hasattr(permutation_element, 'transp_factor'):
            permutation_element.transp_factor = [Permutation([x], size = permutation_element.size) for x in
                                                 permutation_element.transpositions( )[::-1]]

        if len(permutation_element.transp_factor) > 1:
            self.d[tuple(permutation_element._array_form)] = np.linalg.multi_dot([self.group.simple_transpositions[x] for x in permutation_element.transp_factor])
        elif len(permutation_element.transp_factor) == 0:
            self.d[tuple(permutation_element._array_form)] = np.eye(len(self.group.standart_tableauxes))


if __name__ == "__main__":
    S5 = SymmetricGroup(5)
    A5 = AlternatingGroup(5)
    repres = Representation(S5, get_exact_output = True)

    #print(repres.find_remain_elems( ))
    print(repres.find_remain_elems_parallel())
