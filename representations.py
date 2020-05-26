import tableauxes
import partitions

import itertools
import numpy as np

from sympy import Rational, expand, Matrix, pretty, I
from sympy.combinatorics import Permutation
from functools import reduce
from typing import Iterable, Tuple, Union, Generator
from numpy.linalg import multi_dot

from multiprocessing import Manager, Pool as ProcessPool


cycle_format = lambda permutation: "\nPermutation" + "".join(map(lambda x: str(tuple(x)), permutation.cyclic_form))
np.set_printoptions(precision = 3, suppress = True, linewidth = 100)

def get_matrix(matrix: np.ndarray, show_exact: bool = False) -> str:
    if not show_exact:
        try:
            return '\n' + pretty(Matrix((np.array(matrix, dtype = float).round(decimals = 3))), wrap_line = False, use_unicode = True)
        except:
            pass
    return '\n' + pretty(expand(Matrix(matrix)), wrap_line=False, use_unicode=True, full_prec=False)


class FiniteGroup:
    def __init__(self, N: int, shape: partitions.partition, exact_output: bool = False):
        self.n = N
        self._shape = shape
        self.show_exact_values = exact_output

        self.elements = {}
        self.simple_transpositions = {}


    def __call__(self, *args, **kwargs):
        pass

    def __iter__(self) -> Generator:
        yield from self.elements.items()

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.n})"

    def __getitem__(self, element: Permutation) -> np.ndarray:
        return self.elements[element]

    def _calculated_elemets(self) -> Generator:
        not_calculated_elems = []
        yield repr(self) + ' - ' + str(self._shape)
        for perm, repres in self:
            try:
                if not repres:
                    not_calculated_elems.append(perm)
            except:
                yield '\n' + cycle_format(perm) + get_matrix(repres, self.show_exact_values)
        yield f"\n\n{len(not_calculated_elems)} elements not calculated!\nNot calculated elements:\n"
        yield ' '.join(map(cycle_format, not_calculated_elems))


class SymmetricGroup(FiniteGroup):
    def __init__(self, N: int, shape: partitions.partition, exact_output: bool = False, **kwargs) -> None:
        assert isinstance(N, int) and N > 0
        super(SymmetricGroup, self).__init__(N, shape, exact_output)

        self._tableauxes = tableauxes.standart_tableauxes(self._shape)
        self.standart_tableauxes = self._tableauxes.standart_tableauxes
        self.matrix_size = len(self.standart_tableauxes)
        self.elements = dict.fromkeys(self._elements_gen(), None) #element: representation

    def _elements_gen(self) -> Generator:
        yield from (Permutation((0,) + x, size=self.n + 1) for x in itertools.permutations(range(1, self.n + 1)))

    def find_simple_transps_reprs(self) -> None:
        n = self.matrix_size
        conjugate = self._tableauxes.get_conjugated_tablues( )
        axial_distances = self._tableauxes.find_axial_distances( )
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

            el = Permutation(x, x + 1, size = self.n + 1)
            self.simple_transpositions[el] = m

            if not isinstance(self, AlternatingGroup):
                self.elements[el] = m
            elif x == 1:
                p12 = m
            else:
                self.elements[Permutation(1, 2, size=self.n + 1) * el] = p12.dot(m)

    def find_remain_simple_transps(self) -> None:
        def calculate_transposition(x: int, y: int) -> np.ndarray:
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
            for j in range(i + 2, self.n + 1):
                if not (Permutation(i, j, size = self.n + 1) in self.simple_transpositions):
                    m = calculate_transposition(i, j)
                    el = Permutation(i, j, size = self.n + 1)
                    self.simple_transpositions[el] = m
                    if el in self.elements:
                        self.elements[el] = m


class AlternatingGroup(SymmetricGroup):

    def _elements_gen(self) -> Generator:
        yield from (Permutation((0,) + x, size = self.n + 1) for x in itertools.permutations(range(1, self.n + 1))
                    if Permutation((0,) + x, size = self.n + 1).is_even)

    @property
    def multiplier(self) -> Tuple[float, bool]:
        temp = [2 * lambda_i - 2 * (i + 1) + 1 for i, lambda_i in enumerate(self._shape)]
        power = (reduce(lambda x,y: x*y if y > 0 else x, temp) - 1) // 2
        multiplier = I**power
        is_complex = power % 2 != 0
        return multiplier, is_complex

    def find_simple_transps_reprs(self) -> None:
        if self.n <= 4 or not self._shape.is_self_conjugated:
            super(AlternatingGroup, self).find_simple_transps_reprs()
            return

        n = len(self.standart_tableauxes)
        half = n//2
        conjugate = self._tableauxes.get_restricted_conjugated_tablues()
        axial_distances = self._tableauxes.find_axial_distances( )
        mult, is_complex = self.multiplier
        datatype = Rational if self.show_exact_values else np.complex if is_complex else np.float

        for x in range(1, self.n):
            m = np.zeros((n, n), datatype)
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
                        if index < half:
                            m[index, block] = non_diag * mult
                            m[block, index] = non_diag * mult ** (-1)
                        else:
                            m[index, block] = non_diag * mult ** (-1)
                            m[block, index] = non_diag * mult

            if x == 1:
                perm12 = Permutation(1, 2, size=self.n + 1)
                p12 = m
                self.simple_transpositions[perm12] = m
            else:
                perm = Permutation(x, x + 1, size = self.n + 1)
                if not state:
                    self.simple_transpositions[perm] = m
                    self.elements[perm12 * perm] = p12.dot(m)
                else:
                    self.simple_transpositions[perm] = m.dot(self.simple_transpositions[perm12])
                    self.elements[perm12 * perm] = m


class Representation:
    _all_elems = False
    def __init__(self, _group: FiniteGroup, n: int = 0, shape: Iterable = None, get_exact_output: bool = False) -> None:
        if shape:
            assert isinstance(shape, (list, tuple)), 'Shape must be List or Tuple!'
            self._shape = partitions.partition(partition = shape)
        elif n:
            _Partitions = partitions.partitions(n, issubclass(_group, AlternatingGroup))
            print(_Partitions)
            self._shape = _Partitions[int(input("\nInput partition index: "))]

        self._n = self._shape.n
        self.group = _group(N = self._n, shape = self._shape, exact_output = get_exact_output)

    def __getitem__(self, elem: Permutation) -> Union[None, np.ndarray]:
        '''Returns matrix represents element'''
        elem = Permutation(elem.array_form, size = self.group.n + 1)
        if not self.group.simple_transpositions:
            self.group.find_simple_transps_reprs( )

        assert elem in self.group.elements, f"Element {elem} not in {repr(self.group)}"
        self.compute_element(elem)
        return self.group.elements[elem]

    def compute_element(self, elem: Permutation) -> None:
        '''Finds element and modify group dict state'''
        try:
            self.group.elements.get(elem)
        except ValueError:
            return

        transp_factor = [Permutation([x], size = self.group.n + 1) for x in elem.transpositions( )[::-1]]
        if transp_factor:
            try:
                self.group.elements[elem] = np.linalg.multi_dot([self.group.simple_transpositions[x] for x in transp_factor])
            except ValueError:
                self.group.elements[elem] = self.group.simple_transpositions[elem]
        else:
            self.group.elements[elem] = np.eye(len(self.group.standart_tableauxes))

    def compute_element_parallel(self, elem: Permutation) -> None:
        '''Finds element and drop it into shared memory dict'''
        transp_factor = [Permutation([x], size = elem.size) for x in elem.transpositions( )[::-1]]
        if transp_factor:
            try:
                self._shared_dict[elem] = np.linalg.multi_dot([self.group.simple_transpositions[x] for x in transp_factor])
            except ValueError:
                self._shared_dict[elem] = self.group.simple_transpositions[elem]
        else:
            self._shared_dict[elem] = np.eye(len(self.group.standart_tableauxes))

    def find_simple_transps_reprs(self) -> None:
        '''Finds simple traspositions (k, k + i + 1)'''
        self.group.simple_transpositions = {}
        self.group.find_simple_transps_reprs()
        self.group.find_remain_simple_transps()

    def find_remain_elems(self, parallel_type: str = '') -> None:
        '''Finds all group elemets'''
        if not parallel_type:
            tuple(self.compute_element(elem) for elem in self.group.elements)
        else:
            if parallel_type == 'process' or parallel_type == 'p':
                self._shared_dict = Manager().dict()
                pool = ProcessPool(processes = 5)
            else:
                raise AttributeError(f'{parallel_type} in not supported!')
            _pool = pool.map(self.compute_element_parallel, self.group.elements)
            self.group.elements = {**self._shared_dict}
        self._all_elems = True

    def write_to_file(self, filename: str = '', console_output: bool = False) -> Tuple[str, str]:
        '''Writes to file finded elements'''
        if not filename:
            filename = repr(self.group) + '-' + str(self._shape)
            filename += ' --float' if not self.group.show_exact_values else ' --exact'
            filename += ' --all' if self._all_elems else ' --gens_only'
        output_string = self.output
        with open(filename + '.txt', 'w', encoding='utf8') as f:
            f.write(output_string)
        return output_string if console_output else '', filename

    @property
    def output(self):
        """Returns all calculated matricies in string format"""
        return ''.join(self.group._calculated_elemets())


if __name__ == "__main__":
    representation = Representation(AlternatingGroup, shape = (3, 1, 1), get_exact_output =True)
    representation.find_simple_transps_reprs( )
    #representation.find_remain_elems()
    print(representation.output)
    #representation.write_to_file()
    # m = representation[Permutation([(1,2), (7, 8)])]
    # print(get_matrix(m, True))


