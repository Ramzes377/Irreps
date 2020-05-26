from partitions import partition
from sympy import sqrt, Rational
from typing import Iterable, List, TypeVar, Tuple, Dict, Generator, Union

T = TypeVar('T', bound='tablue')

class tablue:
    """Represents Young's tableau which can be generate by sequence of integers and shape.

    Examples:

    print(Tablue(tablue=(3, 1, 2, 4, 5), shape = partition(3, 2)))
    3| 1| 2|
    4| 5

    print(Tablue(tablue=range(1, 10), shape = partition(4, 3, 2)))
    1| 2| 3| 4|
    5| 6| 7|
    8| 9
    """
    def __init__(self, *, shape: partition, tablue: Iterable) -> None: #taking only keyword arguments
        assert isinstance(shape, partition), 'Shape must be Partition'
        self.shape = shape
        self.tablue = tuple(tablue)
        self.n = self.shape.n
        assert len(self.tablue) == self.n, "Tablue must have same number of symbols with sum(shape)"
        self.conjugate = {}

    def __getitem__(self, item: Union[int, slice]) -> Union[int, list]:
        return self.tablue[item]

    def __str__(self) -> str:
        return "|\n".join( "| ".join(str(number) for number in string) for string in self.get_mll) + "\n"

    def __eq__(self, other: T) -> bool:
        return self.shape == other.shape and tuple(self.tablue) == tuple(other.tablue)

    @property
    def get_mll(self) -> List[List]:
        if not hasattr(self, '_tablue_as_mll'):
            get_slice = lambda i: slice(sum(self.shape[:i]), sum(self.shape[:i + 1]))
            self._tablue_as_mll = [list(self.tablue[get_slice(i)]) for i, lenght in enumerate(self.shape)]
        return self._tablue_as_mll

    @property
    def is_standart(self) -> bool:
        if not hasattr(self, '_is_standart'):
            t = self.shape.template
            for i in range(self.n):
                case1, case2 = t.get(i), t.get(str(i))
                if (case1 and self[i] > self[case1]) or (case2 and self[i] > self[case2]):
                    self._is_standart = False
                    break
            else:
                self._is_standart = True
        return self._is_standart

    @property
    def grid(self) -> Dict[int, Tuple[int, int]]:
        if not hasattr(self, '_grid'):
            self._grid = {}
            mll = self.get_mll
            for i in range(len(mll)):
                level = mll[i]
                for j in range(len(level)):
                    self._grid[level[j]] = (i, j)
        return self._grid

    def restrict(self, restriction_length: int = 2) -> T:
        def reflect_and_concat(grid, remainder):
            new_grid = {}
            for x, pos in grid.items():
                new_grid[x] = (pos[1], pos[0])
            for x, pos in remainder:
                new_grid[x] = pos
            return new_grid
        grid = self.grid
        remainder = [(i, grid[i]) for i in range(self.n - restriction_length + 1, self.n + 1)]
        reflected_grid = reflect_and_concat(grid, remainder)
        seq = [key for key, value in sorted(sorted(reflected_grid.items( ), key = lambda x: x[0]), key = lambda x: x[1])]
        new_tablue = tablue(shape = self.shape, tablue = seq)
        return new_tablue

    def conjugated(self, x: int, y: int = 0) -> T:
        if not self.conjugate.get(x):
            y = x + 1 if not y else y
            _tablue = list(self.tablue)
            x_index, y_index = _tablue.index(x), _tablue.index(y)
            _tablue[x_index], _tablue[y_index] = y, x
            t = tablue(shape = self.shape, tablue = _tablue)
            self.conjugate[x] = t
        return self.conjugate.get(x)

    def axial_distance(self, x: int, y: int) -> int:
        assert 1 <= x <= self.n and 1 <= y <= self.n
        tablue = self.get_mll
        x_row = y_row = x_column = y_column = None
        for index, row in enumerate(tablue):
            if not x_row is None and not y_row  is None:
                break
            if x in row:
                x_row, x_column = index, row.index(x)
            if y in row:
                y_row, y_column = index, row.index(y)
        return x_column + y_row - x_row - y_column

    def get_block_values(self, x: int, y: int = 0) -> Tuple[float, float]:
        y = x + 1 if not y else y
        axl_dist = self.axial_distance(x, y)
        diagonal = Rational(1, axl_dist)
        non_diagonal = sqrt(1 - Rational(1, axl_dist * axl_dist))
        return diagonal, non_diagonal


class standart_tableauxes:
    """Helps to find all oriented standart Young's tableauexes for given partition.

    Examples:

    print(Standart_tableauxes(shape = (2, 1)))
    2 standart Young's frames:
    1| 2|
    3

    1| 3|
    2

    print(Standart_tableauxes(shape = (3, 1, 1)))
    6 standart Young's frames:
    1| 2| 3|
    4|
    5

    1| 2| 4|
    3|
    5

    1| 3| 4|
    2|
    5

    1| 2| 5|
    3|
    4

    1| 3| 5|
    2|
    4

    1| 4| 5|
    2|
    3
    """
    def __init__(self, shape: Iterable) -> None:
        assert isinstance(shape, (partition, list, tuple))
        self._shape = shape if isinstance(shape, partition) else partition(partition = shape)
        self._n = self._shape.n
        self.conjugated_tablues = {}

    def __iter__(self) -> Generator:
        assert self.standart_tableauxes, 'List of standart tableauxes is empty!'
        yield from self.standart_tableauxes

    def __str__(self) -> str:
        return f"\n{len(self.standart_tableauxes)} standart Young's frames:\n" + "\n".join(map(str, self))

    @property
    def standart_tableauxes(self) -> List[T]:
        '''Find standart Young tableauxes for given partition\n'''
        if not hasattr(self, '_standart_tableauxes'):
            t = tablue(shape =  self._shape, tablue = range(1, self._n + 1))
            elems = [t.get_mll.index(row) for index in range(1, self._n + 1) for row in t.get_mll if index in row]
            self._standart_tableauxes = (self.t_from_seq(t) for t in self._helper_unique_permutations(elems) if self.t_from_seq(t).is_standart)
            for i in range(1, self._n):
                self._standart_tableauxes = sorted(self._standart_tableauxes, key = lambda f: f[i], reverse = True)
        return self._standart_tableauxes

    def t_from_seq(self, seq: tuple) -> T:
        yamanuchi = dict(sorted(zip(range(1, len(seq) + 1), seq), key = lambda x: x[1]))
        t = tablue(shape = self._shape, tablue = yamanuchi.keys( ))
        t.yamanuchi = yamanuchi
        return t

    def get_conjugated_tablues(self) -> Dict[int, Dict[int, int]]:
        self.conjugated_tablues[1] = {}
        for n in range(2, self._n):
            temp = {}
            for i, t in enumerate(self.standart_tableauxes):
                try:
                    idx = self.standart_tableauxes[:i].index(t.conjugated(n))
                    temp[idx] = i
                except ValueError:
                    continue
            self.conjugated_tablues[n] = temp
        return self.conjugated_tablues

    def find_axial_distances(self) -> Dict[int, Dict[int, Tuple[float, float]]]:
        result = {}
        for x in range(1, self._n):
            result[x] = {t: self.standart_tableauxes[t].get_block_values(x) for t in range(len(self.standart_tableauxes))}
        return result

    def get_restricted_conjugated_tablues(self) -> Dict[int, Dict[int, Union[None, int]]]:

        self.get_conjugated_tablues()

        def get_restricred_index(tablue):
            restricted_and_flipped = tablue.restrict(2)
            try:
                index = self.standart_tableauxes.index(restricted_and_flipped)
            except ValueError:
                restricted_and_flipped = tablue.restrict(1)
                try:
                    index = self.standart_tableauxes.index(restricted_and_flipped)
                except ValueError:
                    index = None
            return index

        n = len(self.standart_tableauxes)//2
        result = {}
        for idx in range(1, self._n - 2):
            result[idx] = {x: y for x, y in self.conjugated_tablues[idx].items( )}

        for idx in range(self._n - 2, self._n):
            temp = {}
            for x, y in self.conjugated_tablues[idx].items( ):
                if x < n:
                    if y < n:
                        temp[x] = y
                    else:
                        tablue = self.standart_tableauxes[x]
                        index = get_restricred_index(tablue)
                        if not index is None and index < n and index > x:
                            temp[x] = None
                            temp[str(x)] = index

                        oposite_tablue = self.standart_tableauxes[y]
                        index = get_restricred_index(oposite_tablue)
                        if not index is None and index >= n and index < y:
                            temp[y] = None
                            temp[str(y)] = index
                elif y > n:
                    temp[x] = y
            if not temp: #maybe wrong sorted basis tablues
                for x, y in self.conjugated_tablues[idx].items():
                    if x < n:
                        if y < n:
                            temp[x] = y
                        else:
                            tablue = self.standart_tableauxes[x]
                            index = get_restricred_index(tablue)
                            if not index is None and index > x:
                                temp[x] = None
                                temp[str(x)] = index

                            oposite_tablue = self.standart_tableauxes[y]
                            index = get_restricred_index(oposite_tablue)
                            if not index is None and index < y:
                                temp[y] = None
                                temp[str(y)] = index
                    elif y > n:
                        temp[x] = y

            result[idx] = temp
        self.conjugated_tablues = result
        return self.conjugated_tablues

    def _helper_unique_permutations(self, elements: tuple) -> Generator:
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
        yield from unique_permutations(elements)


if __name__ == "__main__":
    y = standart_tableauxes(shape = (3, 3, 3))
    print(y)
    print(y.get_conjugated_tablues())
    print(y.get_restricted_conjugated_tablues())



''' 




'''