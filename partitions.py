from typing import TypeVar, Dict, Union, Generator

P = TypeVar('P', bound='partition')

class partition:
    """Represents partition of integer number also it calls Young diagram
    p = partition(3, 2, 1)
    print(p)
    print(repr(p))

    3, 2, 1
    partition(3, 2, 1)

    #Use with generator
    print(partition(*range(5, 1, -1)))
    5, 4, 3, 2

    #Use generator without unpacking it
    print(partition(partition = range(5, 1, -1)))
    5, 4, 3, 2
    """
    def __init__(self, *elems, **kwargs) -> None:
        """Partition can be generate with sequence of integers"""
        self.sequence = tuple(kwargs['partition']) if kwargs.get('partition') else elems

    def __iter__(self) -> Generator:
        """Supports iterate by partition"""
        yield from self.sequence

    def __str__(self) -> str:
        return ", ".join(map(str, self))

    def __repr__(self) -> str:
        return 'partition(' + self.__str__() + ')'

    def __getitem__(self, item: Union[int, slice]) -> Union[int, list]:
        return self.sequence[item]

    def __eq__(self, other: P) -> bool:
        return self.sequence == other.sequence

    @property
    def n(self) -> int:
        if not hasattr(self, '_n'):
            self._n = sum(self.sequence)
        return self._n

    @property
    def is_self_conjugated(self) -> bool:
        '''Returns statement partition is self-conjugated.'''
        if not hasattr(self, '_is_self_conjugated'):
            def conjugated_helper(p, depth=0):
                if not p: #if p is empty -> shape self conjugated
                    return True
                elif self.sequence[depth] != len(p): #elif shape[depth] != len(p) -> shape not self conjugated
                    return False
                return conjugated_helper([x - 1 for x in p if x > 1], depth + 1) #else go recursive| cut one layer of p and again call function
            self._is_self_conjugated = conjugated_helper(self.sequence)
        return self._is_self_conjugated

    @property
    def template(self) -> Dict[Union[int, str], int]:
        """Specific template of partition. Helping calcualting does partition is standart"""
        if not hasattr(self, '_template'):
            index = 0
            self._template = {}
            n = sum(self.sequence)
            for i, step in enumerate(self):
                max = index + step - 1
                max_next = max + self.sequence[i + 1] if i < len(self.sequence) - 1 else 0
                for j in range(index, index + step):
                    if j < max:
                        self._template[j] = j + 1
                    if j + step < n and j + step <= max_next:
                        self._template[str(j)] = j + step
                index += step
        return self._template


class partitions:
    '''Find a partitions of natural integer.

    Examples:

    print(*list(Partitions(2)))
    2
    1, 1

    print(*list(Partitions(4)))
    4
    3, 1
    2, 2
    2, 1, 1
    1, 1, 1, 1
    '''

    def __init__(self, n: int, isselfconjugated: bool = False) -> None:
        assert  isinstance(n, int) and n > 0, "n must be natural number!"
        self._n = n
        self.check = isselfconjugated

    def __getitem__(self, index: Union[int, slice]) -> Union[int, list]:
        if not self.partitions:
            for i, partition in enumerate(self):
                if i == index:
                    return partition
        return self.partitions.get(index)

    def __iter__(self) -> Generator:
        def gen(remaining, temp = tuple( )):
            if remaining == 0:
                yield partition(partition = temp)
            elif remaining > 0:
                for step in range(remaining if not temp else temp[-1], 0, -1):
                    yield from gen(remaining - step, temp + (step,))
        yield from gen(self._n)

    def __str__(self) -> str:
        self.partitions = {}
        part_max_len = self._n * 3 + 1
        def handler(item):
            self.partitions[item[0]] = item[1]
            if self.check:
                return f'{item[0]:<6}{str(item[1]):<{part_max_len}}{item[1].is_self_conjugated}'
            return f'{item[0]:<6}{str(item[1]):>{part_max_len}}'
        result = f'{"№":<5}{"Partition"}'
        result += f'{"Splitable":>{part_max_len}}\n' if self.check else '\n'
        result += '\n'.join(map(handler, enumerate(self)))
        return result



if __name__ == "__main__":
    n = int(input("Input integer number: "))
    print(partitions(n, True))
