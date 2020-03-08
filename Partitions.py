from functools import reduce
from sympy import I

def my_cache(function_to_decorate):
    function_cache = {}
    def wrapper(key):
        if function_cache.get(key):
            return function_cache.get(key)
        else:
            function_cache[key] = function_to_decorate(key)
            return function_cache[key]
    return wrapper

@my_cache
def quantity(n):
    if n == 0:
        return 1
    t = 0
    for k in range(1, n + 1):
        q1 = quantity(n - (3 * k * k - k) // 2)
        q2 = quantity(n - (3 * k * k + k) // 2)
        t += (-1) ** (k + 1) * (q1 + q2)
    return t

class Partition:
    def __init__(self, partition, find_conjugate = False):
        assert isinstance(partition, (tuple, list)), "Partition must be tuple or list type!"
        self.sequence = tuple(partition)
        self.n = sum(self.sequence)
        if find_conjugate:
            self.is_self_conj = self.is_self_conjugated()

    def __iter__(self):
        for i in self.sequence:
            yield i

    def __str__(self):
        return ", ".join(map(str, self))

    # def __len__(self):
    #     return len(self.sequence)

    def is_self_conjugated(self):
        '''Returns statement partition is self-conjugated'''
        p = list(self.sequence[::])
        conjugated = []
        index = 0
        while p:
            level = len(p)
            if self.sequence[index] != level:
                return False
            conjugated.append(level)
            p = [x - 1 for x in p if x - 1 > 0]
            index += 1
        return True

    def conjugate_partition(self):
        '''Transposing the Young diagram about its main diagonal.\nExamples:\n[5,5,3,1] --> [4,3,3,2,2]\n[6,5,3,1] --> [4,3,3,2,2,1]\n[2,1]     --> [2,1]\nAnd returns if diagram is self-conjugated'''
        if not self.conjugated:
            p = list(self.sequence[::])
            conjugate = []
            while p:
                level = len(p)
                conjugate.append(level)
                p = [x - 1 for x in p if x - 1 > 0]
            self.conjugated = Partition(conjugate)
        return self.conjugated

    def multiplier(self):
        if not hasattr(self, 'mult'):
            temp = [2 * lambda_i - 2 * (i + 1) + 1 for i, lambda_i in enumerate(self)]
            power = reduce(lambda x,y: x*y if y > 0 else x, temp)
            power = (power - 1)//2
            self.mult = I**power
            self.is_complex = power % 2 != 0
        return self.mult, self.is_complex

    def template(self):
        if not hasattr(self, 'tempt'):
            index = 0
            self.tempt = {}
            n = sum(self.sequence)
            for i, step in enumerate(self):
                max = index + step - 1
                max_next = max + self.sequence[i + 1] if i < len(self.sequence) - 1 else 0
                for j in range(index, index + step):
                    if j < max:
                        self.tempt[j] = j + 1
                    if j + step < n and j + step <= max_next:
                        self.tempt[str(j)] = j + step
                index += step
        return self.tempt

    def __eq__(self, other):
        return self.sequence == other.sequence

class Partitions:
    '''Find a partitions of integer.\nExamples:\n 2 --> [[2], [1,1]]\n4 --> [[4], [3,1], [2,2], [2,1,1], [1,1,1,1]]\n'''
    partitions = {}
    def __init__(self, n: int, isselfconjugated = False):
        assert  isinstance(n, int) and n > 0, "n must be natural number!"
        self._n = n
        self.check = isselfconjugated

    def quantity(self):
        '''Returns exact quantity of partitions of integer n'''
        return quantity(self._n)

    def approximate_quantity(self):
        '''Returns approximate quantity of partitions of integer n'''
        from math import pi, exp, sqrt, pow
        from decimal import Decimal, getcontext
        getcontext().prec = 45
        c = Decimal(self._n) - Decimal(1)/Decimal(24)
        n1 = Decimal(pi)/(Decimal(sqrt(6))*c) - Decimal(1)/Decimal(2* pow(c, 1.5))
        n2 = Decimal(pi*sqrt(2*c/3)).exp()
        d = Decimal(2*pi*sqrt(2))
        return int(n1*n2/d)
        #return float(exp(pi*sqrt(2*self._n/3))/(4*self._n*sqrt(3)))

    def __getitem__(self, index):
        return self.partitions.get(index)

    def get_dict(self):
        l = list(self)
        n = len(l)
        self.partitions = dict(zip(range(n), map(lambda x: x.list, l)))
        return self.partitions

    def __iter__(self):
        def gen(remaining, temp = tuple( )):
            if remaining == 0:
                yield Partition(temp, self.check)
            elif remaining > 0:
                for step in range(remaining if not temp else temp[-1], 0, -1):
                    yield from gen(remaining - step, temp + (step,))
        return gen(self._n)

    def show(self):
        print(f"№\tPartitions\tSplitable" if self.check else f"№\tPartitions")
        for index, Item in enumerate(self):
            print(f"{index}\t{Item}\t{Item.is_self_conj}" if self.check else f"{index}\t{Item}")
            self.partitions[index] = Item


    def __str__(self):
        from texttable import Texttable
        table = Texttable(0)
        table.set_deco(Texttable.HEADER)
        table.add_row(["№", "Partitions", "Splitable"] if self.check else ["№", "Partitions"])
        for index, Item in enumerate(self):
            table.add_row([index, Item, Item.is_self_conj] if self.check else [index, Item])
            self.partitions[index] = Item
        return table.draw()

if __name__ == "__main__":
    from time import time
    n = int(input("Input integer number: "))
    print(Partitions(n, True))
    # for i in range(1, 100):
    #     t = time()
    #     p = Partitions(i)
    #     #print(p)
    #     q1, q2 = quantity(i), p.approximate_quantity()
    #     print(i, q1, q2, abs(q2 - q1))#, abs(q2 - q1)/q1)
    #input("End...")



"""
def exp_form(self):
    '''[6, 4, 4, 2, 1] -> [1, 0, 2, 0, 1, 1]
'''
    assert len(self.sequence) > 0
    self.exp_form = [self.sequence.count(x) for x in range(self.sequence[0], 0, -1)]
    return self.exp_form

def log_form(self, exp_form):
    '''[1, 0, 2, 0, 1, 1] -> [6, 4, 4, 2, 1]
'''
    assert isinstance(exp_form, (tuple, list))
    n = len(exp_form)
    partition = []
    for sub, element in enumerate(exp_form):
        partition.extend([n-sub]*element)
    return Partition(partition)
"""