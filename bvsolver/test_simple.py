import secrets
from bvsolver import *


def test_simple() -> None:
    solver = Solver([1] * 4)
    a, b, c, d = solver.bitvectors()

    zeros: list[BitVector] = []
    # a ^ 1 == 0
    zeros.append(a ^ 1)
    # a ^ b ^ 1 == 0
    zeros.append(a ^ b ^ 1)
    # a ^ c == 0
    zeros.append(a ^ c)
    # a ^ d ^ 1 == 0
    zeros.append(a ^ d ^ 1)
    # We can solve: [a, b, c, d] == [1, 0, 1, 0]
    solutions = list(solver.solve(zeros))
    assert len(solutions) == 1
    for solution in solutions:
        assert (a ^ 1).get(solution) == 0
        assert a.get(solution) == 1
        assert b.get(solution) == 0
        assert c.get(solution) == 1
        assert d.get(solution) == 0

    # We can solve: [a, b, c, d] == [1, 0, 1, 0]
    # a ^ b ^ c ^ d == 0
    zeros.append(a ^ b ^ c ^ d)
    # c ^ d ^ 1 == 0
    zeros.append(c ^ d ^ 1)
    # d == 0
    zeros.append(d)
    solutions = list(solver.solve(zeros))
    assert len(solutions) == 1
    for solution in solutions:
        assert a.get(solution) == 1
        assert b.get(solution) == 0
        assert c.get(solution) == 1
        assert d.get(solution) == 0


def test_multibit() -> None:
    solver = Solver([8])
    (a,) = solver.bitvectors()

    # a ^ 0x34 == 0
    solutions = list(solver.solve([a ^ 0x34]))
    assert len(solutions) == 1
    for solution in solutions:
        assert a.get(solution) == 0x34


def test_multi_solutions() -> None:
    solver = Solver([1] * 4)
    a, b, c, d = solver.bitvectors()

    zeros: list[BitVector] = []
    # a ^ 1 == 0
    zeros.append(a ^ 1)
    # a ^ b ^ 1 == 0
    zeros.append(a ^ b ^ 1)
    # a ^ c == 0
    zeros.append(a ^ c)
    # We can solve: [a, b, c, d] == [1, 0, 1, 0 or 1]
    solutions = list(solver.solve(zeros))
    assert len(solutions) == 2
    last = None
    for solution in solutions:
        print(solution)
        assert a.get(solution) == 1
        assert b.get(solution) == 0
        assert c.get(solution) == 1
        if last is None:
            last = d.get(solution)
        else:
            assert d.get(solution) == 1 - last

    zeros = []
    # a ^ 1 == 0
    zeros.append(a ^ 1)
    # a ^ b ^ 1 == 0
    zeros.append(a ^ b ^ 1)
    # c ^ d == 0
    zeros.append(c ^ d)
    # We can solve: [a, b, c, d] == [1, 0, 0, 0] or [1, 0, 1, 1]
    solutions = list(solver.solve(zeros))
    assert len(solutions) == 2
    for solution in solutions:
        assert a.get(solution) == 1
        assert b.get(solution) == 0
        assert c.get(solution) == d.get(solution)

    # no constraints
    solutions = list(solver.solve([]))
    assert len(solutions) == 16
    assert len(set(solutions)) == 16  # all unique
    for solution in solutions:
        print(solution)


def test_no_solutions() -> None:
    solver = Solver([1] * 2)
    a, b = solver.bitvectors()

    zeros: list[BitVector] = []
    # a ^ 1 == 0
    zeros.append(a ^ 1)
    # a == 0
    zeros.append(a)
    solutions = list(solver.solve(zeros))
    assert len(solutions) == 0


def test_lfsr() -> None:
    width = 32
    initial_state = secrets.randbits(width)
    lfsr = FibonacciLFSR(width, secrets.randbits(width), initial_state)
    known = [lfsr.gen() for i in range(1024)]

    # solve
    solver = Solver([lfsr._width])
    (state,) = solver.bitvectors()
    lfsr_recover = FibonacciLFSR(lfsr._width, lfsr._poly, state)
    computed: list[BitVector] = [lfsr_recover.gen() for i in range(len(known))]
    solutions = list(
        solver.solve([actual ^ expected for expected, actual in zip(known, computed)])
    )
    assert len(solutions) == 1
    state_recover = state.get(solutions[0])
    assert state_recover == initial_state


def test_equiv() -> None:
    rng = random.Random(0)
    rng_ours = CPythonRandom.from_cpython_random(rng)
    assert [rng.getrandbits(32) for i in range(1000)] == [
        rng_ours.getrandbits(32) for i in range(1000)
    ]


def inner(bits, all_count, known_count) -> None:
    rng = random.Random(0)
    init_state = rng.getstate()
    all = [rng.getrandbits(bits) for i in range(all_count)]
    known = all[:known_count]

    # solve
    solver = Solver([32] * N)
    state = solver.bitvectors()
    rand_recover = CPythonRandom(
        state, 624
    )  # newly created Random instance has its index equal to 624
    computed: list[BitVector] = [
        rand_recover.getrandbits(bits) for i in range(len(known))
    ]
    solutions = list(
        solver.solve(
            [actual ^ expected for expected, actual in zip(known, computed)]
            # newly created Random instance has its first element equal to 0x80000000
            # see init_by_array in _randommodule.c
            + [state[0] ^ 0x80000000]
        )
    )
    assert len(solutions) == 1

    # verify
    state_recover = [s.get(solutions[0]) for s in state]
    assert known == [actual.get(solutions[0]) for actual in computed]
    assert list(init_state[1][:-1]) == state_recover
    rng = CPythonRandom(state_recover, 624)
    assert [rng.getrandbits(bits) for i in range(all_count)] == all
    rng = CPythonRandom(state_recover, 624).to_cpython_random()
    assert [rng.getrandbits(bits) for i in range(all_count)] == all


def test_random_32() -> None:
    inner(32, 1000, 624)


def test_random_16() -> None:
    inner(16, 2000, 1500)


def test_random_1() -> None:
    inner(1, 20000, 15000)
