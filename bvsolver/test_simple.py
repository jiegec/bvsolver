import bvsolver


def test_simple():
    solver = bvsolver.Solver([1] * 4)
    a, b, c, d = solver.bitvectors()

    zeros: list[bvsolver.BitVector] = []
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
        assert a.get(solution) == 1
        assert b.get(solution) == 0
        assert c.get(solution) == 1
        assert d.get(solution) == 0


def test_multibit():
    solver = bvsolver.Solver([8])
    (a,) = solver.bitvectors()

    # a ^ 0x34 == 0
    solutions = list(solver.solve([a ^ 0x34]))
    assert len(solutions) == 1
    for solution in solutions:
        assert a.get(solution) == 0x34


def test_multi_solutions():
    solver = bvsolver.Solver([1] * 4)
    a, b, c, d = solver.bitvectors()

    zeros: list[bvsolver.BitVector] = []
    # a ^ 1 == 0
    zeros.append(a ^ 1)
    # a ^ b ^ 1 == 0
    zeros.append(a ^ b ^ 1)
    # a ^ c == 0
    zeros.append(a ^ c)
    # We can solve: [a, b, c, d] == [1, 0, 1, 0 or 1]
    solutions = list(solver.solve(zeros))
    assert len(solutions) == 2
    for solution in solutions:
        assert a.get(solution) == 1
        assert b.get(solution) == 0
        assert c.get(solution) == 1

    # no constraints
    solutions = list(solver.solve([]))
    assert len(solutions) == 16
    for solution in solutions:
        print(solution)


def test_no_solutions():
    solver = bvsolver.Solver([1] * 2)
    a, b = solver.bitvectors()

    zeros: list[bvsolver.BitVector] = []
    # a ^ 1 == 0
    zeros.append(a ^ 1)
    # a == 0
    zeros.append(a)
    solutions = list(solver.solve(zeros))
    assert len(solutions) == 0
