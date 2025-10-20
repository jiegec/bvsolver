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
