from pwn import process
import bvsolver
import os
import tempfile

code = """
for i=1,10000 do print(math.random(0)) end
"""

with tempfile.NamedTemporaryFile("w", suffix=".lua") as f:
    f.write(code)
    f.flush()

    p = process(["lua", f"{f.name}"])
    all_nums = [
        # convert signed to unsigned
        int(p.recvline().decode().strip()) & 0xFFFFFFFFFFFFFFFF
        for i in range(10000)
    ]
    known = all_nums[:4]

    # create equations
    solver = bvsolver.Solver([64] * 4)
    state = solver.bitvectors()

    # lua uses xoshiro256**
    rng = bvsolver.Xoshiro256ss(
        state,
    )

    zeros = []
    for value in known:
        zeros.append(rng.gen() ^ bvsolver.Xoshiro256ss.backward(value))

    # find solutions
    for sol in solver.solve(zeros):
        print("Found solution")
        state_recover = [s.get(sol) for s in state]

        # recreate random number generator
        rng_recover = bvsolver.Xoshiro256ss(state_recover)
        recovered = [
            bvsolver.Xoshiro256ss.forward(rng_recover.gen())
            for i in range(len(all_nums))
        ]
        assert all_nums == recovered
