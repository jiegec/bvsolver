from pwn import process
import bvsolver
import os
import tempfile

code = """
#include <random>
#include <iostream>
int main() {
    std::random_device rd;
    std::mt19937 gen(rd());
    for (int i = 0;i < 10000;i++) {
        std::cout << gen() << std::endl;
    }
    return 0;
}
"""

with tempfile.NamedTemporaryFile("w", suffix=".cpp") as f:
    f.write(code)
    f.flush()

    os.system(f"g++ -O2 {f.name} -o {f.name}.bin")

    p = process([f"{f.name}.bin"])
    all_nums = [int(p.recvline().decode().strip()) for i in range(10000)]
    known = all_nums[:624]

    # create equations
    solver = bvsolver.Solver([32] * 624)
    state = solver.bitvectors()

    # cpp mt19937 is same as cpython, except for its seeding
    rng = bvsolver.CPythonRandom(
        state,
        624,
    )

    zeros = []
    for value in known:
        zeros.append(rng.getrandbits(32) ^ value)

    # find solutions
    for sol in solver.solve(zeros):
        print("Found solution")
        state_recover = [s.get(sol) for s in state]

        # recreate random number generator
        rng_recover = bvsolver.CPythonRandom(state_recover, 624).to_cpython_random()
        assert all_nums == [rng_recover.getrandbits(32) for i in range(len(all_nums))]

        # multiple solutions, choose one
        break
