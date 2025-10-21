from pwn import process
import bvsolver
import os
import tempfile

code = """
#include <random>
#include <iostream>
int main() {
    std::random_device rd;
    auto seed = rd();
    std::cout << seed << std::endl;
    std::mt19937_64 gen(seed);
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
    seed = int(p.recvline().decode().strip())
    all_nums = [int(p.recvline().decode().strip()) for i in range(10000)]
    # verify our seeder
    test_rng = bvsolver.MT19937_64([0] * 312, 312)
    test_rng.seed(seed)
    for i in range(10000):
        gen = test_rng.gen()
        assert all_nums[i] == gen

    known = all_nums[:312]

    # create equations
    solver = bvsolver.Solver([64] * 312)
    state = solver.bitvectors()

    # cpp mt19937_64
    rng = bvsolver.MT19937_64(
        state,
        312,
    )

    zeros = []
    for value in known:
        zeros.append(rng.gen() ^ value)

    # find solutions
    for sol in solver.solve(zeros):
        print("Found solution")
        state_recover = [s.get(sol) for s in state]

        # recreate random number generator
        rng_recover = bvsolver.MT19937_64(state_recover, 312)
        assert all_nums == [rng_recover.gen() for i in range(len(all_nums))]

        # multiple solutions, choose one
        break
