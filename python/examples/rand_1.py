import bvsolver
import random

# assume we got some random numbers from the remote
all_nums = [random.getrandbits(1) for i in range(30000)]
known = all_nums[:20000]

# create equations
solver = bvsolver.Solver([32] * 624)
state = solver.bitvectors()

rng = bvsolver.CPythonRandom(
    state,
    # newly created Random instance has its index equal to 624
    624,
)
# find solutions
for sol in solver.solve(
    [rng.getrandbits(1) ^ value for value in known]
    # newly created Random instance has its first element equal to 0x80000000
    # see init_by_array in _randommodule.c
    + [state[0] ^ 0x80000000]
):
    print("Found solution")
    state_recover = [s.get(sol) for s in state]

    # recreate random number generator
    rng_recover = bvsolver.CPythonRandom(state_recover, 624).to_cpython_random()
    for i in range(len(all_nums)):
        assert rng_recover.getrandbits(1) == all_nums[i]
