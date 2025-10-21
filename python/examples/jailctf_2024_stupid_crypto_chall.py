import bvsolver
from random import random, seed
from ast import literal_eval

# https://github.com/jailctf/challenges-2024/tree/master/stupid-crypto-chall

"""
a, b = input('> ').split(" ")
if len(b) > 3346:
    print('nuh uh')
    exit()
seed(literal_eval(b))
total_prog = "".join([("".join(map(chr, range(32, 127))) + '\n')[int(random() * 96)] for i in range(literal_eval(a))])
eval(total_prog)
"""

code = "breakpoint()"
alphabet = "".join(map(chr, range(32, 127))) + "\n"
assert len(alphabet) == 96
target = [alphabet.index(ch) for ch in code]
print(target)

# recover random.random() and corresponding a & b
# create equations
solver = bvsolver.Solver([32] * 624)
state = solver.bitvectors()

rng = bvsolver.CPythonRandom(
    state,
    # newly created Random instance has its index equal to 624
    624,
)

# newly created Random instance has its first element equal to 0x80000000
# see init_by_array in _randommodule.c
zeros = [state[0] ^ 0x80000000]
for value in target:
    # see _random_Random_random_impl in _randommodule.c
    s = int((value + 0.5) / 96 * 9007199254740992.0)
    a = s // 67108864
    b = s % 67108864
    zeros.append((rng.getrandbits(32) >> 5) ^ a)
    zeros.append((rng.getrandbits(32) >> 6) ^ b)

# find solutions
for sol in solver.solve(zeros):
    print("Found solution")
    state_recover = [s.get(sol) for s in state]

    # recreate random number generator
    rng_recover = bvsolver.CPythonRandom(state_recover, 624).to_cpython_random()
    # verify
    for i in range(len(target)):
        r = rng_recover.random()
        assert int(r * 96) == target[i]

    # run code
    rng_seed = bvsolver.CPythonRandom.recover_seed(
        bvsolver.CPythonRandom(state_recover, 624).to_cpython_random()
    )
    seed(rng_seed)
    total_prog = "".join(
        [
            ("".join(map(chr, range(32, 127))) + "\n")[int(random() * 96)]
            for i in range(len(target))
        ]
    )
    print(total_prog)
    assert total_prog == "breakpoint()"

    break
