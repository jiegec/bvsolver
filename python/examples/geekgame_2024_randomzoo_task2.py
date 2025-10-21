from pwn import *
import bvsolver
import tempfile

# https://github.com/PKU-GeekGame/geekgame-4th/tree/master/official_writeup/algo-randomzoo
code = """
import random
flag = b"flag{fake_flag_for_local_testing}"
for i in range(2**64):
    print(random.getrandbits(32)+flag[i%len(flag)])
    input()
"""

with tempfile.NamedTemporaryFile("w", suffix=".py") as f:
    f.write(code)
    f.flush()

    # sometimes it overflows so no solution is found
    # loop until solved
    while True:
        p = process(["python3", f.name])

        witness = 2000
        known = []
        for i in range(witness):
            known.append(int(p.recvline().decode().strip()))
            p.sendline()

        # create equations
        solver = bvsolver.Solver([32] * 624)
        state = solver.bitvectors()

        rng = bvsolver.CPythonRandom(
            state,
            # newly created Random instance has its index equal to 624
            624,
        )
        # drop the lowest 20 bits, assuming adding flag does not overflow to upper bits
        shift = 20
        # find solutions
        for sol in solver.solve(
            [(rng.getrandbits(32) >> shift) ^ (value >> shift) for value in known]
            # newly created Random instance has its first element equal to 0x80000000
            # see init_by_array in _randommodule.c
            + [state[0] ^ 0x80000000]
        ):
            print("Found solution")
            state_recover = [s.get(sol) for s in state]

            # recreate random number generator
            rng_recover = bvsolver.CPythonRandom(state_recover, 624).to_cpython_random()
            flag = ""
            for i in range(50):
                flag += chr(known[i] - rng_recover.getrandbits(32))
            print(flag)

            exit(0)
