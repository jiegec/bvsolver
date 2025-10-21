from pwn import process
import bvsolver
import tempfile

# https://github.com/PKU-GeekGame/geekgame-3rd/blob/279868d814fd8371075e06cf705730c54a3e2a13/official_writeup/prob08-cookie/attachment/prob08-server.py
code = """
from random import Random
import secrets

the_void = Random(secrets.randbits(256))
smol_cookie = b"flag{fake_flag_for_local_testing}"
words = b'\\0' * 2500 + smol_cookie

def xor_arrays(a, b, *args):
    if args:
        return xor_arrays(a, xor_arrays(b, *args))
    return bytes([x ^ y for x, y in zip(a, b)])

ancient_words = xor_arrays(words, the_void.randbytes(len(words)))
print(ancient_words.hex())
"""

with tempfile.NamedTemporaryFile("w", suffix=".py") as f:
    f.write(code)
    f.flush()

    p = process(["python3", f.name])
    d = bytes.fromhex(p.recvline().decode().strip())

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
    for i in range(2500 // 4):
        gen = rng.getrandbits(32)
        # little endian
        zeros.append((gen & 0xFF) ^ d[i * 4 + 0])
        zeros.append(((gen >> 8) & 0xFF) ^ d[i * 4 + 1])
        zeros.append(((gen >> 16) & 0xFF) ^ d[i * 4 + 2])
        zeros.append(((gen >> 24) & 0xFF) ^ d[i * 4 + 3])

    # find solutions
    for sol in solver.solve(zeros):
        print("Found solution")
        state_recover = [s.get(sol) for s in state]

        # recreate random number generator
        rng_recover = bvsolver.CPythonRandom(state_recover, 624).to_cpython_random()
        words = bytearray()
        for i in range(len(d) // 4):
            gen_recover = rng_recover.getrandbits(32)
            # little endian
            words.append((gen_recover & 0xFF) ^ d[i * 4 + 0])
            words.append(((gen_recover >> 8) & 0xFF) ^ d[i * 4 + 1])
            words.append(((gen_recover >> 16) & 0xFF) ^ d[i * 4 + 2])
            words.append(((gen_recover >> 24) & 0xFF) ^ d[i * 4 + 3])
        print(words[2500:])
