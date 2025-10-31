# bvsolver

Solve equations on bitvectors over GF(2).

Supported operations (BitVector means a bitvector of unknown bits to solve):

1. BitVector ^ BitVector
2. BitVector ^ int
3. BitVector & int
4. 1-bit BitVector * int
5. BitVector >> int
6. BitVector << int

Heavily inspired by [gf2bv](https://github.com/maple3142/gf2bv).

## Usage

Install the Python package:

```shell
pip3 install git+ssh://git@github.com/jiegec/bvsolver.git
```

See [examples](./python/examples) for usage.
