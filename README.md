# bvsolver

Solve equations on bitvectors over GF(2).

Supported operations (BitVector means a bitvector of unknown bits to solve):

1. BitVector ^ BitVector
2. BitVector ^ int
3. BitVector & int
4. 1-bit BitVector * int
5. BitVector >> int
6. BitVector << int

Heavily inspired by [gf2bv](https://github.com/maple3142/gf2bv). Instead of using m4ri to solve equations, a Rust solver is implemented to accelerate solving. However, the performance is still worse than m4ri, so gf2bv is preferred.

## Usage

Install the Python package:

```shell
pip3 install git+ssh://git@github.com/jiegec/bvsolver.git
```

See [examples](./python/examples) for usage.
