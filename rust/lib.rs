#[pyo3::pymodule]
mod _lib {
    use num_bigint::{BigInt, BigUint, ToBigInt, ToBigUint};
    use pyo3::prelude::*;
    use pyo3::types::PyInt;

    #[pyfunction]
    fn solve(size: usize, zeros: Vec<Bound<'_, PyInt>>) -> PyResult<(BigInt, Vec<BigInt>)> {
        let zeros: Vec<BigUint> = zeros
            .iter()
            .map(|num| num.extract::<BigInt>().unwrap().to_biguint().unwrap())
            .collect();
        // gauss elimination
        let mut done = 0;
        // for each bit, record its row with the bit set to 1
        let mut mapping = vec![BigUint::ZERO; size];
        let one = 1.to_biguint().unwrap();
        for e in zeros {
            if e == one {
                // no solution
                return Ok((-1.to_bigint().unwrap(), vec![]));
            }

            // eliminate with existing rows in reduced row echelon form
            let mut e = e;
            let mut added = false;
            let mut e_no_low = &e ^ (&e & &one);
            while e_no_low != BigUint::ZERO {
                // minus 1: 0b10 -> mapping index 0
                let bit = e_no_low.trailing_zeros().unwrap() as usize - 1;
                if mapping[bit] != BigUint::ZERO {
                    e = &e ^ &mapping[bit];
                    e_no_low = &e ^ (&e & &one);
                } else {
                    mapping[bit] = e.clone();
                    done += 1;
                    added = true;
                    break;
                }
            }

            if !added && e != BigUint::ZERO {
                // eliminate to non-zero, no solution
                return Ok((-1.to_bigint().unwrap(), vec![]));
            }
        }

        // find free variables
        let mut free = vec![];
        for (i, row) in mapping.iter().enumerate() {
            if row == &BigUint::ZERO {
                free.push(i);
            }
        }

        // convert to row reduced echelon form and compute solution if all free variables are zero
        let mut base_res = one.clone(); // constant term
        for i in 0..size {
            for j in (i + 1)..size {
                if mapping[i].bit((j + 1) as u64) {
                    // eliminate
                    mapping[i] = &mapping[i] ^ &mapping[j];
                }
            }
        }
        for (i, row) in mapping.iter().enumerate() {
            // after elimination, compute solution if free variables are zero
            base_res.set_bit(i as u64 + 1, row.bit(0));
        }

        // now all free variables are found
        assert_eq!(done + free.len(), size);

        // compute contribution of each free variable to solution
        let mut contributions = vec![];
        for bit in free {
            let mut c = BigInt::ZERO;
            for (i, row) in mapping.iter().enumerate() {
                if i == bit || row.bit((bit + 1) as u64) {
                    c.set_bit((i + 1) as u64, true);
                }
            }
            contributions.push(c);
        }
        Ok((base_res.to_bigint().unwrap(), contributions))
    }
}
