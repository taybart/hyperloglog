/*
 * Add:
 *  hash input data v with function h
 *  x = h(v)
 *  get the first b bits where b is log2(m)
 *  and add 1 to abtain the register to modify
 *  j = 1 + <x1 x2...xb>2
 *  w is bits of x
 *  w = xb+1 xb+2 ...
 *  roe is the index of the leftmost 1
 *  then take the max of M[j] and roe of w
 *  M[j] = max(M[j], p(w))
 *
 * Count:
 *   sum and get harmonic mean and correct with alpha
 *   Z = 1 / sum(2^-M[j])
 *   estimate = alpha * m^2 * Z
 */

use std::{
    collections::hash_map::DefaultHasher,
    hash::{Hash, Hasher},
};

const TWO_POW_32: f64 = (1_i64 << 32_i64) as f64;
const ERROR_RATE: f64 = 0.008;

pub struct HyperLogLog {
    pub b: usize,
    hasher: DefaultHasher,
    pub m_counters: usize,
    m: Vec<u8>,
    alpha: f64,
}

impl HyperLogLog {
    pub fn new() -> Result<Self, String> {
        /*
         * b, amount of bits to use from hash
         *   b = log2(m)
         * m, number of counter buckets
         *   error_rate = 1.04 / sqrt(m)
         *   m = (1.04 / error_rate)^2
         * therefore,
         *   b = log2((1.04 / error_rate)^2)
         */
        let b = (1.04 / ERROR_RATE).powi(2).log2().ceil() as usize;

        // m, number of counter buckets in M
        let m_counters = 1 << b;

        // using alpha approximation
        // https://en.wikipedia.org/wiki/HyperLogLog#Practical_considerations
        if m_counters < 16 {
            return Err("error rate too high counters are below 16".into());
        }
        let alpha = match m_counters {
            16 => 0.673,
            32 => 0.697,
            64 => 0.709,
            _ => 0.7213 / (1.0 + 1.079 / m_counters as f64),
        };

        Ok(HyperLogLog {
            b,
            hasher: DefaultHasher::new(),
            m: vec![0; m_counters],
            m_counters,
            alpha,
        })
    }

    pub fn add(&mut self, data: &str) {
        let mut hasher = self.hasher.clone();
        data.hash(&mut hasher);
        let hash = hasher.finish();

        // get counter index using the precision bits aka sigma
        let j = (hash >> (64 - self.b)) as usize;
        let estimate_bits = hash & ((1 << (64 - self.b)) - 1);

        // Count the number of leading zeros
        let trailing_zeros = estimate_bits.leading_zeros() as u8 + 1;

        self.m[j] = std::cmp::max(self.m[j], trailing_zeros);
    }

    fn range_correction(&self, estimate: f64, m_counters: f64, empty_counters: usize) -> f64 {
        // magic adjustments from section 4 of the paper
        // https://algo.inria.fr/flajolet/Publications/FlFuGaMe07.pdf
        match estimate {
            // small range correction
            small_range if estimate <= 2.5 * m_counters => {
                if empty_counters > 0 {
                    m_counters * (m_counters / empty_counters as f64).ln()
                } else {
                    small_range
                }
            }
            // intermediate range correction - no correction
            intermediate_range if estimate <= (TWO_POW_32 / 30.0) => intermediate_range,
            // large range correction
            large_range => -TWO_POW_32 * (1.0 - large_range / TWO_POW_32).ln(),
        }
    }

    pub fn count(&self) -> usize {
        let mut z = 0.0;
        let mut empty_counters = 0;
        for &v in self.m.iter() {
            z += 2.0_f64.powf(-(v as f64));
            if v == 0 {
                empty_counters += 1;
            }
        }

        /* sum and get harmonic mean and correct with alpha
         * Z = 1 / sum(2^-M[j])
         * estimate = alpha * m^2 * Z */
        let m_counters = self.m_counters as f64;
        let estimate = self.alpha * m_counters.powi(2) / z;
        self.range_correction(estimate, m_counters, empty_counters) as usize
    }

    pub fn merge(mut self, other: HyperLogLog) {
        for j in 0..self.m.len() {
            self.m[j] = std::cmp::max(self.m[j], other.m[j]);
        }
    }

    pub fn error(self) -> f64 {
        1.04 / ((1 << self.b) as f64).sqrt()
    }
}
