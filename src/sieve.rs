use super::precalculated_values::*;
use std::sync::{Arc, Mutex};
use num::integer::sqrt;
use bit_vec::BitVec;
use std::time::SystemTime;
use threadpool::ThreadPool;
use num_cpus;

pub struct SSoZ {
    p_cnt : usize,
    num : u64,
    twins_cnt : u64,
    primes: Arc<Vec<u64>>,
    kb : usize,
    counts: Arc<Mutex<Vec<u64>>>,
    last_twins: Arc<Mutex<Vec<u64>>>,
    pos: Arc<Vec<u64>>,
    modpg:     u64,
    res_cnt:    usize,
    pairs_cnt:  usize,
    residues:  Vec<u64>,
    res_twins: Vec<u64>,
    res_inv:  Arc<Vec<u64>>,
    bn:        u64,
}

impl SSoZ {
    fn new() -> SSoZ {
        SSoZ {
            p_cnt: 0,
            num: 0,
            twins_cnt: 0,
            primes: Arc::new(Vec::new()),
            kb: 0,
            counts: Arc::new(Mutex::new(Vec::new())),
            last_twins: Arc::new(Mutex::new(Vec::new())),
            pos: Arc::new(Vec::new()),
            modpg: 0,
            res_cnt: 0,
            pairs_cnt: 0,
            residues: Vec::new(),
            res_twins: Vec::new(),
            res_inv: Arc::new(Vec::new()),
            bn: 0
        }
    }
    pub fn largest_twin_prime_before(max : u64) -> u64 {
        let mut ssoz = SSoZ::new();
        let now = SystemTime::now();
        let thread_pool = ThreadPool::new(num_cpus::get());
        ssoz.num = (max - 1) | 1;
        ssoz.select_pg();
        let modpg = ssoz.modpg;
        let mut k = ssoz.num / modpg;
        let mut mod_k = modpg * k;
        let k_max = (ssoz.num - 2) % modpg + 1;
        let b = ssoz.bn * 1024;

        ssoz.soz_pg(sqrt(ssoz.num));
        println!("setup time: {}", {
            match now.elapsed() {
                Ok(elapsed) => {
                    (elapsed.as_secs() * 1000) as f64 + elapsed.subsec_nanos() as f64 / 1_000_000f64},
                Err(e) => {panic!("Timer error {:?}", e)},
            }
        });

        ssoz.twins_cnt = if modpg > 30030u64 {4} else if modpg > 210u64 {3} else {2};
        let now = SystemTime::now();
        let mut index = 0usize;
        while index < ssoz.pairs_cnt {
            thread_pool.execute(|| {
                let ssoz = ssoz;
                SSoZ::twin_sieve(k_max, index, ssoz.kb, ssoz.res_twins[index], ssoz.modpg, ssoz.num, ssoz.p_cnt,
                        ssoz.primes.clone(), ssoz.res_inv.clone(), ssoz.pos.clone(),
                        ssoz.last_twins.clone(), ssoz.counts.clone());
            });
            index += 1;
        }
        thread_pool.join();
        let count = ssoz.counts.lock().unwrap();
        for i in 0..count.len() {
            ssoz.twins_cnt += count[i];
        }
        let last_twins = ssoz.last_twins.lock().unwrap();
        let mut last = 0u64;
        for i in 0..last_twins.len() {
            if last < last_twins[i] {last = last_twins[i];}
        }
        println!("sieve time: {}", {
            match now.elapsed() {
                Ok(elapsed) => {
                    (elapsed.as_secs() * 1000) as f64 + elapsed.subsec_nanos() as f64 / 1_000_000f64},
                Err(e) => {panic!("Timer error {:?}", e)},
            }
        });
        last
    }
    fn select_pg(&mut self) {
        if self.num < 10_000_000u64 {
            self.modpg = PARAMETERS_P5.0;
            self.res_cnt = PARAMETERS_P5.1;
            self.pairs_cnt = PARAMETERS_P5.2;
            self.residues = PARAMETERS_P5.3.to_vec();
            self.res_twins = PARAMETERS_P5.4.to_vec();
            self.res_inv = Arc::new(PARAMETERS_P5.5.to_vec());
            self.bn = 16;
        } else if self.num < 1_100_000_000u64 {
            self.modpg = PARAMETERS_P7.0;
            self.res_cnt = PARAMETERS_P7.1;
            self.pairs_cnt = PARAMETERS_P7.2;
            self.residues = PARAMETERS_P7.3.to_vec();
            self.res_twins = PARAMETERS_P7.4.to_vec();
            self.res_inv = Arc::new(PARAMETERS_P7.5.to_vec());
            self.bn = 32;
        } else if self.num < 35_500_000_000u64 {
            self.modpg = PARAMETERS_P11.0;
            self.res_cnt = PARAMETERS_P11.1;
            self.pairs_cnt = PARAMETERS_P11.2;
            self.residues = PARAMETERS_P11.3.to_vec();
            self.res_twins = PARAMETERS_P11.4.to_vec();
            self.res_inv = Arc::new(PARAMETERS_P11.5.to_vec());
            self.bn = 64;
        } else if self.num < 15_000_000_000_000u64 {
            self.modpg = PARAMETERS_P13.0;
            self.res_cnt = PARAMETERS_P13.1;
            self.pairs_cnt = PARAMETERS_P13.2;
            self.residues = PARAMETERS_P13.3.to_vec();
            self.res_twins = PARAMETERS_P13.4.to_vec();
            self.res_inv = Arc::new(PARAMETERS_P13.5.to_vec());
            if self.num > 7_000_000_000_000u64 { self.bn = 384;}
            else if self.num > 2_500_000_000_000u64 { self.bn = 320;}
            else if self.num > 250_000_000_000u64 { self.bn = 384;}
            else {self.bn = 96;}
        } else {
            self.modpg = PARAMETERS_P17.0;
            self.res_cnt = PARAMETERS_P17.1;
            self.pairs_cnt = PARAMETERS_P17.2;
            self.residues = PARAMETERS_P17.3.to_vec();
            self.res_twins = PARAMETERS_P17.4.to_vec();
            self.res_inv = Arc::new(PARAMETERS_P17.5.to_vec());
            self.bn = 384;
        }
        self.counts = Arc::new(Mutex::new(Vec::with_capacity(self.pairs_cnt)));
        self.last_twins = Arc::new(Mutex::new(Vec::with_capacity(self.pairs_cnt)));
        self.pos = Arc::new(Vec::with_capacity(self.modpg as usize));
        for i in 0..self.res_cnt-1 {
            self.pos[(self.residues[i] - 2) as usize] = i as u64;
        }
    }
    fn soz_pg(&mut self, val: u64) {
        let md = self.modpg;
        let res_cnt = self.res_cnt;

        let num = (val - 1) | 1;
        let mut k = num / md;
        let mut mod_k = md * k;
        let mut r = 0i128;
        while num >= mod_k + self.residues[r as usize] {
            r += 1;
        }
        let max_pcs = k * (res_cnt as u64) + (r as u64);
        let mut primes = BitVec::with_capacity(max_pcs as usize);
        let sqn = sqrt(num);

        mod_k = 0; r = -1; k = 0;

        for i in 0usize..max_pcs as usize {
            if {r += 1; r} == res_cnt as i128 {
                r = 0;
                mod_k += md;
                k += 1;
            }
            if primes[i] {continue;}
            let pmr_r = self.residues[i];
            let prime = mod_k + pmr_r;
            if prime > sqn {break;}
            let prm_step = prime * (res_cnt as u64);
            for ri in &self.residues {
                let prod = pmr_r * (*ri) - 2;
                let mut prm_mult = (k * (prime + *ri) + prod / md) * (res_cnt as u64) + self.pos[(prod % md) as usize];
                while prm_mult < max_pcs {
                    primes.set(prm_mult as usize, true);
                    prm_mult += prm_step;
                }
            }
        }

        mod_k = 0; r = -1;
        for prm in primes {
            if {r += 1; r} == res_cnt as i128{
                r = 0;
                mod_k += md;
            }
            if !prm {
                self.primes.push(mod_k + self.residues[r as usize]);
            }
        }
        self.p_cnt = self.primes.len();
    }
    fn next_p_init(r_hi: u64, modpg: u64, primes: Arc<Vec<u64>>, p_cnt: usize,
                   res_inv: Arc<Vec<u64>>, pos: Arc<Vec<u64>>) -> Vec<u64> {
        let mut next_p = Vec::with_capacity(p_cnt * 2);
        let r_lo = r_hi - 2;
        let (row_lo, row_hi) = (0usize, p_cnt);
        for i in 0..p_cnt {
            let prime = primes[i];
            let k = (prime - 2) / modpg;
            let r = (prime - 2) % modpg + 2;
            let r_inv = res_inv[pos[r as usize - 2] as usize];
            let mut ri = (r_lo * r_inv - 2) % modpg + 2;
            next_p[row_lo + i] = k * (prime + ri) + (r * ri - 2) / modpg;
            ri = (r_hi * r_inv - 2) % modpg + 2;
            next_p[row_hi + i] = k * (prime + ri) + (r * ri - 2) / modpg;
        }
        return next_p;
    }
    fn twin_sieve(k_max: u64, index: usize, kb: usize, r_hi: u64, modpg: u64, num: u64, p_cnt: usize,
                  primes: Arc<Vec<u64>>, res_inv: Arc<Vec<u64>>, pos: Arc<Vec<u64>>,
                  last_twins: Arc<Mutex<Vec<u64>>>, count: Arc<Mutex<Vec<u64>>>) {
        let (mut sum, mut ki, mut kn) = (0u64, 0u64, 0u64);
        let (mut hi_tp, mut upk) = (0u64, 0usize);
        let mut k_hi = 0u64;
        let mut last_tw = 0u64;
        let mut seg: Vec<u8> = Vec::with_capacity(kb);
        let mut next_p = SSoZ::next_p_init(r_hi, modpg,
               primes.clone(), p_cnt, res_inv.clone(), pos.clone());
        while ki < k_max {
            kn = {if (kb as u64) < (k_max - ki) {kb as u64} else {k_max - ki}};
            for b in 0..kn as usize {seg[b] = 0;}
            for i in 0..p_cnt {
                let mut k = next_p[i];
                while k < kn {
                    seg[k as usize] |= 1;
                    k += primes[i];
                }
                next_p[i] = k - kn;

                k = next_p[p_cnt + i];
                while k < kn {
                    seg[k as usize] |= 1;
                    k += primes[i];
                }
                next_p[i + p_cnt] = k - kn;
            }
            let mut cnt = 0u64;
            for k in 0..kn as usize {
                if seg[k] == 0 {cnt += 1;}
            }
            if cnt > 0 {
                sum += cnt;
                for k in 1..=kn as usize {
                    if seg[(kn as usize) - k] == 0 {
                        upk = (kn as usize) - k;
                        break;
                    }
                }
                k_hi = hi_tp;
                hi_tp = ki + upk as u64;
            }
            ki += kb as u64;
        }
        let mut mod_k = hi_tp * modpg;
        if mod_k + r_hi > num {
            for k in 0..upk + 1 {
                if seg[upk - k] == 0 {
                    hi_tp = mod_k + r_hi;
                    if hi_tp <= num {last_tw = hi_tp; break;}
                    sum -= 1;
                }
                mod_k -= modpg;
            }
            last_tw = if r_hi > num {0u64} else {k_hi * modpg + r_hi};
        } else {last_tw = mod_k + r_hi}
        last_twins.lock().unwrap()[index] = last_tw;
        count.lock().unwrap()[index] = sum;
    }
}