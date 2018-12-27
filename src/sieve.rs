use super::precalculated_values::*;
use std::sync::{Arc, Mutex};
use num::integer::sqrt;
use bit_vec::BitVec;
use std::time::SystemTime;
use threadpool::ThreadPool;
use num_cpus;

pub fn largest_twin_prime_before(max : u64) -> (u64, u64) {
    let now = SystemTime::now();

    //Initialize with precompiled data
    let num = (max - 1) | 1;
    let (modpg, res_cnt, pairs_cnt, bn
        , residues, res_twins, res_inv) = select_pg(num);
    let counts = Arc::new(Mutex::new(vec![0u64;pairs_cnt]));
    let last_twins = Arc::new(Mutex::new(vec![0u64;pairs_cnt]));
    let mut pos = vec![0u64;modpg as usize];
    for i in 0..res_cnt-1 {
        pos[(residues[i] - 2) as usize] = i as u64;
    }
    let pos = Arc::new(pos);


    let thread_pool = ThreadPool::new(num_cpus::get());
    let _k = num / modpg;
    let k_max = (num - 2) / modpg + 1;
    let b = bn << 10;
    let kb = (if k_max < b {k_max} else {b}) as usize;
    let mut twins_cnt = if modpg > 30030u64 {4u64} else if modpg > 210u64 {3u64} else {2u64};

    let (primes, p_cnt) =
        soz_pg(sqrt(num), modpg, res_cnt, &residues, pos.clone());

    println!("setup time: {}ms", {
        match now.elapsed() {
            Ok(elapsed) => {
                (elapsed.as_secs() * 1000) as f64 + elapsed.subsec_nanos() as f64 / 1_000_000f64},
            Err(e) => {panic!("Timer error {:?}", e)},
        }
    });


    let now = SystemTime::now();

    let mut index = 0usize;
    while index < pairs_cnt {
        let (k_max_0, index_0, kb_0, r_hi, modpg_0, num_0, p_cnt_0,
            primes_0, res_inv_0, pos_0,
            last_twins_0, counts_0) =
            (k_max, index, kb, res_twins[index], modpg, num, p_cnt,
             primes.clone(), res_inv.clone(), pos.clone(),
             last_twins.clone(), counts.clone());
        thread_pool.execute(move || {
            twin_sieve(k_max_0, index_0, kb_0, r_hi, modpg_0, num_0, p_cnt_0,
                             primes_0, res_inv_0, pos_0, last_twins_0, counts_0);
        });
        index += 1;
    }
    thread_pool.join();
    let count = counts.lock().unwrap();
    for i in 0..count.len() {
        twins_cnt += count[i];
    }
    let last_twins = last_twins.lock().unwrap();
    let mut last = 0u64;
    for i in 0..last_twins.len() {
        if last < last_twins[i] {last = last_twins[i];}
    }
    println!("sieve time: {}ms", {
        match now.elapsed() {
            Ok(elapsed) => {
                (elapsed.as_secs() * 1000) as f64 + elapsed.subsec_nanos() as f64 / 1_000_000f64},
            Err(e) => {panic!("Timer error {:?}", e)},
        }
    });
    (last, twins_cnt)
}

fn twin_sieve(k_max: u64, index: usize, kb: usize, r_hi: u64, modpg: u64, num: u64, p_cnt: usize,
              primes: Arc<Vec<u64>>, res_inv: Arc<Vec<u64>>, pos: Arc<Vec<u64>>,
              last_twins: Arc<Mutex<Vec<u64>>>, count: Arc<Mutex<Vec<u64>>>) {
    let (mut sum, mut ki, mut kn) = (0u64, 0u64, 0u64);
    let (mut hi_tp, mut upk) = (0u64, 0usize);
    let mut k_hi = 0u64;
    let mut last_tw = 0u64;
    let mut seg: Vec<u8> = vec![0u8;kb];
    let mut next_p = next_p_init(r_hi, modpg,
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

fn next_p_init(r_hi: u64, modpg: u64, primes: Arc<Vec<u64>>, p_cnt: usize,
               res_inv: Arc<Vec<u64>>, pos: Arc<Vec<u64>>) -> Vec<u64> {
    let mut next_p = vec![0u64; p_cnt * 2];
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

fn select_pg(num: u64) -> (u64, usize, usize, u64, Vec<u64>, Vec<u64>, Arc<Vec<u64>>) {
    if num < 10_000_000u64 {
         return (PARAMETERS_P5.0 , PARAMETERS_P5.1, PARAMETERS_P5.2, 16, PARAMETERS_P5.3.to_vec(),
                 PARAMETERS_P5.4.to_vec(), Arc::new(PARAMETERS_P5.5.to_vec()));
    }
    if num < 1_100_000_000u64 {
        return (PARAMETERS_P7.0 , PARAMETERS_P7.1, PARAMETERS_P7.2, 32, PARAMETERS_P7.3.to_vec(),
                PARAMETERS_P7.4.to_vec(), Arc::new(PARAMETERS_P7.5.to_vec()));
    }
    if num < 35_500_000_000u64 {
        return (PARAMETERS_P11.0 , PARAMETERS_P11.1, PARAMETERS_P11.2, 64, PARAMETERS_P11.3.to_vec(),
                PARAMETERS_P11.4.to_vec(), Arc::new(PARAMETERS_P11.5.to_vec()));
    }
    if num < 15_000_000_000_000u64 {
        let mut bn = 0u64;
        if num > 7_000_000_000_000u64 { bn = 384;}
        else if num > 2_500_000_000_000u64 { bn = 320;}
        else if num > 250_000_000_000u64 { bn = 384;}
        else {bn = 96;}
        return (PARAMETERS_P13.0 , PARAMETERS_P13.1, PARAMETERS_P13.2, bn, PARAMETERS_P13.3.to_vec(),
                PARAMETERS_P13.4.to_vec(), Arc::new(PARAMETERS_P13.5.to_vec()));
    }
    return (PARAMETERS_P17.0 , PARAMETERS_P17.1, PARAMETERS_P17.2, 384, PARAMETERS_P17.3.to_vec(),
            PARAMETERS_P17.4.to_vec(), Arc::new(PARAMETERS_P17.5.to_vec()));
}

fn soz_pg(val: u64, md: u64, res_cnt: usize, residues: &Vec<u64>, pos: Arc<Vec<u64>>)
    -> (Arc<Vec<u64>>, usize) {

    let num = (val - 1) | 1;
    let mut k = num / md;
    let mut mod_k = md * k;
    let mut r = 0i128;

    let mut prms = Vec::new();
    while num >= mod_k + residues[r as usize] {
        r += 1;
    }
    let max_pcs = k * (res_cnt as u64) + (r as u64);
    let mut primes = BitVec::from_elem(max_pcs as usize, false);
    let sqn = sqrt(num);

    mod_k = 0; r = -1; k = 0;

    for i in 0usize..max_pcs as usize {
        if {r += 1; r} == res_cnt as i128 {
            r = 0;
            mod_k += md;
            k += 1;
        }
        if primes[i] {continue;}
        let pmr_r = residues[i];
        let prime = mod_k + pmr_r;
        if prime > sqn {break;}
        let prm_step = prime * (res_cnt as u64);
        for ri in residues {
            let prod = pmr_r * (*ri) - 2;
            let mut prm_mult = (k * (prime + *ri) + prod / md) * (res_cnt as u64) + pos[(prod % md) as usize];
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
            prms.push(mod_k + residues[r as usize]);
        }
    }
    let len = prms.len();
    (Arc::new(prms), len)
}