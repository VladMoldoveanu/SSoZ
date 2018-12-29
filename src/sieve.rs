use super::precalculated_values::*;
use std::sync::{Arc, Mutex};
use num::integer::sqrt;
use std::time::SystemTime;
use threadpool::ThreadPool;
use num_cpus;
use memsec::memset;

pub fn largest_twin_prime_before(max : usize) -> (usize, usize) {
    let now = SystemTime::now();

    //Initialize with precompiled data
    let num = (max - 1) | 1;
    let (modpg, res_cnt, pairs_cnt, bn
        , residues, res_twins, res_inv) = select_pg(num);
    let counts = Arc::new(Mutex::new(vec![0usize;pairs_cnt]));
    let last_twins = Arc::new(Mutex::new(vec![0usize;pairs_cnt]));
    let mut pos = vec![0usize; modpg];
    for i in 0..res_cnt {
        pos[residues[i] - 2] = i;
    }
    let pos = Arc::new(pos);


    let thread_pool = ThreadPool::new(num_cpus::get());
    let _k = num / modpg;
    let k_max = (num - 2) / modpg + 1;
    let b = bn << 10;
    let kb = if k_max < b {k_max} else {b};
    let mut twins_cnt = if modpg > 30030usize {4usize} else if modpg > 210usize {3usize} else {2usize};

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
    let mut last = 0usize;
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

fn twin_sieve(k_max: usize, index: usize, kb: usize, r_hi: usize, modpg: usize, num: usize, p_cnt: usize,
              primes: Arc<Vec<usize>>, res_inv: Arc<Vec<usize>>, pos: Arc<Vec<usize>>,
              last_twins: Arc<Mutex<Vec<usize>>>, count: Arc<Mutex<Vec<usize>>>) {
    let (mut sum, mut ki) = (0usize, 0usize);
    let (mut hi_tp, mut upk) = (0usize, 0usize);
    let mut k_hi = 0usize;
    let mut seg: Vec<u8> = vec![0u8;kb];
    let mut next_p =
        next_p_init(r_hi, modpg, primes.clone(),p_cnt, res_inv.clone(), pos.clone());

    while ki < k_max {
        let kn = {if kb < (k_max - ki) {kb} else {k_max - ki}};
        unsafe {
            memset(seg.as_mut_ptr(), 0, kn);
        }
        //for b in 0..kn {seg[b] = 0;}
        for (i, prime) in primes.iter().enumerate() {
            let mut k = next_p[i];
            while k < kn {
                seg[k] |= 1;
                k += prime;
            }
            next_p[i] = k - kn;

            k = next_p[p_cnt + i];
            while k < kn {
                seg[k] |= 1;
                k += prime;
            }
            next_p[i + p_cnt] = k - kn;
        }
        let mut cnt = 0usize;
        seg.iter().take(kn).for_each(|&x| if x == 0 {cnt += 1});

        if cnt > 0 {
            sum += cnt;
            for k in 1..=kn {
                if seg[kn - k] == 0 {
                    upk = kn - k;
                    break;
                }
            }
            k_hi = hi_tp;
            hi_tp = ki + upk;
        }
        ki += kb;
    }
    hi_tp = hi_tp * modpg + r_hi;
    if hi_tp > num {
        let mut prev = true;
        for k in 0..upk + 1 {
            if seg[upk - k] == 0 {
                if hi_tp <= num {prev = false; break;}
                sum -= 1;
            }
            hi_tp -= modpg;
        }
        if prev {
            hi_tp = if r_hi > num {0usize} else {k_hi * modpg + r_hi};
        }
    }
    last_twins.lock().unwrap()[index] = hi_tp;
    count.lock().unwrap()[index] = sum;
}

fn next_p_init(r_hi: usize, modpg: usize, primes: Arc<Vec<usize>>, p_cnt: usize,
               res_inv: Arc<Vec<usize>>, pos: Arc<Vec<usize>>) -> Vec<usize> {
    let mut next_p = vec![0usize; p_cnt * 2];
    let r_lo = r_hi - 2;
    let (row_lo, row_hi) = (0usize, p_cnt);
    for (i, prime) in primes.iter().enumerate() {
        let k = (prime - 2) / modpg;
        let r = (prime - 2) % modpg + 2;
        let r_inv = res_inv[pos[r - 2]];
        let mut ri = (r_lo * r_inv - 2) % modpg + 2;
        next_p[row_lo + i] = k * (prime + ri) + (r * ri - 2) / modpg;
        ri = (r_hi * r_inv - 2) % modpg + 2;
        next_p[row_hi + i] = k * (prime + ri) + (r * ri - 2) / modpg;
    }
    return next_p;
}

fn select_pg(num: usize) -> (usize, usize, usize, usize, Vec<usize>, Vec<usize>, Arc<Vec<usize>>) {
    if num < 10_000_000usize {
         return (PARAMETERS_P5.0 , PARAMETERS_P5.1, PARAMETERS_P5.2, 16, PARAMETERS_P5.3.to_vec(),
                 PARAMETERS_P5.4.to_vec(), Arc::new(PARAMETERS_P5.5.to_vec()));
    }
    if num < 1_100_000_000usize {
        return (PARAMETERS_P7.0 , PARAMETERS_P7.1, PARAMETERS_P7.2, 32, PARAMETERS_P7.3.to_vec(),
                PARAMETERS_P7.4.to_vec(), Arc::new(PARAMETERS_P7.5.to_vec()));
    }
    if num < 35_500_000_000usize {
        return (PARAMETERS_P11.0 , PARAMETERS_P11.1, PARAMETERS_P11.2, 64, PARAMETERS_P11.3.to_vec(),
                PARAMETERS_P11.4.to_vec(), Arc::new(PARAMETERS_P11.5.to_vec()));
    }
    if num < 15_000_000_000_000usize {
        let mut bn = {if num > 7_000_000_000_000usize {384}
        else if num > 2_500_000_000_000usize {320}
        else if num > 250_000_000_000usize {196}
        else {96} };
        return (PARAMETERS_P13.0 , PARAMETERS_P13.1, PARAMETERS_P13.2, bn, PARAMETERS_P13.3.to_vec(),
                PARAMETERS_P13.4.to_vec(), Arc::new(PARAMETERS_P13.5.to_vec()));
    }
    return (PARAMETERS_P17.0 , PARAMETERS_P17.1, PARAMETERS_P17.2, 384, PARAMETERS_P17.3.to_vec(),
            PARAMETERS_P17.4.to_vec(), Arc::new(PARAMETERS_P17.5.to_vec()));
}

fn soz_pg(val: usize, md: usize, res_cnt: usize, residues: &Vec<usize>, pos: Arc<Vec<usize>>)
    -> (Arc<Vec<usize>>, usize) {

    let num = (val - 1) | 1;
    let mut k = num / md;
    let mut mod_k = md * k;
    let mut r = 0usize;

    let mut prms = Vec::new();
    while num >= mod_k + residues[r] {
        r += 1;
    }
    let max_pcs = k * res_cnt + r;
    let mut primes = vec![false; max_pcs];
    let sqn = sqrt(num);

    mod_k = 0; r = 0; k = 0;

    for i in 0usize..max_pcs {
        if r == res_cnt {
            r = 0;
            mod_k += md;
            k += 1;
        }
        if primes[i] {r += 1; continue;}
        let pmr_r = residues[i];
        let prime = mod_k + pmr_r;
        if prime > sqn {break;}
        let prm_step = prime * (res_cnt);
        for ri in residues {
            let prod = pmr_r * ri - 2;
            let mut prm_mult = (k * (prime + ri) + prod / md) * res_cnt + pos[prod % md];
            while prm_mult < max_pcs {
                primes[prm_mult] = true;
                prm_mult += prm_step;
            }
        }
        r += 1;
    }

    mod_k = 0; r = 0;
    for prm in primes {
        if r == res_cnt{
            r = 0;
            mod_k += md;
        }
        if !prm {
            prms.push(mod_k + residues[r]);
        }
        r += 1;
    }
    let len = prms.len();
    (Arc::new(prms), len)
}