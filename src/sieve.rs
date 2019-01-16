use super::precalculated_values::*;
use std::sync::{Arc, mpsc};
use std::thread;
use num::integer::sqrt;
use memsec::memzero;

// Returns the largest twin prime (+-1) as well as the number of twin primes before 'max'
pub fn largest_twin_prime_before(max : usize) -> (usize, usize) {

    //Initialize with precompiled data
    let num = (max - 1) | 1;
    let (modpg, res_cnt, pairs_cnt, bn
        , residues, res_twins, res_inv) = select_pg(num);
    let mut pos = vec![0usize; modpg];
    for i in 0..res_cnt {pos[residues[i] - 2] = i;}
    let pos = Arc::new(pos);

    let k_max = (num - 2) / modpg + 1;
    let b = bn << 10;
    let kb = if k_max < b {k_max} else {b};
    let mut twins_cnt = if modpg > 30030usize {4usize} else if modpg > 210usize {3usize} else {2usize};

    let (primes, p_cnt) = soz_pg(sqrt(num), &residues);

    //Create a channel to send the results from the threads
    let (sender, receiver) = mpsc::channel();

    for index in 0..pairs_cnt {
        //Clone the necessary data
        let (k_max_0, kb_0, r_hi, modpg_0, num_0, p_cnt_0,
            primes_0, res_inv_0, pos_0, sender_0) =
            (k_max, kb, res_twins[index], modpg, num, p_cnt,
             primes.clone(), res_inv.clone(), pos.clone(), sender.clone());
        //Execute 'twin_sieve' in parallel, sending the results to the channel
        thread::spawn(move || {
            sender_0.send(twin_sieve(k_max_0, kb_0, r_hi, modpg_0,
                             num_0, p_cnt_0, primes_0, res_inv_0, pos_0))
                .expect("Unable to send to channel");
        });
    }

    //Collect the results from the channel, as they finish
    let mut last = 0usize;
    for _ in 0..pairs_cnt {
        let (l, c) = receiver.recv().expect("Unable to receive from channel");
        twins_cnt += c;
        if l > last {last = l;}
    }
    (last - 1, twins_cnt)
}

// Perform in a thread, the ssoz for a given twin pair, for k_max resgroups.
// First create|init 'next_p' array of 1st prime mults for given twin pair and
// its seg array of KB bytes, which will be gc'd|recovered at end of thread.
// For sieve, mark seg byte to '1' if either twin pair restrack is nonprime,
// for primes mults resgroups, update 'next_p' restrack slices accordingly.
// Then find last twinprime|sum <= num, return it with the total number of twin primes found.
fn twin_sieve(k_max: usize, kb: usize, r_hi: usize, modpg: usize, num: usize, p_cnt: usize,
              primes: Arc<Vec<usize>>, res_inv: Arc<Vec<usize>>, pos: Arc<Vec<usize>>) -> (usize, usize) {
    //Initialize variables
    let (mut sum, mut ki, mut kn) = (0usize, 0usize, kb);
    let (mut hi_tp, mut upk, mut k_max) = (0usize, 0usize, k_max);
    let mut seg = vec![0u8;kb];
    let mut next_p =
        next_p_init(r_hi, modpg, primes.clone(),p_cnt, res_inv.clone(), pos.clone());

    //Consider all resgroup size slices up to k_max
    if (k_max - 1) * modpg + r_hi > num {k_max -= 1}
    while ki < k_max {
        if kb > (k_max - ki) {kn = k_max - ki;}
        unsafe {memzero(seg.as_mut_ptr(), kn);}
        //For each prime, mark the multiples of the twin pair
        for (i, &prime) in primes.iter().enumerate() {
            //lower twin
            let mut k = next_p[i];
            while k < kn {
                seg[k] |= 1;
                k += prime;
            }
            next_p[i] = k - kn;

            //higher twin
            k = next_p[p_cnt + i];
            while k < kn {
                seg[k] |= 1;
                k += prime;
            }
            next_p[i + p_cnt] = k - kn;
        }

        //count the number of twins found
        let mut cnt = 0usize;
        seg.iter().take(kn).for_each(|&x| if x == 0 {cnt += 1});

        //Add the number of twin primes found to the local counter 'sum'
        if cnt > 0 {
            sum += cnt;
            // Save the location of the largest prime
            upk = kn - 1;
            while seg[upk] == 1 {upk -= 1}
            hi_tp = ki + upk;
        }
        ki += kb;
    }
    //Numerate largest twin prime for thread
    hi_tp = if r_hi > num {0usize} else {hi_tp * modpg + r_hi};
    (hi_tp, sum)
}

// Initialize 'next_p' array for given twin pair in res_twins.
// Set each row[j] w/1st prime multiple resgroup for each prime r1..sqrt(N).
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

//Selects the desired precompiled data which matches the sieve limit
fn select_pg(num: usize) -> (usize, usize, usize, usize, Vec<usize>, Vec<usize>, Arc<Vec<usize>>) {
    if num < 100_000usize {
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
        let bn = {if num > 7_000_000_000_000usize {384}
        else if num > 2_500_000_000_000usize {320}
        else if num > 250_000_000_000usize {196}
        else {96} };
        return (PARAMETERS_P13.0 , PARAMETERS_P13.1, PARAMETERS_P13.2, bn, PARAMETERS_P13.3.to_vec(),
                PARAMETERS_P13.4.to_vec(), Arc::new(PARAMETERS_P13.5.to_vec()));
    }
    return (PARAMETERS_P17.0 , PARAMETERS_P17.1, PARAMETERS_P17.2, 384, PARAMETERS_P17.3.to_vec(),
            PARAMETERS_P17.4.to_vec(), Arc::new(PARAMETERS_P17.5.to_vec()));
}

//Computes the primes in r1..sqrt(val) - any algorithm can be used (fast|small)
//Here the SoZ for P5 is used, 'residues' is for selected PG (only need 'residues[0]')
fn soz_pg(val: usize, residues: &Vec<usize>) -> (Arc<Vec<usize>>, usize) {
    //Create parameters for P5
    let md = 30;
    let res_cnt = 8;
    let res: [usize, 8] = [7,11,13,17,19,23,29,31];
    let posn: [usize; 30] = [0,0,0,0,0,0,0,0,0,1,0,2,0,0,0,3,0,4,0,0,0,5,0,0,0,0,0,6,0,7];

    let num = (val - 1) | 1;
    let mut k = num / md;
    let mut mod_k = md * k;
    let mut r = 0usize;

    let mut prms = Vec::new();
    while num >= mod_k + res[r] {r += 1;}
    let max_pcs = k * res_cnt + r;
    let mut primes = vec![false; max_pcs];
    let sqn = sqrt(num);

    mod_k = 0; r = 0; k = 0;

    //For each prime, mark its multiples
    for i in 0usize..max_pcs {
        if r == res_cnt {r = 0; mod_k += md; k += 1;}
        if primes[i] {r += 1; continue;}
        let pmr_r = res[i];
        let prime = mod_k + pmr_r;
        if prime > sqn {break;}
        let prm_step = prime * res_cnt;
        for ri in res {
            //compute resgroup val of 1st prime multiple, then mark all prime multiples up to end of prms
            let prod = pmr_r * ri - 2;
            let mut prm_mult = (k * (prime + ri) + prod / md) * res_cnt + posn[prod % md];
            while prm_mult < max_pcs {primes[prm_mult] = true; prm_mult += prm_step;}
        }
        r += 1;
    }
    //Extract the primes from prms
    mod_k = 0; r = 0;
    for prm in primes {
        if r == res_cnt {r = 0; mod_k += md;}
        if !prm {prms.push(mod_k + res[r]);}
        r += 1;
    }
    //Discard in prms any of modpg's base primes for selected PG
    while prms[0] < residues[0] {prms.remove(0);}
    let len = prms.len();
    (Arc::new(prms), len)
}
