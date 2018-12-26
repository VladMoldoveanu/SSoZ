use num::integer::gcd;
use std::mem::swap;

fn mod_inverse(a0 : u64, b0: u64) -> u64 {
    if b0 == 1 {
        return 1;
    }
    let mut a = a0 as i128;
    let mut b = b0 as i128;
    let mut x0 = 0i128;
    let mut result = 1;

    while a > 1 {
        let q = a / b;
        a = a % b;
        swap(&mut a, &mut b);
        result -= q * x0;
        swap(&mut x0, &mut result);
    }
    if result < 0 {
        result += b0 as i128;
    }
    return result as u64;
}
pub fn gen_pg_parameters(prime: u64) -> String {
    println!("Generating pg parameters for prime: {}", prime);
    let primes = vec![2u64, 3, 5, 7, 11, 13, 17, 19, 23];
    let mut modpg = 1u64;
    //PG modulus
    for prm in primes {
        modpg *= prm;
        if prm == prime {break;}
    }

    let mut residues : Vec<u64> = vec![];
    let mut pc = 5u64;
    let mut inc = 2u64;

    while pc < (modpg >> 1) {
        if gcd(modpg, pc) == 1 {
            residues.push(pc);
            residues.push(modpg - pc);
        }
        pc += inc;
        inc ^= 0b110;
    }
    residues.sort();
    residues.push(modpg - 1);
    residues.push(modpg + 1);
    let res_cnt = residues.len();

    let mut res_twins : Vec<u64> = vec![];
    let mut j = 0usize;
    while j < res_cnt - 1 {
        if residues[j] + 2 == residues[j + 1] {
            res_twins.push(residues[j + 1]);
        }
        j += 1;
    }
    let twin_pairs = res_twins.len();

    let mut inverses : Vec<u64> = vec![];
    for res in &residues {
        inverses.push(mod_inverse(*res, modpg));
    }
    //The usual return value
    //(modpg, res_cnt as u64, twin_pairs as u64, residues, res_twins, inverses)
    let mut formatted = String::new();
    formatted.push_str(&format!("pub static PARAMETERS_P{}: (u64, u64, u64, [u64;{}], [u64;{}], [u64;{}]) = ",
                                        prime, res_cnt, twin_pairs, res_cnt));
    formatted.push_str(&format!("({}", modpg));
    formatted.push_str(&format!(", {}", res_cnt));
    formatted.push_str(&format!(", {}", twin_pairs));
    formatted.push_str(&format!(", {:?}", residues));
    formatted.push_str(&format!(", {:?}", res_twins));
    formatted.push_str(&format!(", {:?});", inverses));
    formatted
}