extern crate SSoZ;
use SSoZ::sieve::largest_twin_prime_before;
use std::time::SystemTime;
use std::env;
use std::process;
use std::io::Write;

fn time(f: fn(usize) -> (usize, usize), x: usize) -> u32 {
    let now = SystemTime::now();
    let _res = f(x);

    match now.elapsed() {
        Ok(e) => {e.subsec_micros() + e.as_secs() as u32 * 1_000_000},
        Err(e) => {panic!("Timer error {:?}", e)},
    }
}

fn main() {
    let usage = "Usage: [sieve limit] [no. tests]";
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        println!("Wrong number of arguments\n{}",usage);
        process::exit(1);
    }
    let n = match args[1].parse::<usize>() {
        Ok(x) => {x},
        Err(e) => {
            println!("First argument is not usize: {:?}\n{}", e, usage);
            process::exit(-1);
        },
    };
    let x = match args[2].parse::<u64>() {
        Ok(x) => {x},
        Err(e) => {
            println!("Second argument is not usize: {:?}\n{}", e, usage);
            process::exit(-1);
        },
    };
    let mut total_time = 0u64;
    let mut min_time = 4294967295u32;
    let mut max_time = 0u32;
    let mut step = 0;
    for i in 0..x {
        if (x <= 10) || ((i * 10 / x) > step) {
            step += 1;
            print!("#{}", step);
            std::io::stdout().flush().unwrap();
        }
        let t = time(largest_twin_prime_before, n);
        if min_time > t {min_time = t;}
        if max_time < t {max_time = t;}
        total_time += t as u64;
    }
    println!();
    println!("Executed {} tests:", x);
    println!("Average time: {}s", (total_time/x) as f64 /1_000_000f64);
    println!("Min time: {}s", min_time as f64 /1_000_000f64);
    println!("Max time: {}s", max_time as f64 /1_000_000f64);
}