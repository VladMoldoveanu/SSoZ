extern crate SSoZ;
use SSoZ::sieve::largest_twin_prime_before;

fn main() {
    let (last, count) = largest_twin_prime_before(10000000000);
    println!("largest twin prime: {} \ntwin primes found: {}", last, count);
}