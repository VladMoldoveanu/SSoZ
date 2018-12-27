extern crate SSoZ;
use SSoZ::sieve::largest_twin_prime_before;

fn main() {
    let (last, count) = largest_twin_prime_before(1000000000);
    println!("largest twin prime: {} \ntwin primes found: {}", last, count);
}