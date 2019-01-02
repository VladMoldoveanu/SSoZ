# SSoZ
This is an implementation of the Twin Prime Segmented Sieve of Zakiya (SSoZ). 
You can find the thread with the discussions, as well the paper about the algorithm and a Nim implementation here: https://users.rust-lang.org/t/a-friendly-challenge-for-rust/23197

This is mostly a translation of the NIM code with a few tweaks. If you have any ideas for further optimisations, let me know.

## Benchmarks
In those benchmarks I compare the results with the Nim implementations as well as a highly optimized traditional sieve from https://primesieve.org/

Laptop with i7-4710HQ (8 threads @2.5GHz), 8GB RAM, Windows 10:

 Max_Sieve | AVG | MIN | MAX | Nim code | primesieve
 :---:|:---:|:---:|:---:|:---:|:---:|
 
 Laptop with i3-6006U (4 threads @2.0GHz), 8GB RAM, Ubuntu 18.04:
 
 Max_Sieve | AVG | MIN | MAX | Nim code | primesieve
 :---:|:---:|:---:|:---:|:---:|:---:|
