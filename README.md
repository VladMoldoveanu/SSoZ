# SSoZ
This is an implementation of the Twin Prime Segmented Sieve of Zakiya (SSoZ). 
You can find the thread with the discussions, as well as the paper about the algorithm and a Nim implementation here: https://users.rust-lang.org/t/a-friendly-challenge-for-rust/23197

This is mostly a translation of the NIM code with a few tweaks. If you have any ideas for further optimisations, let me know.

## Benchmarks
In those benchmarks I compare the results with the Nim implementations as well as a highly optimized traditional sieve from https://primesieve.org/

Laptop with i7-4710HQ (8 threads @2.5GHz), 8GB RAM, Windows 10:

 Max_Sieve | AVG | MIN | MAX | Nim code | primesieve
 :---:|:---:|:---:|:---:|:---:|:---:|
 10^4|0.00025s|<0.0001s|0.001004s|0.001s|0.000s
 10^5|0.000255s|<0.0001s|0.00101s|0.001s|0.000s
 10^6|0.001147s|0.001s|0.002002s|0.002s|0.001s
 10^7|0.0015s|0.001s|0.002003s|0.003s|0.002s
 10^8|0.005926s|0.00098s|0.013s|0.010s|0.007s
 10^9|0.053918s|0.04997s|0.066961s|0.061s|0.065s
 10^10|0.638534s|0.628641s|0.651627s|0.653s|0.697s
 10^11|6.655793s|6.619218s|6.670174s|6.517s|8.141s
 10^12|88.983951s|88.409469s|89.697746s|88.236s|95.763s
 
 Laptop with i3-6006U (4 threads @2.0GHz), 8GB RAM, Ubuntu 18.04:
 
 Max_Sieve | AVG | MIN | MAX | Nim code | primesieve
 :---:|:---:|:---:|:---:|:---:|:---:|
 10^4|0.000098s|0.000047s|0.002594s|0.000s|0.000s
 10^5|0.000119s|0.000065s|0.002627s|0.000s|0.000s
 10^6|0.000625s|0.000416s|0.004324s|0.001s|0.000s
 10^7|0.002595s|0.001663s|0.018036s|0.006s|0.002s
 10^8|0.013957s|0.013104s|0.019407s|0.016s|0.013s
 10^9|0.145579s|0.143175s|0.152799s|0.143s|0.141s
 10^10|1.35965s|1.345758s|1.381346s|1.302s|1.635s
 10^11|15.38504s|15.330496s|15.472683s|14.187s|19.094
