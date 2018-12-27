extern crate num;
extern crate threadpool;
extern crate num_cpus;
extern crate memsec;

mod precalculate;
mod precalculated_values;
pub mod sieve;

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::prelude::*;
    use crate::precalculate::gen_pg_parameters;
    #[test]
    fn build_constants_in_file() {
        let mut file = match File::create("precalculate.txt") {
            Ok(f) => {f},
            Err(_) => {panic!("Could not create file")},
        };
        match file.write_all(gen_pg_parameters(5).as_bytes()) {
            Ok(_) => {println!("Done 5");},
            Err(_) => {panic!("Could not write to file");},
        }
    }
}