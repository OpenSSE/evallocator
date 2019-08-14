mod utils;
pub use crate::utils::*;

mod two_choice_alloc;

// fn main() {
//     let n = 1 << 30;
//     let m = 1 << 20;
//     let max_len = 1 << 17;

//     let allocation = alloc(n, m, max_len, true);

//     let total: usize = allocation.iter().sum();
//     let max_load: usize = allocation.iter().fold(0, |m, x| m.max(*x));

//     println!("Maximum load {}", max_load);
//     println!("Total (with insertion padding) {}", total);
//     println!("Padding overhead {}", (total as f64) / (n as f64) - 1.);
//     println!("Total (with bucket padding) {}", max_load * m);
//     println!(
//         "Complete overhead {}",
//         ((max_load * m) as f64) / (n as f64) - 1.
//     );
//     // println!("{:?}", allocation);
// }

fn main() {
    let n = 1 << 30;
    let m = 1 << 20;
    let max_len = 1 << 17;
    let iterations = 200;

    let results = two_choice_alloc::iterated_experiment(iterations, n, m, max_len, true, true);

    let load_stats = compute_stats(results.iter().map(|x| x.max_load));
    let overhead_stats = compute_stats(results.iter().map(|x| x.size));

    println!("Load {:?}", load_stats);
    println!("Size {:?}", overhead_stats);
}
