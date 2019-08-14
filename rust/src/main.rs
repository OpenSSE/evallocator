#![allow(dead_code)]

extern crate rand;
extern crate rayon;

use rand::prelude::*;
use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};

use indicatif::{ProgressBar, ProgressStyle};

mod utils;
pub use crate::utils::*;

fn alloc_progress<F>(n: usize, m: usize, max_len: usize, mut progress_callback: F) -> Vec<usize>
where
    F: FnMut(usize, usize),
{
    assert!(m > max_len);
    let mut buckets = vec![0; m];
    let mut remaining_elements = n;

    let mut rng = thread_rng();

    while remaining_elements != 0 {
        let l: usize = rng.gen_range(0, max_len.min(remaining_elements)) + 1;

        let n_i = l.next_power_of_two();
        let meta_buckets_counts = m / n_i;

        let b_1: usize = rng.gen_range(0, meta_buckets_counts);
        let b_2: usize = rng.gen_range(0, meta_buckets_counts);

        let slice_1: &[usize] = &buckets[(n_i * b_1)..(n_i * (b_1 + 1))];
        let count_b_1: usize = slice_1.iter().sum();

        let slice_2: &[usize] = &buckets[(n_i * b_2)..(n_i * (b_2 + 1))];
        let count_b_2: usize = slice_2.iter().sum();

        let chosen_bucket = if count_b_1 > count_b_2 { b_2 } else { b_1 };

        for count in buckets.iter_mut().skip(n_i * chosen_bucket).take(n_i) {
            *count += 1;
        }

        remaining_elements -= l;
        progress_callback(n - remaining_elements, l);
    }

    buckets
}

fn alloc(n: usize, m: usize, max_len: usize, show_progress: bool) -> Vec<usize> {
    let pb = ProgressBar::new(n as u64);
    if show_progress {
        pb.set_style(ProgressStyle::default_bar()
        .template("[{elapsed_precise}] {pos}/{len} [{bar:40.cyan/blue}] ({percent}%) | ETA: {eta_precise}")
        .progress_chars("##-"));
    }

    let progress_callback = |_, l: usize| {
        if show_progress {
            pb.inc(l as u64);
        }
    };

    let buckets = alloc_progress(n, m, max_len, progress_callback);

    if show_progress {
        pb.finish_with_message("Done!");
    }

    buckets
}

#[derive(Debug, Clone, Copy)]
struct ExperimentResult {
    size: usize,
    max_load: usize,
}

fn experiment_progress<F>(
    n: usize,
    m: usize,
    max_len: usize,
    progress_callback: F,
) -> ExperimentResult
where
    F: FnMut(usize, usize),
{
    let rand_alloc = alloc_progress(n, m, max_len, progress_callback);

    ExperimentResult {
        size: rand_alloc.iter().sum(),
        max_load: rand_alloc.iter().fold(0, |m, x| m.max(*x)),
    }
}

fn experiment(n: usize, m: usize, max_len: usize, show_progress: bool) -> ExperimentResult {
    let rand_alloc = alloc(n, m, max_len, show_progress);

    ExperimentResult {
        size: rand_alloc.iter().sum(),
        max_load: rand_alloc.iter().fold(0, |m, x| m.max(*x)),
    }
}

fn iterated_experiment(
    iterations: usize,
    n: usize,
    m: usize,
    max_len: usize,
    show_progress: bool,
) -> Vec<ExperimentResult> {
    let elements_pb = ProgressBar::new((iterations * n) as u64);
    if show_progress {
        elements_pb.set_style(ProgressStyle::default_bar()
        .template("[{elapsed_precise}] {msg} [{bar:40.cyan/blue}] ({pos}/{len} elts - {percent}%) | ETA: {eta_precise}")
        .progress_chars("##-"));
        elements_pb.set_draw_delta(100_000);
    }

    let progress_callback = |_, l: usize| {
        if show_progress {
            elements_pb.inc(l as u64);
        }
    };

    let mut iter_completed = AtomicUsize::new(0);

    if show_progress {
        elements_pb.set_position(0);
        elements_pb.set_message(&format!(
            "{}/{} iterations",
            *iter_completed.get_mut(),
            iterations
        ));
    }

    let results: Vec<ExperimentResult> = (0..iterations)
        .into_par_iter()
        .map(|_| {
            let r = experiment_progress(n, m, max_len, progress_callback);

            let previous_count = iter_completed.fetch_add(1, Ordering::SeqCst);
            if show_progress {
                elements_pb.set_message(&format!(
                    "{}/{} iterations",
                    previous_count + 1,
                    iterations
                ));
            }
            r
        })
        .collect();

    if show_progress {
        elements_pb.finish_with_message("Done!");
    }

    results
}

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

    println!(
        "{} allocation iteration with N={}, m={}, max_len={}",
        iterations, n, m, max_len
    );

    let results = iterated_experiment(iterations, n, m, max_len, true);

    let load_stats = compute_stats(results.iter().map(|x| x.max_load));
    let overhead_stats = compute_stats(results.iter().map(|x| x.size));

    println!("Load {:?}", load_stats);
    println!("Size {:?}", overhead_stats);
}
