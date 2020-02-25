mod utils;
pub use crate::utils::*;

extern crate rayon;
use rayon::prelude::*;

mod alloc_algorithm;
use alloc_algorithm::*;

mod blocked_one_choice_alloc;
mod max_flow;
mod one_choice_alloc;
mod two_choice_alloc;

extern crate gnuplot;
// use gnuplot::{Figure, Caption, Color, LineWidth};
use gnuplot::*;

extern crate csv;
extern crate serde_json;

use std::convert::TryInto;
use std::fs::File;
use std::io;
use std::io::{BufReader, BufWriter};
use std::path::Path;

extern crate structopt;

use structopt::StructOpt;

// fn main() {
//     let n = 134217728;
//     let m = 134217728;
//     let max_len = 10000;

//     let allocation = one_choice_alloc::experiment(n, m, max_len, 9, true);

//     println!("Maximum load {}", allocation.max_load);
//     println!("Total (with insertion padding) {}", allocation.size);
//     println!("Padding overhead {}", (allocation.size as f64) / (n as f64) - 1.);
//     println!("Total (with bucket padding) {}", allocation.max_load * m);
//     println!(
//         "Complete overhead {}",
//         ((allocation.max_load * m) as f64) / (n as f64) - 1.
//     );

//     println!("{:?}", allocation.load_modes);
//     println!("{:?}", allocation.overflows);
//     // println!("{:?}", allocation);
// }

// fn main() {
//     let n = 1 << 30;
//     let m = 1 << 20;
//     let max_len = 1 << 17;
//     let iterations = 200;

//     let results = two_choice_alloc::iterated_experiment(iterations, n, m, max_len, true, true);

//     let load_stats = compute_stats(results.iter().map(|x| x.max_load));
//     let overhead_stats = compute_stats(results.iter().map(|x| x.size));

//     println!("Load {:?}", load_stats);
//     println!("Size {:?}", overhead_stats);
// }

#[derive(Debug, StructOpt)]
#[structopt(name = "example", about = "An example of StructOpt usage.")]
struct CliArgs {
    #[structopt(parse(from_os_str), short = "c", long = "config")]
    /// Path to a JSON configuration file. See "example_config.json" for an example
    config_path: std::path::PathBuf,
    #[structopt(
        parse(from_os_str),
        short = "o",
        long = "output",
        default_value = "results"
    )]
    /// Path for the output statistics of the experiments. A JSON and two CSV files (one for the load, the other for the space) will be generated
    output_path: std::path::PathBuf,
    #[structopt(short = "g", long = "gnuplot")]
    gnuplot: bool,
}

fn run_experiments_stats(inputs: &[AllocParams]) -> Vec<AllocStats> {
    let tot_iterations: usize = inputs.iter().map(|p| p.iterations).sum();
    let tot_elements: usize = inputs.iter().map(|p| p.iterations * p.n).sum();

    let iter_completed = std::sync::atomic::AtomicUsize::new(0);

    let pb = indicatif::ProgressBar::new(tot_elements as u64);

    pb.set_style(indicatif::ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] {msg} [{bar:40.cyan/blue}] ({pos}/{len} elts - {percent}%) | ETA: {eta_precise}")
        .progress_chars("##-"));
    pb.set_message(&format!("{}/{} iterations", 0, tot_iterations));
    pb.enable_steady_tick(1000);

    let iteration_progress_callback = |n: usize| {
        let previous_count = iter_completed.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
        pb.set_message(&format!(
            "{}/{} iterations",
            previous_count + 1,
            tot_iterations
        ));
        pb.inc(n as u64);
    };

    let results = inputs
        .into_par_iter()
        .map(|p| {
            (
                p,
                match p.algorithm {
                    AllocAlgorithm::OneChoiceAllocation => one_choice_alloc::iterated_experiment(
                        p.iterations,
                        p.n,
                        p.m,
                        p.max_len,
                        p.overflow_max,
                        false,
                        iteration_progress_callback,
                    ),
                    AllocAlgorithm::BlockedOneChoiceAllocation => {
                        blocked_one_choice_alloc::iterated_experiment(
                            p.iterations,
                            p.n,
                            p.m,
                            p.max_len,
                            p.block_size,
                            p.overflow_max,
                            false,
                            iteration_progress_callback,
                        )
                    }
                    AllocAlgorithm::TwoChoiceAllocation => two_choice_alloc::iterated_experiment(
                        p.iterations,
                        p.n,
                        p.m,
                        p.max_len,
                        p.overflow_max,
                        p.pad_power_2,
                        false,
                        iteration_progress_callback,
                    ),
                    // AllocAlgorithm::MaxFlowAllocation => max_flow::iterated_experiment(
                    //     p.iterations,
                    //     p.n,
                    //     p.m,
                    //     p.max_len,
                    //     p.max_len,
                    //     false,
                    //     iteration_progress_callback,
                    // ),
                },
                // two_choice_alloc::iterated_experiment(
            )
        })
        .map(|(p, results)| {
            let load_stat = compute_stats(results.iter().map(|x| x.max_load));
            AllocStats {
                parameters: *p,
                size: compute_stats(results.iter().map(|x| x.size)),
                load: load_stat,
                load_modes: compute_modes_stat(
                    results.iter().map(|x| &x.load_modes),
                    load_stat.max.try_into().unwrap(),
                ),
                overflows: compute_modes_stat(results.iter().map(|x| &x.overflows), p.overflow_max),
            }
        })
        .collect();

    pb.finish();
    results
}

fn read_config_file<P: AsRef<Path>>(path: P) -> io::Result<Vec<AllocParams>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    // Read the JSON contents of the file as an instance of `User`.
    let params = serde_json::from_reader(reader)?;

    Ok(params)
}

fn plot_load_stats(stats: &[AllocStats]) {
    let x: Vec<usize> = stats.iter().map(|s| s.parameters.n).collect();
    let y: Vec<f64> = stats.iter().map(|s| s.load.mean).collect();
    let err: Vec<f64> = stats.iter().map(|s| s.load.variance).collect();
    let min_loads: Vec<u64> = stats
        .iter()
        .map(|s| s.load.min.try_into().unwrap())
        .collect();
    let max_loads: Vec<u64> = stats
        .iter()
        .map(|s| s.load.max.try_into().unwrap())
        .collect();
    let expected_max_load: Vec<f64> = stats
        .iter()
        .map(|s| 3.0 * (s.parameters.n as f64) / (s.parameters.m as f64))
        .collect();

    let mut fg = Figure::new();
    fg.axes2d()
        .set_x_log(Some(2.0))
        .lines_points(
            &x,
            &y,
            &[
                Caption("Load"),
                PointSymbol('+'),
                LineWidth(1.0),
                // LineStyle(Dash),
                Color("blue"),
            ],
        )
        .y_error_bars(
            &x,
            &y,
            &err,
            &[Caption("Variance"), LineWidth(1.0), Color("red")],
        )
        .lines_points(
            &x,
            &min_loads,
            &[
                Caption("Min"),
                PointSymbol('+'),
                LineWidth(1.0),
                LineStyle(Dash),
                Color("red"),
            ],
        )
        .lines_points(
            &x,
            &max_loads,
            &[
                Caption("Min"),
                PointSymbol('+'),
                LineWidth(1.0),
                LineStyle(Dash),
                Color("red"),
            ],
        )
        .lines_points(
            &x,
            &expected_max_load,
            &[
                Caption("Expected max load: 3n/m"),
                // PointSymbol('+'),
                LineWidth(1.0),
                // LineStyle(Dash),
                Color("green"),
            ],
        );
    fg.show();
}

fn write_load_stats_csv<P: AsRef<Path>>(stats: &[AllocStats], path: P) -> io::Result<()> {
    let f = File::create(path)?;
    let writer = BufWriter::new(f);

    let mut wtr = csv::Writer::from_writer(writer);
    stats
        .iter()
        .for_each(|x| wtr.serialize((x.parameters, x.load)).unwrap());

    wtr.flush()?;
    Ok(())
}

fn write_size_stats_csv<P: AsRef<Path>>(stats: &[AllocStats], path: P) -> io::Result<()> {
    let f = File::create(path)?;
    let writer = BufWriter::new(f);

    let mut wtr = csv::Writer::from_writer(writer);
    stats
        .iter()
        .for_each(|x| wtr.serialize((x.parameters, x.size)).unwrap());

    wtr.flush()?;
    Ok(())
}

fn write_stats_json<P: AsRef<Path>>(stats: &[AllocStats], path: P) -> io::Result<()> {
    let f_json = File::create(path)?;
    serde_json::to_writer_pretty(&f_json, &stats)?;
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = CliArgs::from_args();

    println!("{:?}", args);

    let inputs = read_config_file(args.config_path)?;

    let stats = run_experiments_stats(&inputs);

    write_load_stats_csv(&stats, args.output_path.with_extension("load.csv"))?;
    write_size_stats_csv(&stats, args.output_path.with_extension("size.csv"))?;
    write_stats_json(&stats, args.output_path.with_extension("json"))?;

    if args.gnuplot {
        plot_load_stats(&stats);
    }

    Ok(())
}
