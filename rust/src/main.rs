mod utils;
pub use crate::utils::*;

extern crate rayon;
use rayon::prelude::*;

mod one_choice_alloc;
mod two_choice_alloc;

extern crate gnuplot;
// use gnuplot::{Figure, Caption, Color, LineWidth};
use gnuplot::*;

use serde::{Deserialize, Serialize};
extern crate csv;
extern crate serde_json;

use std::fs::File;
use std::io;
use std::io::{BufReader, BufWriter};
use std::path::Path;

extern crate structopt;

use structopt::StructOpt;

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

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
struct AllocParams {
    pub n: usize,
    pub m: usize,
    pub max_len: usize,
    pub pad_power_2: bool,
    pub iterations: usize,
}

#[derive(Debug, Clone, Copy, Serialize)]
struct AllocStats {
    pub parameters: AllocParams,
    pub size: utils::Stats,
    pub load: utils::Stats,
}

#[derive(Debug, StructOpt)]
#[structopt(name = "example", about = "An example of StructOpt usage.")]
struct CliArgs {
    #[structopt(parse(from_os_str), short = "c", long = "config")]
    config_path: std::path::PathBuf,
    #[structopt(parse(from_os_str), short = "o", long = "output")]
    output_path: std::path::PathBuf,
    #[structopt(short = "g", long = "gnuplot")]
    gnuplot: bool,
}

fn run_experiments_stats(inputs: &Vec<AllocParams>) -> Vec<AllocStats> {
    inputs
        .into_par_iter()
        .map(|p| {
            (
                p,
                // two_choice_alloc::iterated_experiment(
                one_choice_alloc::iterated_experiment(
                    p.iterations,
                    p.n,
                    p.m,
                    p.max_len,
                    p.pad_power_2,
                    false,
                ),
            )
        })
        .map(|(p, results)| AllocStats {
            parameters: *p,
            size: compute_stats(results.iter().map(|x| x.size)),
            load: compute_stats(results.iter().map(|x| x.max_load)),
        })
        .collect()
}

fn plot_load_stats(stats: &Vec<AllocStats>) {
    let x: Vec<usize> = stats.iter().map(|s| s.parameters.n).collect();
    let y: Vec<f64> = stats.iter().map(|s| s.load.mean).collect();
    let err: Vec<f64> = stats.iter().map(|s| s.load.variance).collect();
    let min_loads: Vec<usize> = stats.iter().map(|s| s.load.min).collect();
    let max_loads: Vec<usize> = stats.iter().map(|s| s.load.max).collect();
    let expected_max_load: Vec<f64> = stats
        .iter()
        .map(|s| 4.0 * (s.parameters.n as f64) / (s.parameters.m as f64))
        .collect();

    // println!("{:?}",x);
    // println!("{:?}",n_m_ratio);

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
                Caption("Expected max load: 4n/m"),
                // PointSymbol('+'),
                LineWidth(1.0),
                // LineStyle(Dash),
                Color("green"),
            ],
        );
    fg.show();
}

fn read_config_file<P: AsRef<Path>>(path: P) -> io::Result<Vec<AllocParams>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    // Read the JSON contents of the file as an instance of `User`.
    let params = serde_json::from_reader(reader)?;

    Ok(params)
}

fn write_load_stats_csv<P: AsRef<Path>>(stats: &Vec<AllocStats>, path: P) -> io::Result<()> {
    let f = File::create(path)?;
    let writer = BufWriter::new(f);

    let mut wtr = csv::Writer::from_writer(writer);
    stats
        .iter()
        .for_each(|x| wtr.serialize((x.parameters, x.load)).unwrap());

    wtr.flush()?;
    Ok(())
}

fn write_size_stats_csv<P: AsRef<Path>>(stats: &Vec<AllocStats>, path: P) -> io::Result<()> {
    let f = File::create(path)?;
    let writer = BufWriter::new(f);

    let mut wtr = csv::Writer::from_writer(writer);
    stats
        .iter()
        .for_each(|x| wtr.serialize((x.parameters, x.size)).unwrap());

    wtr.flush()?;
    Ok(())
}

fn write_stats_json<P: AsRef<Path>>(stats: &Vec<AllocStats>, path: P) -> io::Result<()> {
    let f_json = File::create(path)?;
    serde_json::to_writer_pretty(&f_json, &stats)?;
    Ok(())
}

fn main() {
    let args = CliArgs::from_args();

    println!("{:?}", args);

    let inputs = read_config_file(args.config_path).unwrap();

    let stats = run_experiments_stats(&inputs);

    write_load_stats_csv(&stats, args.output_path.with_extension("load.csv")).unwrap();
    write_size_stats_csv(&stats, args.output_path.with_extension("size.csv")).unwrap();
    write_stats_json(&stats, args.output_path.with_extension("json")).unwrap();

    if args.gnuplot {
        plot_load_stats(&stats);
    }
}
