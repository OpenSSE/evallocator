extern crate rand;
use rand::prelude::*;

extern crate zipf;
use zipf::ZipfDistribution;

extern crate gnuplot;
use gnuplot::*;

extern crate rayon;
use rayon::prelude::*;

extern crate structopt;

use structopt::StructOpt;

mod utils;

/// The parameters of a Zipf distribution
#[derive(Debug, Clone, Copy, StructOpt)]
struct ZipfDistributionParameter {
    /// The number of elements in the distribution (domain)
    #[structopt(short = "n")]
    n: usize,
    /// The parameter of the Zipf distribution
    #[structopt(short = "a", long = "alpha")]
    alpha: f64,
}

#[derive(Debug, Clone, Copy, StructOpt)]
struct ZipfSamplingParameter {
    #[structopt(flatten)]
    distrib_params: ZipfDistributionParameter,
    /// The number of samples
    #[structopt(short = "s", long = "samples")]
    samples: usize,
}

#[derive(Debug, StructOpt)]
struct PlotCliArgs {
    #[structopt(flatten)]
    sampling_params: ZipfSamplingParameter,
    /// The modular parameter of the distribution
    #[structopt(short = "m", long = "modulus")]
    modulus: usize,
}

#[derive(Debug, StructOpt)]
struct TruncatedZipfExpCliArgs {
    #[structopt(flatten)]
    sampling_params: ZipfSamplingParameter,
    /// The modular parameter of the distribution
    #[structopt(short = "", long = "iterations")]
    iterations: usize,
    /// The modular parameter of the distribution
    #[structopt(short = "m", long = "moduli", required = true, min_values = 1)]
    moduli: Vec<usize>,
}

#[derive(Debug, StructOpt)]
#[structopt(name = "modular_zipf", about = "Sample a modular Zipf distribution")]
enum ZipfArgs {
    #[structopt(name = "plot")]
    Plot(PlotCliArgs),
    #[structopt(name = "sample")]
    Sample(TruncatedZipfExpCliArgs),
}

fn sample_zipf(sample_params: ZipfSamplingParameter) -> Vec<usize> {
    let mut bins =
        Vec::<std::sync::atomic::AtomicUsize>::with_capacity(sample_params.distrib_params.n);

    for _ in 0..sample_params.distrib_params.n {
        bins.push(std::sync::atomic::AtomicUsize::new(0));
    }

    let zipf_distr = ZipfDistribution::new(
        sample_params.distrib_params.n,
        sample_params.distrib_params.alpha,
    )
    .unwrap();
    // let zipf_iter = zipf_distr.sample_iter(rng);

    (0..sample_params.samples).into_par_iter().for_each(|_| {
        let x = zipf_distr.sample(&mut thread_rng());
        bins[x - 1].fetch_add(1, std::sync::atomic::Ordering::SeqCst);
    });

    bins.into_iter().map(|x| x.into_inner()).collect()
}

fn plot_histogram(hist: &[usize]) {
    let mut fg = Figure::new();
    fg.axes2d().boxes(0..=hist.len(), hist, &[]);
    fg.show();
}

fn plot_stats(stats: &[(usize, utils::ModeStats)]) {
    let mut fg = Figure::new();
    let x: Vec<usize> = stats.into_iter().map(|(i, _)| *i).collect();
    let mean: Vec<f64> = stats
        .into_iter()
        .map(|(_, utils::ModeStats(_, _, mean, _))| *mean)
        .collect();
    let max: Vec<usize> = stats
        .into_iter()
        .map(|(_, utils::ModeStats(_, max, _, _))| *max)
        .collect();
    let min: Vec<usize> = stats
        .into_iter()
        .map(|(_, utils::ModeStats(min, _, _, _))| *min)
        .collect();

    fg.axes2d().lines_points(&x ,&max,
            &[
                Caption("Max"),
                PointSymbol('+'),
                LineWidth(1.0),
                LineStyle(Dash),
                Color("red"),
            ]).lines_points(&x ,&min,
            &[
                Caption("Min"),
                PointSymbol('+'),
                LineWidth(1.0),
                LineStyle(Dash),
                Color("green"),
            ]).lines_points(&x ,&mean,
            &[
                Caption("Mean"),
                PointSymbol('+'),
                LineWidth(1.0),
                Color("blue"),
            ]);
    fg.show();
}

fn multimod(x: usize, moduli: &[usize]) -> Vec<usize> {
    moduli.iter().map(|m| x % m).collect()
}

fn vec_add(vt: Vec<usize>, vx: &[usize]) -> Vec<usize> {
    if vt.len() != vx.len() {
        panic!("Incompatible vector length");
    }
    vt.iter().zip(vx.iter()).map(|(v, x)| v + x).collect()
}

fn compute_remaining(bins: &[usize], moduli: &[usize]) -> Vec<usize> {
    bins.iter().fold(vec![0usize; moduli.len()], |acc, &b| {
        vec_add(acc, &multimod(b, moduli))
    })
}

fn cli_plot(args: PlotCliArgs) {
    let bins = sample_zipf(args.sampling_params);
    let h: Vec<usize> = bins.into_iter().map(|x| (x) % args.modulus).collect();
    plot_histogram(&h);
}

fn cli_sample(args: TruncatedZipfExpCliArgs) {
    let samples: Vec<Vec<usize>> = (0..args.iterations)
        .into_iter()
        .map(|_| {
            let bins = sample_zipf(args.sampling_params);

            compute_remaining(&bins, &args.moduli)
        })
        .collect();

    // println!("{:?}", samples);

    let stats = crate::utils::compute_modes_stat(samples.iter(), args.moduli.len() - 1);

    let printable_res: Vec<(usize, utils::ModeStats)> =
        args.moduli.into_iter().zip(stats.into_iter()).collect();
    println!("{:?}", printable_res);

    plot_stats(&printable_res);
}

fn main() {
    let args = ZipfArgs::from_args();
    println!("{:?}", args);

    match args {
        ZipfArgs::Plot(plot_args) => cli_plot(plot_args),
        ZipfArgs::Sample(sample_args) => cli_sample(sample_args),
    }
}
