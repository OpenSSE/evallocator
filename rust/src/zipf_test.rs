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

#[derive(Debug, StructOpt)]
#[structopt(name = "modular_zipf", about = "Sample a modular Zipf distribution")]
struct CliArgs {
    /// The number of elements in the distribution (domain)
    #[structopt(short = "n")]
    n: usize,
    /// The parameter of the Zipf distribution
    #[structopt(short = "a", long = "alpha")]
    alpha: f64,
    /// The modular parameter of the distribution
    #[structopt(short = "m", long = "mod")]
    m: usize,
    /// The number of samples
    #[structopt(short = "s", long = "samples")]
    n_samples: usize,
}
fn plot_histogram(hist: &[usize]) {
    let mut fg = Figure::new();
    fg.axes2d().boxes(0..=hist.len(), hist, &[]);
    fg.show();
}

fn main() {
    let args = CliArgs::from_args();
    println!("{:?}", args);

    let mut bins = Vec::<std::sync::atomic::AtomicUsize>::with_capacity(args.n);

    for _ in 0..args.n {
        bins.push(std::sync::atomic::AtomicUsize::new(0));
    }

    let zipf_distr = ZipfDistribution::new(args.n, args.alpha).unwrap();
    // let zipf_iter = zipf_distr.sample_iter(rng);

    (0..args.n_samples).into_par_iter().for_each(|_| {
        let x = zipf_distr.sample(&mut thread_rng());
        bins[x - 1].fetch_add(1, std::sync::atomic::Ordering::SeqCst);
    });

    let h: Vec<usize> = bins
        .into_iter()
        .map(|x| (x.into_inner()) % args.m)
        .collect();
    plot_histogram(&h);
}
