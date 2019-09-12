extern crate rand;
use rand::prelude::*;

extern crate zipf;
use zipf::ZipfDistribution;

extern crate gnuplot;
use gnuplot::*;

extern crate rayon;
use rayon::prelude::*;

fn plot_histogram(hist: &[usize]) {
    let mut fg = Figure::new();
    fg.axes2d().boxes(0..=hist.len(), hist, &[]);
    fg.show();
}

fn main() {
    let alpha = 1.2f64;
    let n_elts = 1000;
    let max_count: usize = 100_000;
    let n_samples: usize = 1_000_000_000;
    let mut bins = Vec::<std::sync::atomic::AtomicUsize>::with_capacity(n_elts);

    for _ in 0..n_elts {
        bins.push(std::sync::atomic::AtomicUsize::new(0));
    }

    let rng = thread_rng();

    let zipf_distr = ZipfDistribution::new(n_elts, alpha).unwrap();
    // let zipf_iter = zipf_distr.sample_iter(rng);

    (0..n_samples).into_par_iter().for_each(|_| {
        let x = zipf_distr.sample(&mut thread_rng());
        bins[x - 1].fetch_add(1, std::sync::atomic::Ordering::SeqCst);
    });

    let h: Vec<usize> = bins
        .into_iter()
        .map(|x| (x.into_inner()) % max_count)
        .collect();
    plot_histogram(&h);
}
