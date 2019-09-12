extern crate rand;
use rand::prelude::*;

extern crate zipf;
use zipf::ZipfDistribution;

extern crate gnuplot;
use gnuplot::*;

fn plot_histogram(hist: &[usize]) {
    let mut fg = Figure::new();
    fg.axes2d().boxes(0..=hist.len(), hist, &[]);
    fg.show();
}

fn main() {
    let alpha = 1.2f64;
    let n_elts = 1000;
    let max_count: usize = 1000;
    let n_samples: usize = 1_000_000;
    let mut bins = vec![0usize; n_elts];

    let rng = thread_rng();

    let zipf_distr = ZipfDistribution::new(n_elts, alpha).unwrap();
    let zipf_iter = zipf_distr.sample_iter(rng);

    for x in zipf_iter.take(n_samples) {
        if bins[x - 1] == max_count - 1 {
            bins[x - 1] = 0;
        } else {
            bins[x - 1] += 1;
        }
    }

    plot_histogram(&bins);
}
