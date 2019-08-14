mod utils;
pub use crate::utils::*;

mod two_choice_alloc;

extern crate gnuplot;
// use gnuplot::{Figure, Caption, Color, LineWidth};
use gnuplot::*;

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

#[derive(Debug, Clone, Copy)]
struct AllocParams {
    pub n: usize,
    pub m: usize,
    pub max_len: usize,
    pub pad_power_2: bool,
}

fn main() {
    // let n_list: Vec<usize> = vec![
    //     1 << 21,
    //     1 << 22,
    //     1 << 23,
    //     1 << 24,
    //     1 << 25,
    //     1 << 26,
    //     1 << 27,
    // ];
    let n_list: Vec<usize> = vec![
        1 << 11,
        1 << 12,
        1 << 13,
        1 << 14,
        1 << 15,
        1 << 16,
        1 << 17,
        1 << 18,
        1 << 19,
    ];
    // let m = 1 << 10;
    // let max_len = 1 << 7;
    let iterations = 20;

    // let inputs: Vec<AllocParams> = n_list
    //     .into_iter()
    //     .map(|n| AllocParams {
    //         n: n,
    //         m: 1 << 10,
    //         max_len: 1 << 7,
    //         pad_power_2: true,
    //     })
    //     .collect();

    let inputs: Vec<AllocParams> = (24..30)
        .map(|i| {
            let n = 1 << i;
            let m = (((n as f64) / (i as f64)).ceil() as usize).next_power_of_two();
            // println!("{}",m);
            AllocParams {
                n: n,
                m: m,
                max_len: 1 << 17,
                pad_power_2: true,
            }
        })
        .collect();

    let load_stats: Vec<(AllocParams, Stats)> = inputs
        .into_iter()
        .map(|p| {
            (
                p,
                two_choice_alloc::iterated_experiment(
                    iterations,
                    p.n,
                    p.m,
                    p.max_len,
                    p.pad_power_2,
                    true,
                ),
            )
        })
        .map(|(p, results)| (p, compute_stats(results.iter().map(|x| x.max_load))))
        .collect();

    let x: Vec<usize> = load_stats.iter().map(|(p, _)| p.n).collect();
    let y: Vec<f64> = load_stats.iter().map(|(_, stat)| stat.mean).collect();
    let err: Vec<f64> = load_stats.iter().map(|(_, stat)| stat.variance).collect();
    let min_loads: Vec<usize> = load_stats.iter().map(|(_, stat)| stat.min).collect();
    let max_loads: Vec<usize> = load_stats.iter().map(|(_, stat)| stat.max).collect();
    let n_m_ratio: Vec<f64> = load_stats.iter().map(|(p, _)| (p.n as f64)/(p.m as f64)).collect();


    // println!("{:?}",x);
    println!("{:?}",n_m_ratio);
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
        ).lines_points(
            &x,
            &n_m_ratio,
            &[
                Caption("4n/m"),
                // PointSymbol('+'),
                LineWidth(1.0),
                // LineStyle(Dash),
                Color("green"),
            ],
        );
    fg.show();
}
