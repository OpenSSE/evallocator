use serde::Serialize;

pub fn msb(x: u64) -> u8 {
    let bits = 64u8;
    if x != 0 {
        (bits - 1 - (x.leading_zeros() as u8))
    } else {
        0u8
    }
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct Stats {
    pub count: usize,
    pub max: usize,
    pub min: usize,
    pub mean: f64,
    pub variance: f64,
}

pub fn compute_stats<I>(iter: I) -> Stats
where
    I: Iterator<Item = usize>,
{
    let mut sum: u128 = 0;
    let mut sum_square: u128 = 0;
    let mut count: usize = 0;
    let mut min = usize::max_value();
    let mut max = 0usize;

    for v in iter {
        sum += v as u128;
        sum_square += (v as u128) * (v as u128);
        count += 1;
        min = min.min(v);
        max = max.max(v);
    }

    let mean = if count > 0 {
        (sum as f64) / (count as f64)
    } else {
        0.0
    };

    let variance = if count > 0 {
        ((sum_square as f64) / (count as f64)) - (mean * mean)
    } else {
        0.0
    };

    Stats {
        count,
        max,
        min,
        mean,
        variance,
    }
}

pub fn compute_modes<I>(iter:I, max:usize) -> Vec<usize>
where
    I: Iterator<Item = usize>
{
    let mut modes =  vec![0usize;max+1];
    for v in iter {
        modes[v] += 1;
    }
    
    modes
}

pub fn compute_modes_stat<'a,I>(iter: I, max:usize) -> Vec<f64>
where
    I: Iterator<Item = &'a Vec<usize>>,
{
    let mut modes =  vec![0.0;max+1];
    let mut count = 0;

    for m in iter {
        for i in  0..m.len() {
            modes[i] += m[i] as f64;
        }
        count += 1;
    }

    for i in  0..modes.len() {
        modes[i] /= count as f64;
    }

    modes
}


pub trait Mean<A = Self>: Sized {
    fn mean<I: Iterator<Item = A>>(iter: I) -> f64;
}

macro_rules! unsigned_integer_mean {
    (@impls  $($a:ty)*) => ($(
        impl Mean for $a {
            fn mean<I: Iterator<Item=$a>>(iter: I) -> f64 {
                let mut sum: u128 = 0;
                let mut count: usize = 0;

                for v in iter {
                    sum += v as u128;
                    count += 1;
                }
                if count > 0 {
                    (sum as f64) / (count as f64)
                } else {
                    0.0
                }
            }
        }

        impl<'a> Mean<&'a $a> for $a {
            fn mean<I: Iterator<Item=&'a $a>>(iter: I) -> f64 {
                let mut sum: u128 = 0;
                let mut count: usize = 0;

                for v in iter {
                    sum += *v as u128;
                    count += 1;
                }
                if count > 0 {
                    (sum as f64) / (count as f64)
                } else {
                    0.0
                }
            }
        }
    )*);
    ($($a:ty)*) => (
        unsigned_integer_mean!(@impls $($a)+);
    );
}

macro_rules! signed_integer_mean {
    (@impls $($a:ty)*) => ($(
        impl Mean for $a {
            fn mean<I: Iterator<Item=$a>>(iter: I) -> f64 {
                let mut sum: i128 = 0;
                let mut count: usize = 0;

                for v in iter {
                    sum += v as i128;
                    count += 1;
                }
                if count > 0 {
                    (sum as f64) / (count as f64)
                } else {
                    0.0
                }
            }
        }

        impl<'a> Mean<&'a $a> for $a {
            fn mean<I: Iterator<Item=&'a $a>>(iter: I) -> f64 {
                let mut sum: i128 = 0;
                let mut count: usize = 0;

                for v in iter {
                    sum += *v as i128;
                    count += 1;
                }
                if count > 0 {
                    (sum as f64) / (count as f64)
                } else {
                    0.0
                }
            }
        }
    )*);
    ($($a:ty)*) => (
        signed_integer_mean!(@impls $($a)+);
    );
}

macro_rules! float_mean {
    (@impls $($a:ty)*) => ($(
        impl Mean for $a {
            fn mean<I: Iterator<Item=$a>>(iter: I) -> f64 {
                let mut sum: f64 = 0.0;
                let mut count: usize = 0;

                for v in iter {
                    sum += v as f64;
                    count += 1;
                }
                if count > 0 {
                    sum / (count as f64)
                } else {
                    0.0
                }
            }
        }

        impl<'a> Mean<&'a $a> for $a {
            fn mean<I: Iterator<Item=&'a $a>>(iter: I) -> f64 {
                let mut sum: f64 = 0.0;
                let mut count: usize = 0;

                for v in iter {
                    sum += *v as f64;
                    count += 1;
                }
                if count > 0 {
                    sum / (count as f64)
                } else {
                    0.0
                }
            }
        }
    )*);
    ($($a:ty)*) => (
        float_mean!(@impls $($a)+);
    );
}

unsigned_integer_mean! { u8 u16 u32 u64 u128 usize }
signed_integer_mean! { i8 i16 i32 i64 i128 isize }
float_mean! { f32 f64 }
