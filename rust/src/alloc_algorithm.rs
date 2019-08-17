use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum AllocAlgorithm {
    OneChoiceAllocation,
    TwoChoiceAllocation,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct AllocParams {
    pub n: usize,
    pub m: usize,
    pub max_len: usize,
    pub algorithm: AllocAlgorithm,
    pub pad_power_2: bool,
    pub iterations: usize,
}

#[derive(Debug, Clone, Copy)]
pub struct ExperimentResult {
    pub size: usize,
    pub max_load: usize,
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct AllocStats {
    pub parameters: AllocParams,
    pub size: crate::utils::Stats,
    pub load: crate::utils::Stats,
}
