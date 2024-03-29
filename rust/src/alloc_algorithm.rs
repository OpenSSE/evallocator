use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum AllocAlgorithm {
    OneChoiceAllocation,
    BlockedOneChoiceAllocation,
    TwoChoiceAllocation,
    // MaxFlowAllocation,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct AllocParams {
    pub n: usize,
    pub m: usize,
    pub max_len: usize,
    pub overflow_max: usize,
    pub algorithm: AllocAlgorithm,
    pub pad_power_2: bool,
    pub block_size: usize,
    pub iterations: usize,
}

#[derive(Debug, Clone)]
pub struct ExperimentResult {
    pub size: usize,
    pub max_load: usize,
    pub load_modes: Vec<usize>,
    pub overflows: Vec<usize>,
}

#[derive(Debug, Clone, Serialize)]
pub struct AllocStats {
    pub parameters: AllocParams,
    pub size: crate::utils::Stats,
    pub load: crate::utils::Stats,
    pub load_modes: Vec<crate::utils::ModeStats>,
    pub overflows: Vec<crate::utils::ModeStats>,
}
