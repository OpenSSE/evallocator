# Evallocator

Evaluate the performance of several allocation algorithms.
Implementation in Python and Rust exist.

> **_NOTE:_**  The Python code only implements the two choice allocation algorithm, and this implementation is incorrect.

## Usage (Rust)

You will need a functional rust installation. See the [Rust's documentation](<https://www.rust-lang.org/tools/install>) for details on how to install.

Move to the `rust` directory (with `cd rust`) before running the next commands.

To compile (in dev mode):
```
cargo build
```

To run (in release mode):
```
cargo run --release -- -c input_configuration.json -o output_file
```
This will use the JSON configuration file `input_configuration.json` and run several experiments, whose results will be written to `output_file.load.csv`, `output_file.size.csv` and `output_file.json`. The JSON file contains all the results of the experiments, while the CSV files contain, respectively, the statics of the load of the buckets and the size of the resulting storage.

The configuration file is a list of the following key-value pairs:
```json
{
    "n": 268435456,
    "m": 9586980,
    "max_len": 1000,
    "algorithm": "OneChoiceAllocation",
    "pad_power_2": true,
    "iterations": 10
  }
```
where `n` is the number of elements to insert, `m` the number of buckets, `max_len` the maximum list length, `algorithm` the allocation algorithm (for now, only the `"OneChoiceAllocation"` and `"TwoChoiceAllocation"` algorithms are supported), `pad_power_2` forces the lists' length to be powers of two (only used in the two choice allocation algorithm), and `iterations` is the number of experiments to run with these parameters.
The file [`example_config.json`](rust/example_config.json
) gives an example of such configuration file.


You can add the `-g` option to display the load results directly with gnuplot.