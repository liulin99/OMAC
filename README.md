# Online Mode-Adaptive Conformal Prediction for Time Series with Distribution Shift

Code and data for the OMAC experiments on Electricity, Traffic, simulated glucose and a three-mode switching autoregressive process. OMAC issues the union of all retained intervals, monitors recent misses, and fits a new predictor and threshold using a fresh post-alarm buffer.

## Reproduce the paper

From the repository root, use Python 3.12 and install the dependencies.

```sh
python -m pip install -r requirements.txt
python -m jupyterlab
```

Open each notebook and choose **Run All Cells**. The paper inputs are included, so the default experiments can run offline after installation. The notebooks write the figure panels and CSV tables in the working directory. In the synthetic notebook, the method labeled **Proposed** is OMAC.

| Result | File to run | Output |
|---|---|---|
| Electricity and Traffic, Table I and Figures 1, 2, 5 and 6 | `OMAC_real_data.ipynb` | `real_world_cp_comparison.csv` and plots |
| Synthetic results, Table II and Figures 3, 4, 7 and 8 | `OMAC_synthetic_data.ipynb` | `synthetic_overall_results.csv` and plots |
| Holm-adjusted two-sample KS values, Table III | Run both notebooks | `real_world_omac_detection_window_tests.csv` and `synthetic_detection_window_tests.csv` |
| Short-warm-up comparison | Real-data notebook, last cell | `short_warmup_cp_comparison.csv` |
| Detection delays and post-adaptation coverage | Synthetic notebook | `synthetic_detection_results.csv` and `synthetic_post_adaptation_results.csv` |

For the additional comparisons and stationary false-alarm experiments in the camera-ready discussion, run:

```sh
python run_comparisons.py --output results
```

The script writes the following files. Reference outputs from verification are included under `reference_results/`.

- `comparisons.csv` reports coverage, mean width and mean 90% interval score. For OMAC, coverage and width use the actual union, while the interval score uses its outer endpoints. The matched newest-only and refitted ACI comparisons use the original OMAC refit schedule. A second newest-only run uses its own alarms.
- `alarms.csv` records zero-based detection indices. `buffer_labels.csv` counts the physical labels in the detection and post-alarm windows of the two synthetic series.
- `stationary_null.csv` and `fitted_stationary_null.csv` report the overlapping-window sensitivity grid for a fixed 90% interval and a fitted OMAC predictor, respectively. A false detection is any alarm on a stationary stream. The random seeds are included in the script.

## Data and simulators

| Input | Observations | Contents |
|---|---:|---|
| `data/hourly_aggregates.npz` | Electricity 21,212 and Traffic 14,204 | Hourly series in the paper's scoring units, with the timestamp of each series' first observation. Electricity sums client consumption in millions of kWh. Traffic averages sensor occupancy in percentage points. |
| `three_mode_switching_dataset.csv` | 7,000 | Observations, physical modes and segment indices. |
| `glucose_simulation_hourly.csv` | 8,761 | Hourly observations for a 365-day mouse-model simulation, including hour zero, with meals, physical modes, true glucose and noisy glucose measurements. |

To download the original Electricity and Traffic data instead, install `gluonts==0.17.0` in the notebook environment and set `USE_ARCHIVED_AGGREGATES = False` in the first real-data cell before running it. This path calls `get_dataset("electricity")` and `get_dataset("traffic")` and computes the aggregates from the downloaded series.

To regenerate the switching AR series, write to a new filename:

```sh
python generate_three_mode_switching.py regenerated_three_mode.csv
```

The generator uses NumPy seed 42 and reproduces the included CSV byte for byte. To generate glucose using MATLAB, run the following from this directory:

```matlab
data = simulate_year_modes_full("generated_glucose", 365);
```

The MATLAB driver contains the physiological model functions and writes CSV, MAT and PNG outputs into the specified folder. The supplied glucose CSV can be used directly for every analysis.
