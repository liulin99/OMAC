"""Regenerate the 7,000-point OMAC switching AR dataset using NumPy."""
import argparse
import csv
from pathlib import Path

import numpy as np


def generate(output_path):
    means = {1: 0.0, 2: 3.0, 3: -2.0}
    ar = {1: 0.6, 2: 0.3, 3: 0.8}
    schedule = (1, 2, 1, 3, 2, 3, 1)
    noise_sd = 0.5
    segment_length = 1000
    rng = np.random.default_rng(42)
    previous = 0.0
    step = 0

    with output_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("time", "observation", "mode", "segment"))
        for segment, mode in enumerate(schedule, start=1):
            mean = means[mode]
            phi = ar[mode]
            for _ in range(segment_length):
                noise = noise_sd * rng.standard_normal()
                if step == 0:
                    value = mean + noise / np.sqrt(1 - phi**2)
                else:
                    value = mean + phi * (previous - mean) + noise
                step += 1
                writer.writerow((step, f"{value:.10f}", mode, segment))
                previous = value


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", nargs="?", type=Path,
                        default=Path("three_mode_switching_dataset.csv"),
                        help="CSV path to create (must not already exist)")
    generate(parser.parse_args().output)
