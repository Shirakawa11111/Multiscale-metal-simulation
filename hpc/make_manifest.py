"""Generate the HPC matrix manifest (name|args per line).

Matrix (57 runs, each 1 CPU, ~20-60 min with scipy/numpy FFT at 512^2):
  poly: r x strain-rate x seed  -> rate/temperature mechanism maps
  quad: r                       -> multiplication-threshold vs quench depth
  cyc:  amplitude x seed        -> cyclic recovery/ratcheting boundary
"""

lines = []

for r in (-0.35, -0.30, -0.25, -0.20, -0.15):
    for relax in (100, 400, 1600):
        for seed in (7, 11, 13):
            name = f"poly_r{r}_x{relax}_s{seed}"
            lines.append(
                f"{name}|--kind poly --r {r} --relax {relax} --seed {seed} "
                f"--strain-to 0.16 --out results_hpc/{name}")

for r in (-0.35, -0.30, -0.25, -0.20):
    name = f"quad_r{r}"
    lines.append(f"{name}|--kind quad --r {r} --relax 400 "
                 f"--strain-to 0.20 --out results_hpc/{name}")

for amp in (0.005, 0.01, 0.02, 0.03):
    for seed in (7, 11):
        name = f"cyc_a{amp}_s{seed}"
        lines.append(f"{name}|--kind cyc --amp {amp} --seed {seed} "
                     f"--out results_hpc/{name}")

with open(__file__.replace("make_manifest.py", "manifest.txt"), "w") as f:
    f.write("\n".join(lines) + "\n")
print(f"{len(lines)} runs")
