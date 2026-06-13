# Monte Carlo Percolation Model for MWCNT/PU Piezoresistive Nanocomposites

Monte Carlo simulation of the electrical percolation threshold in multi-walled carbon nanotube (MWCNT) / polyurethane (PU) nanocomposites, for both 3D bulk and 2D spray-deposited (sandwich) architectures. Accompanies the manuscript *"Predictive interphase model for PU/MWCNT piezoresistive nanocomposites"* (submitted to *npj Flexible Electronics*).

## What it does

Nanotubes are modelled as randomly positioned, randomly oriented capped cylinders (length L = 1000 nm, diameter D = 20 nm, from TEM). Two tubes are treated as electrically connected when their minimum surface-to-surface separation is within a tunnelling cutoff (default 3 nm). A configuration percolates when a single connected cluster spans the simulation cell. Percolation probability is averaged over multiple random realisations.

- **3D bulk:** percolation vs. filler loading (wt%); transition near the measured EPT of 1.8 wt%.
- **2D plane:** percolation vs. dimensionless areal density n_c·L^2; transition near the stick-percolation value 5.64.

## Requirements

- Python 3.8+
- numpy, scipy, matplotlib

Install:
```
python -m pip install numpy scipy matplotlib
```

## Usage

```
python run_percolation_full.py
```

Prints the 3D and 2D percolation sweeps to the console and saves `percolation_figure.png` and `percolation_figure.pdf` to the working directory.

## Parameters (top of the script)

| Symbol | Meaning | Default |
|--------|---------|---------|
| `L`     | nanotube length (nm)        | 1000 |
| `D`     | nanotube outer diameter (nm)| 20   |
| `DELTA` | tunnelling connection cutoff, surface-to-surface (nm) | 3 |

The random seed is fixed (`default_rng(42)`) for reproducibility.

## Citation

If you use this code, please cite the associated manuscript and this repository (Zenodo DOI to be added on release).

## License

MIT License (see LICENSE).
