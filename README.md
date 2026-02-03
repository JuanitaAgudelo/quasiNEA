# quasiNEA ☄️
## Quasi-Analytic NEA Distribution

<p align="center">
  <img src="logo_quasiNEAs.png" alt="Quasi-analytic NEA distribution cover" width="800">
</p>

This repository contains the code and resources developed for my thesis focused on the  
**analytical and quasi-analytical study of the distribution of Near-Earth Asteroids (NEAs)**.

The work combines theoretical phase-space formulations, numerical integration, and Monte Carlo simulations  
to investigate the statistical structure of NEA orbital elements and their observational counterparts.

## 🔭 Scientific context

Near-Earth Asteroids (NEAs) constitute a dynamically complex population shaped by gravitational perturbations,
resonances, and observational selection effects. Understanding their distribution in orbital-element space
is essential for impact risk assessment, survey optimization, and population modeling.

This thesis explores a quasi-analytic formulation of the NEA distribution function, focusing on:

- The phase-space probability density \( P(X) \)
- Marginal distributions derived from analytical integration
- Numerical validation via Monte Carlo experiments
- Empirical calibration using multiple observational catalogs

---

## 📂 Repository structure:
```text
quasiNEA/
│
├── P_X/
│   ├── results/
│   │   └── Numerical results from the integration of P(X)
│   │
│   ├── marginal_distributions.ipynb
│   │   └── Numerical integration of marginal distributions derived from P(X)
│   │
│   └── MonteCarloExp.ipynb
│       └── Monte Carlo experiments for validating the quasi-analytic formulation
│
├── datos/
|   └──kernels/
│   │   └── Required kernels for spiceypy routines
│   │
│   └── Observational datasets:
│       - CNEOS fireball catalog
│       - JPL Small-Body Database (SBDB) queries
│       - FRIPON atmospheric entry catalog
│
├── orbital_elements_cmnd/
│   ├── products/
│   │   └── Fitted cumulative magnitude number distributions (CMND)
│   │
│   └── multimin_NEAs.ipynb
│       └── CMND fitting procedure and related diagnostic plots
│
├── earth_impactors/
│   ├── orbit_reconstruction.ipynb
│       └── Background time integration to reconstruct heliocentric orbits of impactors 
│
├── utils/
│   ├── Core utilities, classes, and dataclasses
│   ├── Orbital transformations and Jacobians
│   └── Numerical integration and probability tools
│
├── multimin.py
│   └── Standalone CMND optimization package
│
└── README.md
```

## 🧮 Methodology overview

The core methodology of this thesis consists of:

1. Defining a quasi-analytic probability density function \( P(X) \) in orbital-element space
2. Performing numerical integration to obtain marginal distributions
3. Validating analytical expectations using Monte Carlo sampling
4. Calibrating magnitude distributions through CMND fitting
5. Comparing theoretical predictions with observational data


---
## 🧪 Numerical experiments
Two complementary numerical approaches are implemented:

### Phase-space integration
Numerical integration of \( P(X) \) is performed to compute marginal distributions over selected
orbital parameters. Results are stored in `P_X/results`.

### Monte Carlo simulations
Monte Carlo experiments are used to sample the phase space and validate the analytical structure
of the model, enabling direct comparison between numerical and theoretical distributions.

## 📊 Data sources

This work uses publicly available observational datasets, including:

- [NASA CNEOS Fireball Dataset](https://cneos.jpl.nasa.gov/fireballs/)
- [JPL Small-Body Database](https://ssd.jpl.nasa.gov/tools/sbdb_query.html)
- [FRIPON](https://fireball.fripon.org/list_pipeline.php) 

All datasets are stored locally in the `datos/` directory and processed consistently
across all experiments.

## ⚙️ Tools & Dependencies

The code is written in Python and has been tested with the following environment:

- Python 3.11+
- [`NumPy`](https://numpy.org/)
- [`Pandas`](https://pandas.pydata.org/)
- [`SciPy`](https://scipy.org/)
- [`Matplotlib`](https://matplotlib.org/)

To install dependencies:

```bash
pip install -r requirements.txt
```

## 📜 License and citation
This repository is intended for academic and research purposes. If you use this work or part of the code, please cite the associated thesis and this repository.

