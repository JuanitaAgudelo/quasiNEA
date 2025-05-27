# quasiNEA ☄️
## Quasi Analytic NEA Distribution

![orbitas_asteroides_zoom](https://github.com/user-attachments/assets/cd3f7e69-742f-4516-8395-6e5fa881b60c)

This repository contains the code and resources developed for my thesis focused on the **analytical and quasi-analytical study of the distribution of Near-Earth Asteroids (NEAs)** within the Solar System.

The central aim is to investigate and model the distribution of NEAs using a combination of theoretical methods, empirical datasets, and computational optimization techniques. This work leverages data from public space databases and incorporates statistical modeling to better understand the orbital behavior and potential observational biases in NEA data.

---

## 📓 Main Notebooks: 

- 🔁 **Back-in-time integration of fireball impact data** to compare impact trajectories with modeled distributions.
- 📊 **Fitting the NEA distribution** from JPL Small-Body Database using the `multimin` optimization package.
- 🌌 **Fitting the NEOPOP model distribution** (Near-Earth Object Population Observation Program) for comparative analysis.
- 🔍 **Magnitude analysis** to identify and quantify observational biases in the distribution of orbital elements.
- 📈 **Visualization of modeled vs observed distributions** to assess model fit and potential biases.

---

## 📁 Repository Structure
```
.
├── data/                  # Raw and processed datasets (NEAs, fireballs) 
├── figures/               # Output plots and visualizations
├── src/                   # Supplementary Python scripts and utilities
├── docs/                  # Thesis notes and reference material
├── integration.ipynb      # Main analysis notebook
├── multimin_NEAs.ipynb    # Main analysis notebook
├── visualizations.ipynb   # Main analysis notebook
└── magnitude.ipynb        # Main analysis notebook
```

---

## ⚙️ Tools & Dependencies

This project uses the following Python libraries and tools:

- Python 3.9+
- [`NumPy`](https://numpy.org/)
- [`Pandas`](https://pandas.pydata.org/)
- [`SciPy`](https://scipy.org/)
- [`Matplotlib`](https://matplotlib.org/)
- [`multimin`](https://pypi.org/project/multimin/) – for numerical optimization
- Data sources:
  - [JPL Small-Body Database](https://ssd.jpl.nasa.gov/tools/sbdb_query.html)
  - [NASA CNEOS Fireball Dataset](https://cneos.jpl.nasa.gov/fireballs/)
  - NEOPOP synthetic population data 

To install dependencies:

```bash
pip install -r requirements.txt
