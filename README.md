
# COVID-19 Modeling in Croatia 🇭🇷

This repository contains code and data for modeling the spread of COVID-19 in Croatia, developed as part of the postgraduate course *Mathematical Models in Infectious Diseases* at the National and Kapodistrian University of Athens (NKUA), during the 2024–2025 academic year.

##  Project Overview

The goal is to analyze and simulate the progression of the COVID-19 epidemic using both deterministic and stochastic models:

- Epidemic curve construction
- Estimation of the effective reproduction number (Rt)
- SEIR simulation models (baseline & intervention scenarios)
- Age-structured transmission modeling
- Stochastic chain binomial model (with MCMC)
- Change point detection in transmission rate

##  Repository Contents

| File                            | Description                                                  |
|---------------------------------|--------------------------------------------------------------|
| `Croatia_COVID19_Model_Bitinis.R` | Main R script with fully annotated code and simulations      |
| `data.xlsx`                     | Dataset used for analysis (source: ECDC)                     |


##  Tools & Packages

- `R` (main language)
- Key packages: `ggplot2`, `dplyr`, `deSolve`, `EpiEstim`, `R2OpenBUGS`, `coda`, `bayesplot`
- External software: **OpenBUGS** for MCMC simulation

##  Methodologies Implemented

- Rt estimation using **Cori et al. (2013)**
- SEIR modeling under:
  - No interventions (baseline)
  - Social distancing (reduced contact rate)
  - Age-structured population
- Stochastic modeling with Bayesian estimation and change points

##  Author

**Ilias Bitinis**  
M.Sc. student in Biostatistics & Health Data Science, Medical University of NKUA, Athens Greece||

B.Sc.  in Mathematics, University of Patras, Patras Greece||

[LinkedIn Profile](https://www.linkedin.com/in/ilias-bitinis-77b158260)

---

> For educational and research purposes only. Data source: European Centre for Disease Prevention and Control (ECDC).
