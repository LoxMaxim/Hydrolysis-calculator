# Glycylglycine (GG) Kinetics Model — Shiny Application 🧪

This repository contains an interactive **R Shiny application** for analyzing the **hydrolysis** and **cyclization kinetics** of **glycylglycine (GG)** into **cyclo-glycylglycine (cGG)** and **glycine (G)**. The application performs nonlinear kinetic fitting, optional robust weighting, glycine-only exponential fitting, and exports structured datasets.

---

## Model Overview ⚛️

The kinetic scheme modeled is:

* **GG → G** with rate constant $k_1$
* **GG → cGG** with rate constant $k_2$
* **$k_{sum} = k_1 + k_2$**

Background baselines for cGG and G are fitted. All rate constants are **exponentiated** internally to maintain positivity during optimization.


---

## Features ✨

* Interactive text-based data input
* Nonlinear model fitting using **nls.lm** (from the **minpack.lm** package)
* Optional **residual weighting**
* Residual diagnostics per species
* Glycine-only exponential fit for comparison
* Export of parameters, fitted data, and high-resolution glycine curves

---

## Bibliography 📚

<details>
<summary>Citations for the R packages</summary>

### shiny
Chang W, Cheng J, Allaire J, Sievert C, Schloerke B, Aden-Buie G, Xie Y, Allen J, McPherson J, Dipert A, Borges B (2025). *shiny: Web Application Framework for R*. R package version 1.11.1.9001, https://github.com/rstudio/shiny.

### minpack.lm
Elzhov TV, Mullen KM, Spiess A-N, Bolker B (2023). *minpack.lm: R Interface to the Levenberg-Marquardt Nonlinear Least-Squares Algorithm Found in MINPACK, Plus Support for Bounds*. R package version 1.2-4, https://CRAN.R-project.org/package=minpack.lm.
</details>
  
## License

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
