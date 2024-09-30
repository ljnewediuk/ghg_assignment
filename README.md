# ghg_assignment

Author: Levi Newediuk

This project tests whether wetland methane concentrations differ by surrounding land use. The code begins with data cleaning, focusing on formatting dates and times and sampling intervals. The rest of the code is focused on model fitting and visualization.

All code was run in R v4.3.2 on macOS Monterey v12.7.3.

Required packages and versions:

* tidyverse v2.0.0
* brms v2.20.4
* bayesplot v1.10.1
* performance v0.10.9
* bayestestR v0.13.2

Note that brms must be installed in the following steps:

* Configuring the C++ toolchain
* Installing RStan and Stan
* Installing brms

OS-specific instructions are available at https://learnb4ss.github.io/learnB4SS/articles/install-brms.html

The data required to run the code is intentionally excluded from the project.

Scripts include:

* 01 - Cleaning and visualizing the data
* 02 - Plotting the autocorrelation function to visualize autocorrelation in the time series
* 03 - Modeling parameter selection, fitting, predicting from the posterior predictive distribution, and plotting the model

