# An R Shiny App for 'Predicting targets and costs for feral-cat reduction on large islands using stochastic population models'

This is a simple interactive app for exploring the effects of cat culling on feral cat populations through time. 

The model is described in the original paper by [Kathryn Venning et al. (2021)](https://conbio.onlinelibrary.wiley.com/doi/10.1111/csp2.448) with the source code [located here](https://github.com/KathrynVenning/FeralCatEradication). 

The simulation shown here assumes high harvest for first 2 years, and constant proportional harvest in remaining years. 

In the app, slider bars control: 
- initial culling proportion in the first two years
- constant proportion culled in subsequent years
- the number of replicate simulations to view stochasticity
- the time horizon in years (with yearly timesteps)

## Installation

After installing R, install the following packages with the command:

```r
install.packages(c("shiny", "tidyverse"))
```

Then run the `run_app.R` script. 

