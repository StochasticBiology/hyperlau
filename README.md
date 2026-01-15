# hyperlau
R library version of HyperLAU

Very early build, little testing. Install with `remotes::install_github("StochasticBiology/hyperlau")`

`curate.uncertain.tree` curates a tree with uncertain observations on the tips (phrased as character strings containing "0", "1", "?"). `HyperLAU` runs HyperLAU, with a mandatory first argument of observations (as a numeric matrix, with 0, 1, 2 (2 corresponding to "?"). It will also take initial states, bootstrap number, and other parameters. It returns a named list where $Dynamics is the interesting content (different bootstrap resamples are concatenated).

_hyperinf_ should now include, and provide, R HyperLAU functionality, running inference natively and plotting output. Examples below.

Here's a test run

```
remotes::install_github("StochasticBiology/hyperlau")
library(ape)
library(hyperlau)

set.seed(4)
my.tree = rtree(4)
my.df = data.frame(label = my.tree$tip.label,
                   obs = c("001", "0??", "?10", "??1"))

c.utree = curate.uncertain.tree(my.tree, my.df)
fit = HyperLAU(c.utree$dests, c.utree$srcs)
fit$Dynamics

# for plotting functions, we can use the hyperinf environment
remotes::install_github("StochasticBiology/hyperinf")
library(hyperinf)

plot_hyperinf(fit)

# and we can hopefully run HyperLAU natively from hyperinf, including uncertainty
fit.inf = hyperinf(my.df, my.tree, nboot = 5)
plot_hyperinf(fit.inf)

# natively in hyperinf, with and without guaranteeing transition independence
fit.1 = hyperinf(my.df, my.tree)
fit.2 = hyperinf(my.df, my.tree, independent.transitions = FALSE)
library(ggpubr)
ggarrange(plot_hyperinf_data(my.df, my.tree),
          plot_hyperinf(fit.1),
          plot_hyperinf(fit.2),
          nrow=1)

# uncertainty and visualisation via the bootstrap
fit.1 = hyperinf(my.df, nboot = 5)
ggarrange(plot_hyperinf_data(my.df),
          plot_hyperinf(fit.1),
          plot_hyperinf(fit.1, uncertainty = TRUE),
          nrow=1)
```
