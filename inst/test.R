library(readxl)
library(dplyr)
library(tidyr)
library(phytools)
library(hyperinf)
library(ggpubr)

# small test set for HyperLAU
set.seed(4)
my.tree = rtree(4)
my.df = data.frame(label = my.tree$tip.label,
                   obs = c("001", "0??", "?10", "??1"))

if(FALSE) {
  # do the HyperLAU analysis in hyperlau library, not hyperinf
c.utree = curate.uncertain.tree(my.tree, my.df)
fit = HyperLAU(c.utree$dests, c.utree$srcs)

c.utree = curate.uncertain.tree(my.tree, my.df,
                                independent.transitions = FALSE)
fit = HyperLAU(c.utree$dests, c.utree$srcs, Xbootstrap = 2)
}

# natively in hyperinf, with and without guaranteeing transition independence
fit.1 = hyperinf(my.df, my.tree)
fit.2 = hyperinf(my.df, my.tree, independent.transitions = FALSE)
ggarrange(plot_hyperinf_data(my.df, my.tree),
          plot_hyperinf(fit.1),
          plot_hyperinf(fit.2),
          nrow=1)

# uncertainty and visualisation via the bootstrap
fit.1 = hyperinf(my.df, Xboot = 5)
ggarrange(plot_hyperinf_data(my.df),
plot_hyperinf(fit.1),
plot_hyperinf(fit.1, uncertainty = TRUE),
nrow=1)

# uncertainty and visualisation via the bootstrap with HyperHMM
m1 = matrix(rep(c(0,0,1,
                 0,1,1,
                 1,1,1), 1), byrow=TRUE, ncol=3)
m2 = matrix(rep(c(0,0,1,
                  0,1,1,
                  1,1,1), 5), byrow=TRUE, ncol=3)
fit.hmm1 = hyperinf(m1, nboot = 5)
fit.hmm2 = hyperinf(m2, nboot = 5)
ggarrange(plot_hyperinf_data(m2),
          plot_hyperinf(fit.hmm1, uncertainty = TRUE),
          plot_hyperinf(fit.hmm2, uncertainty = TRUE),
          nrow=1)

# get profiles of drug resistance/susceptibility
# Supplementary Table 4 of Casali et al., https://www.nature.com/articles/ng.2878.s3
# https://static-content.springer.com/esm/art%3A10.1038%2Fng.2878/MediaObjects/41588_2014_BFng2878_MOESM35_ESM.xls
o.df = read_excel("41588_2014_BFng2878_MOESM35_ESM.xls")

# extract isolate ID and resistance profiles to our ten drugs (discarding mutation info)
# remove any incomplete profiles and recast as binary strings
col.interest = c("Isolate", "INH", "RIF", "PZA", "EMB", "STR", "AMI", "CAP", "MOX", "OFL", "PRO")
col.interest = c("Isolate", "INH", "RIF", "PZA", "EMB", "STR")
df = o.df[,which(colnames(o.df) %in% col.interest)]
missing.rows = unique(which(df == ".", arr.ind=TRUE)[,1])
#final.df = df[-missing.rows,]
final.df = df %>%
  mutate(across(-Isolate, ~ case_when(
    . == "R" ~ "1",
    . == "S" ~ "0",
    . == "." ~ "?")
  )) %>%
  unite("concatenated", -1, sep = "", remove = TRUE) 
final.df = as.data.frame(final.df)
final.df$Isolate = as.character(final.df$Isolate)
final.df$concatenated = as.character(final.df$concatenated)
my.data = final.df 
colnames(my.data)[1] = "label"

# get phylogeny linking isolates
# Supplementary Data Set 1 of Casali et al., https://www.nature.com/articles/ng.2878.s3
# https://static-content.springer.com/esm/art%3A10.1038%2Fng.2878/MediaObjects/41588_2014_BFng2878_MOESM34_ESM.txt
my.tree = read.tree("41588_2014_BFng2878_MOESM34_ESM.txt")
my.rooted.tree = my.tree

fit.tb = hyperinf(my.data, my.tree)

if(FALSE) {
prepped = curate.uncertain.tree(my.tree, my.data)
fit.tb = HyperLAU(prepped$dests, Xinitialstates = prepped$srcs, Xboot = 5)
}
plot_hyperinf(fit.tb)
