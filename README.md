# XCU

This repository contains the codes for the in-silico model that validates the hypothesis for X-chromosome upregulation. 
To generate the data, the script "dataGen.R" needs to be run. It needs the "func.R", "init.R" and "plot_theme.R" scripts in the working directory.

R version : 4.1.2
Package dependencies: tidyverse, ggthemes, compiler, rlang, ggstatsplot

The main parameters of the simulation are as follows:

Recruitment probability ranges in XaXi genes:
Upregulated Genes : 0.3-0.8
Downregulated Genes: 0.0-0.5
Unchanged Genes : 0.2-0.5

Recruitment probabilty ranges for all XaXa genes is 0.2-0.5.
number of factors near each gene is initiated at multiple levels ranging from 100 to 1000. The threshold level of recruitment for turning the gene on is set to 70. The expression probability of genes is sigmoidally dependent on the recruitment levels. Each gene has a degradation probability of 0.35. Simulations are run for 1000 time steps. 

Author: [Kishore Hari](https://github.com/askhari139/)
