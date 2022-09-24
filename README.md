# GLSP-gamma: Sparse Bayesian inference on gamma-distributed observations using shape-scale inverse-gamma mixtures

This repository provides R code implementing MCMC algorithms for scaled beta (SB) prior and inverse rescaled beta (IRB), as proposed by the following paper.

Hamura, Y., Onizuka, T., Hashimoto, S. and Sugasawa, S. (2022). Sparse Bayesian inference on gamma-distributed observations using shape-scale inverse-gamma mixtures. (https://arxiv.org/abs/2203.08440)

The repository includes the following files.

* `Gamma-GLSP-MCMC.R` : Implementation of MCMC algorithms 
* `oneshot.R`: one-shot One-shot simulation study 
* `Example-gene.R`: Example to for applying the proposed method  
* `prostatedata.RData`: Gene expression dataset used in `Example-gene.R`
