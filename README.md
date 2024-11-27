### CBB914 Graphical Models for Biological Data

**EEGGO**: bayesian **E**stimation of **E**ffect-**G**ene **G**uide **O**utcomes 

Last Updated: 11.26.2024 

Authors: Jason Mohabir (jtm98), Ed Ned Neddy ()

Description: Hierarchical Gaussian Mixture Model with Latent Guide Potency

The following paper describes the biological/experimental problem being addressed:

```
Gasperini M, Hill AJ, McFaline-Figueroa JL, et al. (2019) A Genome-wide Framework for
Mapping Gene Regulation via Cellular Genetic Screens. Cell 176(1-2):377-390.e19.
doi:10.1016/j.cell.2018.11.029
https://pubmed.ncbi.nlm.nih.gov/30612741/
```

---

## Project Model 

---

## Project Description 

We will be using simulated data, which is in the shared directory on DCC:
- data1.txt
- true-betas1.txt

The file data1.txt is to be used as the input to the model.
The true-betas file is NOT to be used as input to the model, but rather only for evaluating the accuracy of the model’s estimates—i.e., to compare your estimated betas to the true betas.

The columns in data1.txt are as follows:
- geneID : identifier of a gene
- enhancerID : a putative enhancer being targeted by CRISPR; -1 means no targeting (control)
- guideID : a guide RNA targeting this enhancer; -1 means no targeting (control)
- cellID : the cell receiving this guide
- expression : the measured expression level for this gene in this cell

The experiment, in a nutshell is:
* Each cell receives either 0 or 1 guide; a guide ID of -1 denotes no guide in the cell
(control).
* Some guides are functional, and some are non-functional. Which guides are functional
is not known a priori.
* Each enhancer is targeted by 5 guides, but each cell only gets one guide (or no guide for
the control cells).
* The goal is to estimate the effect size (b) of each enhancer (not each guide!) on each
gene. This estimate should aggregate information across the guides targeting that
enhancer.
* A b of 1 means no effect; a b<1 means that perturbing the enhancer reduces gene
expression by a proportion given by b. Thus, Y = bX, where X is the control expression
level, b is the effect size, and Y is the treatment expression level when the enhancer is
perturbed by a functional guide. Remember that nonfunctional guides do not perturb
the expression; unfortunately, which guides are functional is not know a priori.

Your goal is to use STAN or JAGS to implement a Bayesian hierarchical model to estimate the b
values (one for each enhancer-gene pair), aggregating information across the guides targeting
an enhancer while accounting for the unknown functional status of individual guides. Run your
model on the data and use the true-betas file to evaluate your predictions.

---
