### CBB914 Graphical Models for Biological Data

**EEGGO**: bayesian **E**stimation of **E**nhancer-**G**ene **G**uide **O**utcomes 

Last Updated: 11.26.2024 

Authors: Jason Mohabir (jtm98), Edward Moseley

We implement a Gaussian Mixture Model with Latent Guide Potency in STAN. 

The following paper describes the biological/experimental problem being addressed:

```
Gasperini M, Hill AJ, McFaline-Figueroa JL, et al. (2019) A Genome-wide Framework for
Mapping Gene Regulation via Cellular Genetic Screens. Cell 176(1-2):377-390.e19.
doi:10.1016/j.cell.2018.11.029
https://pubmed.ncbi.nlm.nih.gov/30612741/
```

---

## Project Model 

![Plate diagram of Final Model](https://github.com/JasonMohabir/eeggo/blob/main/plate_diagram_final.png?raw=true)

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

---

## Repository Content 

This repository contains:
- `jupyter` Colab notebooks used for the development of the model
- `Rstan` implementation of the baseline & final model
- `cmdstanpy` driver for running models in batch 
