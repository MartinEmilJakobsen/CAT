# Structure Learning for Directed Trees  {.tabset}
Repository for implementations, simulations and illustrations for "Structure Learning for Directed Trees" @ https://arxiv.org/abs/2108.08871. 


## Implementation of CAT 

The function `CAT` estimates a minimum edge weight directed tree using either a Gaussian or Entropy score function.
See the paper https://arxiv.org/abs/2108.08871 for further details.

### Requirements
To run the CAT implementation the following R-package dependencies are needed: `tidyverse` `magrittr` `mgcv`  `furrr`  `Matrix` `RBGL`.  
In your R-session run the following code chunk to setup and initialize the CAT implementation:

```R
install.packages("pacman")
pacman::p_load(tidyverse,magrittr,mgcv,furrr,Matrix,RBGL,devtools)
source_url("https://https://raw.githubusercontent.com/MartinEmilJakobsen/CAT/main/CAT.R")
```

### Usage

>`CAT(data, opt.args)`

##### Arguments
* `data`: a data.frame object with row-wise observations and columns representing system variables.
* `noise`: takes values "G" for Gaussian score and "NA" for entropy score. Default is "G".
* `workers`: number of threads to run the conditional expectation fitting on. Default is 1.
* `numbasisfnct`: number of basis functions for spline regression. Default is NULL which uses default values from the mgcv-package.
* `pvalCutoff`: a pruning pval cutoff threshold (between zero and 1). Detault is 1, which means no pruining.
* `crossfit`: TRUE or FALSE. Specifies wheteher edge weights are to be crossfitted or not. Default is FALSE for no crossfitting.


## Implementation of Testing Procedure

The function `HypothesisTest` allows the user to test hypotheses about whether the underlying causal directed tree contains certain substructure restrictions. See the paper https://arxiv.org/abs/2108.08871 for further details.

### Requirements
To run our hypothesis testing procedure you need to have both Python3 and R installed.

Python3 needs to have the access to networkx library.
Run the following command in your terminal
```bash
pip install networkx
```
To run our implementation the following R-package dependencies are needed: `tidyverse` `magrittr` `mgcv`  `furrr`  `Matrix` `RBGL`.  
In your R-session run the following code chunk to setup and initialize and initialize the hypothesis testing implementation:

```R
install.packages("pacman")
pacman::p_load(tidyverse,magrittr,mgcv,furrr,Matrix,devtools)
source_url("https://https://raw.githubusercontent.com/MartinEmilJakobsen/CAT/main/CAT.R")
```

### Usage

>`HypothesisTest(data, hypInput, opt.args)`

##### Arguments
* `data`: data.frame object with row-wise observations and columns representing system variables.
* `hypInput`: data.frame object containing the substructure hypothesis to test. See format below.
* `level`: sets sigfinicance level of the test. Default is 0.05
* `queryscheme`: type of query scheme method. Takes values "asymptotic" or "exact". Default is "asymptotic".
* `numbasisfnct`: number of basis functions for spline regression. Default is NULL which uses default values from the mgcv-package.

##### Value
0/1 output. 1 implies that the null-hypothesis is rejected.

### Example

`hypInput` should be a data.frame object with three columns `hypInput  = data.frame(from=*,to=*,hyp=*)`
which contains row-wise simple substructure restrictions. The collections of the simple hypotheses is then tested as a combined substructure hypothesis.

There is three types of simple hypotheses:

* Present edge hypothesis: hypothesis that the edge from V1 to V2 is present in the causal graph. 
    * Format: `from=c("V1"),to=c("V2"),hyp=c("1")`
* Missing edge hypothesis: hypothesis that the edge from V3 to V5 is not present in the causal graph. 
    * Format: `from=c("V3"),to=c("V5"),hyp=c("-1")`
* Root node hypothesis: Hypothesis that a node V1 is a root node in the causal graph. 
    * Format: `from=c("V1"),to=c("V1"),hyp=c("r")`
    

Suppose, for example, that we have observational data of a 6 dimensional system
```R
> data
# A tibble: 500 x 6
       V1     V2    V3     V4    V5     V6
    <dbl>  <dbl> <dbl>  <dbl> <dbl>  <dbl>
 1 -1.06  -0.729 0.381 -0.682 0.452  0.470
 2 -0.880 -0.646 0.617 -0.560 0.404 -0.392
 3  0.567 -1.96  0.481 -0.513 0.653  0.578
 4  2.25  -0.917 0.394 -0.524 0.951  0.339
 5 -0.561 -1.84  0.306 -0.775 0.336  0.575
 6  0.170 -1.72  0.474 -0.680 0.803  0.163
 7  1.38  -0.873 0.628 -0.707 0.638 -0.324
 8 -0.115 -1.61  0.606 -0.403 0.722  0.440
 9  1.40  -0.764 0.331 -0.923 0.354  1.08 
10  1.55  -1.20  0.533 -0.824 0.500  0.252
# ... with 490 more rows
```
    
We can test the substructure hypothesis that the causal graph satisfies all of the above simple substructure restrictions by combining them in a single `hypInput` and running the `HypothesisTest` function:

```R
hypInput = data.frame(from = c("V1","V3","V1"),
                      to   = c("V2","V5","V1"),
                      hyp  = c("1","-1","r"))
HypothesisTest(data,hypInput)
```

## Replicate Experiments in Paper

This repository also contains all R scripts for conducting the experiments referenced in the paper. In /Data the experiment data used in the paper is present. Thus, to simply rerun analysis and illustration generation on this data start from step 2) below and ignore comments about data-pointers.

If you want to replicate all experiments then start from step 1) and do not ignore comments about data-pointer as the subsequent analysis and illustration generation has to point to the newly generated data.


### Requirements
To run all our experiments and analysis script run the following code chunk to install the R-package dependencies that are needed.

```R
install.packages("pacman")
pacman::p_load(tidyverse,facetscales,stringr,ggpubr,igraph,gridExtra,DescTools,IndepTest,mgcv,furrr,RBGL,pcalg,SID,Matrix,readxl,xtable,latex2exp,glmnet,mboost,kableExtra,Rcpp,devtools)
install_url('https://cran.r-project.org/src/contrib/Archive/CAM/CAM_1.0.tar.gz')
install.packages("https://www.bnlearn.com/releases/bnlearn_latest.tar.gz", repos = NULL, type = "source")
devtools::install_github('rlebron-bioinfo/gnlearn')
```


### 1) Data generation 


To generate the experiment data: Run Script.R or manually run the following scripts. 

1. **CAMvsTree_GP_Gaussian.R**
    * *This script runs the Gaussian Process Gaussian Noise Experiment*
2. **CAMvsTree_GP_Gaussian_CamScores.R**
    * *This script runs the Gaussian Process Gaussian Noise Experiment but CAT uses the CAM scores*
3. **CAMvsTree_GP_Gaussian_Crossfit.R**
    * *This script runs the Gaussian Process Gaussian Noise Experiment with the additional method of CAT.G with cross-fitted edge weights*    
4. **CAMvsTree_GP_NonGaussian.R**
    * *This script runs the Gaussian Process Non-Gaussian Noise Experiment*
5. **IdentifiabilityConstantBi.R**
    * *This scripts runs the bivariate identifiability gap experiment*
6. **IdentifiabilityConstantMulti.R**
    * *This scripts runs the multivariate identifiability gap experiment*
6. **Testing_Multivariate.R**
    * *This scripts runs the multivariate hypothesis testing experiment*

###

### 2) Data Analysis and Generate Illustrations

All analysis and illustrations on experiment data, excluding hypothesis testing data, reported in the paper is generated by the following script

1. **AnalysisOfData.R**

All analysis and illustrations on hypothesis testing experiment data reported in the paper is generated by the following script

2. **AnalysisOfData_Testing.R**

Illustrations of the Gaussian process samplepaths are generated by the script

3. **SamplePathIllustrations.R**

Empirical application experiment of the Sachs et al. (2005) data set.

4. **SachsData.R**


Note that data-pointers in these scripts have to be adjusted to if used on new experiment data.


### Version control
The folder CAM_1.0 contains the R-package CAM used.
