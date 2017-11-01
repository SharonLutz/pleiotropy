## Pleiotropy
The pleiotropy R package provides two approaches to formally test for pleiotropy with a single SNP (pleiotropySNP) or a region (pleiotropyGENE). These approaches depend on permuting the genetic region of interest and comparing the set of observed p-values to the set of permuted p-values in relation to the origin either using the Hausdorff metric or a cut-off based approach.

#### Installation
```
install.packages("devtools") # devtools must be installed first
install.packages("SKAT") # SKAT must be installed first

devtools::install_github("SharonLutz/pleiotropy")
```
#### Example 1
For a given SNP (i.e. X), one can test if this SNP is associated with 2 normally distributed phenotypes (i.e. Y) adjusting for one covariate (i.e. Z) in the given example dataset called dataS. The code below runs this analysis.
```
library(pleiotropy)
?pleiotropySNP # For details on how to use this function to test for pleiotropy with a SNP

data("dataS")
X <- dataS[,1]
Y <- dataS[,2:3]
Z <- dataS[,4]
pleiotropySNP(X,Y,Ydist=c("gaussian","gaussian"),Z,covariates=TRUE) 

```

#### Output 1
For this analysis, we have the following output. We can see that this SNP is jointly associated with both phenotypes.

```
$cutoffPvalue
[1] 0.0064

$hausdorffPvalue
[1] 0.0012
```

#### Example 2
For a given gene or collection of SNPs (i.e. X), one can test if this region is associated with 2 normally distributed phenotypes (i.e. Y) adjusting for one covariate (i.e. Z) in the given example dataset called dataG. The code below runs this analysis. 

```
library(pleiotropy)
?pleiotropyGENE # For details on how to use this function to test for pleiotropy with a gene

data("dataG")
X <- dataG[,1:5]
Y <- dataG[,6:7]
Z <- dataG[,8]
pleiotropyGENE(X,Y,Ydist=c("C","C"),Z,covariates=TRUE) # Adjusting for covariates Z
```

#### Output 2
For this analysis, we have the following output. We can see that this region is jointly associated with both phenotypes.

```
$cutoffPvalue
[1] 0.041

$hausdorffPvalue
[1] 0.036
```

#### Warning
As stated in the paper, the cut-off based approach is more robust and should be used over the Hausdorff based approach if there is any disagreement.

#### Reference
**Lutz SM**, Fingerlin T, John E Hokanson, Lange C. (2016) A General Approach to Testing for Pleiotropy with Rare and Common Variants. *Genetic Epidemiology*. 41(2):163-170.

