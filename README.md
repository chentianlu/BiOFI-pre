BiOFI
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Introduction

<!-- badges: start -->

<!-- badges: end -->

Mounting evidences from fast-growing omics data have demonstrated the close association between metabolome and microbiome and the diverse roles of their interplay in disease and health. However, inferring reliable and driving associations faces multiple statistical challenges, given the inconsistent results from different methods, the many-to-many inter-/intra-associations among features, the insufficiently integrated analysis on composition and function, and the vast statistically significant but biologically unproved candidates. 
Here we developed an R package, BiOFI, for metabolome and microbiome (Bi-omics) feature identification. It is able to screen out key/driving features (metabolites, microbes, and functions) by an integrated importance score (IIS), combining sub-scores derived from difference, correlation, abundance, and network analysis. It is able to identify microbe-function-metabolite chains and to rank association pairs within specific functions. Rich figures and tables are provided. 
In summary, our strategy and R package provide an easy and powerful tool for the identification of robust features and correlations from big dual-omics data.


## Installation

You can install the development version of BiOFI from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("chentianlu/BiOFI")
```

## Example

# `MicPathmatch`

A funciton that is used to get pathways correspondent with the microbes
in your own 16S or Metagenome data.The result includes two parts that
one contains related detailed strain information and pathways of
microbes, and the other only contains names of microbes you want to
predict and related pathways.

An example of MicPathmatch:

``` r
library(BiOFI)
data('Microbe_path_example')
data('pathlist')
Micrest=MicPathmatch(Microbe_path_example,mic_pathways=pathlist)
```

# `Meta2pathway`

A function that can transfer metabolite abundance into related pathway
data.you can get data containing metabolite abundance and pathway that
metabolites belong to.  
Note: Please use C number in KEGG database to replace the metabolite
names.For example, C00002 means ATP or Adenosine 5’-triphosphate in
KEGG. Try your best to match the C number of metabolites for predicting
pathway that metabolites belong to more accurately. if the metabolites
do not have C numbers, keep the original names.

An example of Meta2pathway:

``` r
data('Meta2pathway_example')
data('KEGG')
Metares=Meta2pathway(Meta2pathway_example,KEGG)
```

# `IIS`

A function that calculates and selects the key nodes,the formula of IIS
is: IIS= nDFS×DFS + nDS×DS+nES×ES+nAS×AS; nDFS,nDS,nES,nAS indict the
weight value of DFS,DS,ES and AS, and nDFS+nDS+nES+nAS=1; DFS means
Difference Score; DS means Degree Score; ES means Edge Score; AS means
Abundance Score.

An example of IIS:

``` r
data('Microbe_example')
data('Meta_example')
data('groupInfo')
data('confounderData')
data('MMpathway3')
IIScore=IIS(Microbe_example,Meta_example,confounderData,groupInfo,MMpath=MMpathway3)
```

# `MMfunc`

A function that selects microbe and metabolites on the same pathway
based on the data of Nodescore.The result will be saved as HTML file.

An example of MMfunc:

``` r
data('Microbe_example')
data('Meta_example')
data('groupInfo')
data('confounderData')
data('MMpathway3')
data('MMpathway')
data('pathlist')
IIScore=IIS(Microbe_example,Meta_example,confounderData,groupInfo,MMpath=MMpathway3)
meta_mic2func=MMfunc(IIScore,MM2path=MMpathway,mic_pathways=pathlist)
```

# `TarNet`

A function used to plot a network that consists of microbes and
metabolites in same pathway.This function offers a network of microbes
and metabolites in same pathways.

An example of TarNet:

``` r
data('Microbe_example')
data('Meta_example')
data('groupInfo')
data('confounderData')
data('MMpathway3')
data('MMpathway')
data('pathlist')
IIScore=IIS(Microbe_example,Meta_example,confounderData,groupInfo,MMpath=MMpathway3)
meta_mic2func=MMfunc(IIScore,MM2path=MMpathway,mic_pathways=pathlist)
TarnetPairs <- TarNet(IIScore,meta_mic2func)
```

# `IPS`

A function used to calculate the importance pair score (IPS) of a
network based on the return of the TarNet function.

An example of IPS:

``` r
data('Microbe_example')
data('Meta_example')
data('groupInfo')
data('confounderData')
data('MMpathway3')
data('MMpathway')
data('pathlist')
IIScore=IIS(Microbe_example,Meta_example,confounderData,groupInfo,MMpath=MMpathway3)
meta_mic2func=MMfunc(IIScore,MM2path=MMpathway,mic_pathways=pathlist)
TarnetPairs <- TarNet(IIScore,meta_mic2func)
IPSres <- IPS(TarnetPairs,IIScore)
```
