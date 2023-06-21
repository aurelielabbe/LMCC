---
title: "LMCC: a Linear Model of Coregionalization with informative Covariates"
author: "Melina Ribaud"
date: "2023-04-04"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Some real data sets}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::knitr}
editor_options: 
  chunk_output_type: console
---




## Introduction

In this vignette, we provide a small example of how to apply the functions of the package *lmcc*. The data used in this example is a subset of the methylation dataset described in the paper "Charting a dynamic DNA methylation landscape of the human genome" from Ziller et al. (2013). Methylation was obtained by whole genome bisulfite sequencing across different cell and tissue types. We focus here on 175 sites and 15 samples obtained from two cell types defining the binary covariate used in the analysis: human embryonic stem cells (ES) and primary cells. 

## Install and load package 

Here are the instructions to load and install the packages.


```r
install.packages("lmcc")
install.packages("knitr")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("ggthemes")
install.packages("tidyr")
```

```r
suppressPackageStartupMessages(library(lmcc))
library(knitr)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(DiceEval)
#> Le chargement a nécessité le package : DiceKriging
library(tidyr)
```

## Load data

The dataset is available in the package *lmcc* and can be downloaded as follows: 


```r
data(data_methyl)
Y = data_methyl$Y
sites = data_methyl$sites
X = data_methyl$X
colnames(X) = c("ES-cell","primary-cell")
K = nrow(Y)
N = ncol(Y)
```


## Simulation of missing data 

Here, we simulate missing values on samples 2 to 8 and 10 to 15 for the sites 20 to 100 and 120 to 150.  


```r
k_star = (1:K)[-c(1,9)]
n_star = c(20:100,120:150)
Y_obs = Y
Y_obs[k_star,n_star] = NA
```


## Model fitting

The LMCC model is initialized with the function *svd_lmcc*, fitted with the function  *fit_lmcc* and missing values are predicted with the function *pred_lmcc*. 



```r
obj_lmcc = svd_lmcc(Y_obs,sites,X,tol_eig = 1e-6)
obj_lmcc = fit_lmcc(obj_lmcc)
Y_pred = pred_lmcc(obj_lmcc)
```

Since missing values were simulated in this small example, one can compute the RMSE or R2 criteria using pre-defined functions.


```r
ind_na = is.na(Y_obs)
R2(Y[ind_na],Y_pred[ind_na])
#> [1] 0.8583698
RMSE(Y[ind_na],Y_pred[ind_na])
#> [1] 0.1152161
```

## Visualize results

We choose here two samples (sample 2 and sample 10) to visualize the true methylation values  and the corresponding predicted values. 


```r
proc_a = 2
proc_b = 10
```




```r
colorTitle = "black"
sizeTitle = 15
formeTitle = "bold.italic"
colorAxe = "black"
sizeAxe = 10
formeAxe = "bold"
textSize = 15
Title = "Two samples"
low = "#349be8"
high = "#cc0000"
point_size = 2 
size_point_graph = c(2,3)
x = X

df = data.table(sites = sites,
                  Y_A = Y[proc_a, ],
                  Y_A_obs = Y_pred[proc_a, ],
                  Y_B = Y[proc_b, ],
                  Y_B_obs = Y_pred[proc_b, ])
  
  
df = data.table(df %>% pivot_longer(!sites, names_to = "origine", values_to = "values"))
df = data.table(df, Type = factor(rep(c("Pred","True"),nrow(df)/2),labels = c("True","Pred")), 
                  X = as.factor(rep(rep(c(x[proc_a],x[proc_b]),each=2),nrow(df)/4)))
    
  p = ggplot(df, aes(
    x = sites,
    y = values,
    color = X,
    shape = Type,
    size = Type
  )) +
    geom_point(alpha = 1) +
    scale_shape_manual(values = c(1, 15), name=c("Values")) + 
    scale_size_manual(guide="none", values = size_point_graph) +
    guides(color = guide_legend(override.aes = list(size=point_size)),
           shape = guide_legend(override.aes = list(size=point_size))) +
    scale_color_manual(values = c(low,high),
                       name=c("Cell type"),
                       breaks=c("0", "1"),
                       labels=c("ES", "primary"))+

    labs(title = Title,
         x = "Sites",
         y = "Gaussian processes")+
    theme_minimal() +
    theme(
      text = element_text(size=textSize),
      plot.title = element_text(
        hjust = 0.5,
        color = colorTitle,
        size = sizeTitle,
        face = formeTitle
      ),
      axis.title.x = element_text(
        color = colorAxe,
        size = sizeAxe,
        face = formeAxe
      ),
      axis.title.y = element_text(
        color = colorAxe,
        size = sizeAxe,
        face = "bold"
      )
    )
  
p
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png)
