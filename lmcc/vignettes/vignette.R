## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  install.packages("lmcc")
#  install.packages("knitr")
#  install.packages("ggplot2")
#  install.packages("ggrepel")
#  install.packages("ggthemes")
#  install.packages("tidyr")

## ---- echo=TRUE, results='asis', warning=FALSE--------------------------------
suppressPackageStartupMessages(library(lmcc))
library(knitr)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(DiceEval)
library(tidyr)

## -----------------------------------------------------------------------------
data(data_methyl)
Y = data_methyl$Y
sites = data_methyl$sites
X = data_methyl$X
colnames(X) = c("ES-cell","primary-cell")
K = nrow(Y)
N = ncol(Y)

## -----------------------------------------------------------------------------
k_star = (1:K)[-c(1,9)]
n_star = c(20:100,120:150)
Y_obs = Y
Y_obs[k_star,n_star] = NA

## -----------------------------------------------------------------------------
obj_lmcc = svd_lmcc(Y_obs,sites,X,tol_eig = 1e-6)
obj_lmcc = fit_lmcc(obj_lmcc)
Y_pred = pred_lmcc(obj_lmcc)

## -----------------------------------------------------------------------------
ind_na = is.na(Y_obs)
R2(Y[ind_na],Y_pred[ind_na])
RMSE(Y[ind_na],Y_pred[ind_na])

## -----------------------------------------------------------------------------
proc_a = 2
proc_b = 10

## ----fig.width=5.5------------------------------------------------------------
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


