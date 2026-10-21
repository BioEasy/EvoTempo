## global.R
## EvoTempo: Global Environment & Dependencies

library(shiny)
library(bslib)
library(promises)
library(future)
library(shinyjs)

library(ape)        
library(ade4)       
library(seqinr)     
library(coda)       

library(dplyr)
library(ggplot2)
library(ggtree)     
library(plotly)     
library(crosstalk)  

options(shiny.maxRequestSize = 100 * 1024^2)

plan(multisession)

invisible(sapply(list.files("R", pattern = "\\.R$", full.names = TRUE), source))
