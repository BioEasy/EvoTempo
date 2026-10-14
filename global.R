## global.R
## EvoTempo: Global Environment & Dependencies

# 1. 核心运行库
library(shiny)
library(bslib)
library(promises)
library(future)
library(shinyjs)

# 2. 进化与统计分析核心算法库
library(ape)        
library(ade4)       
library(seqinr)     
library(coda)       

# 3. 数据处理与可视化库
library(dplyr)
library(ggplot2)
library(ggtree)     
library(plotly)     
library(crosstalk)  # 用于多图表联动

options(shiny.maxRequestSize = 100 * 1024^2)

# 设置异步执行计划，避免高强度计算（如置换检验）卡顿前端
plan(multisession)

# 4. 加载 R 目录下的所有子模块和底层工具函数
# 确保文件处于正确的相对路径
invisible(sapply(list.files("R", pattern = "\\.R$", full.names = TRUE), source))