## R/mod05_About.R
## EvoTempo: Module 05 - About

library(shiny)
library(bslib)

# ---------------------------------------------------------
# UI 端：About 页面布局
# ---------------------------------------------------------
mod05_About_ui <- function(id) {
  ns <- NS(id)
  
  fluidPage(
    # 限制最大宽度并居中，让阅读体验更像官方文档
    div(class = "container mt-4 mb-5", style = "max-width: 900px;",
        
        # 1. Software Introduction (About)
        card(
          class = "border-primary mb-4 shadow-sm",
          card_header(class = "bg-primary text-white font-weight-bold", " About EvoTempo"),
          card_body(
            p("EvoTempo is a comprehensive, interactive R/Shiny web application designed for temporal signal evaluation prior to phylodynamic analysis."),
            p("Bridging the gap between exploratory data analysis and strict Bayesian hypothesis testing, EvoTempo provides an all-in-one platform for:"),
            tags$ul(
              tags$li(tags$strong("Mantel Tests:"), " Assessing data collinearity and confounding factors before analysis."),
              tags$li(tags$strong("Root-to-Tip Regression:"), " Interactive visualization of global and local molecular clocks."),
              tags$li(tags$strong("Date Randomization Test:"), " One-click XML permutation and automated log parsing for strict CR1/CR2 evaluation."),
              tags$li(tags$strong("BETS Calculator:"), " Automated Bayesian Evaluation of Temporal Signal via Marginal Likelihoods.")
            )
          )
        ),
        
        # 2. Docs & Github
        card(
          class = "border-primary mb-4 shadow-sm",
          card_header(class = "bg-primary text-white font-weight-bold", " Docs & GitHub"),
          card_body(
            p("EvoTempo is an open-source project. For detailed tutorials, user manuals, and source code, please visit our repositories:"),
            div(class = "d-flex gap-2",
                tags$a(href = "https://github.com/BioEasy/EvoTempo", target = "_blank", class = "btn btn-outline-primary", icon("github"), " View on GitHub"),
                tags$a(href = "Manual.pdf", target = "_blank", class = "btn btn-outline-primary", icon("book-open"), " Read Manual")
            )
          )
        ),
        
        # 3. Core Algorithms & References
        card(
          class = "border-primary mb-4 shadow-sm",
          card_header(class = "bg-primary text-white font-weight-bold", " Core Algorithms & References"),
          card_body(
          p("If you use EvoTempo in your research, please consider citing the software along with the underlying algorithms:"),
            tags$ul(
              tags$li(tags$strong("TempEst / RTT Regression: "), "Rambaut, A., et al. (2016). Exploring the temporal structure of heterochronous sequences using TempEst. ", tags$em("Virus Evolution"), ", 2(1), vew007."),
              tags$li(tags$strong("DRT Criteria (CR1/CR2): "), "Duchêne, S., et al. (2015). The performance of the date-randomization test in phylogenetic trait analyses. ", tags$em("Molecular Biology and Evolution"), ", 32(4), 1065-1074."),
              tags$li(tags$strong("Confounding Theory & Mantel Test: "),  "Murray, G. G. R., et al. (2016). The effect of genetic structure on molecular dating and tests for temporal signal. ", tags$em("Methods in Ecology and Evolution"), ", 7, 80–89."),
              tags$li(tags$strong("Cluster-based DRT Algorithm: "), "Rieux, A., & Khatchikian, C. E. (2017). TipDatingBeast: an R package to assist the implementation of phylogenetic tip-dating tests using BEAST. ", tags$em("Molecular Ecology Resources"), ", 17(4), 608-613."),
              tags$li(tags$strong("BETS Framework: "), "Duchêne, S., et al. (2020). Bayesian evaluation of temporal signal in measurably evolving populations. ", tags$em("Molecular Biology and Evolution"), ", 37(11), 3363-3379."),
              tags$li(tags$strong("Tree Visualization: "), "Yu, G., et al. (2017). ggtree: an R package for visualization and annotation of phylogenetic trees. ", tags$em("Methods in Ecology and Evolution"), ", 8(1), 28-36.")
            )
          )
        ),
        
        # 4. Contact & License (并排布局)
        layout_columns(
          col_widths = c(6, 6),
          card(
            class = "border-primary shadow-sm",
            card_header(class = "bg-primary text-white font-weight-bold", " Contact"),
            card_body(
              p("For bug reports, please reach out to us:"),
              p(tags$strong("Email: "), tags$a(href = "mailto:raindy@fafu.edu.cn", class = "text-primary", "raindy@fafu.edu.cn (F. Gao)"), tags$br(), tags$strong("Or："), tags$a(href = "mailto:lienhui@fafu.edu.cn", class = "text-primary", "lienhui@fafu.edu.cn (E. Li)"))
            )
          ),
          card(
            class = "border-primary shadow-sm",
            card_header(class = "bg-primary text-white font-weight-bold", " License"),
            card_body(
              p("EvoTempo is released under the ", tags$strong("MIT License"), "."),
              p("You are free to use, modify, and distribute this software, provided that proper attribution is given.")
            )
          )
        )
    )
  )
}

# ---------------------------------------------------------
# Server 端
# ---------------------------------------------------------
mod05_About_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    # About 页面为纯静态展示，暂不需要 Server 逻辑
  })
}