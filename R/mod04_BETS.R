## R/mod05_BETS.R
## EvoTempo: Module 05 - BETS Calculator (Bayesian Evaluation of Temporal Signal)

library(shiny)
library(bslib)
library(stringr)
library(ggplot2)
library(plotly)

# ---------------------------------------------------------
# 底层 MLE 自动提取函数
# ---------------------------------------------------------
extract_mle_from_log <- function(filepath) {
  tryCatch({
    # 读取文件，限制最大读取行数防止误传几个G的 MCMC log 导致卡死
    lines <- readLines(filepath, n = 5000, warn = FALSE)
    
    # 正则匹配 BEAST1 和 BEAST2 常见的 Path Sampling / Stepping Stone 结果输出
    # 匹配目标: "Marginal likelihood: -1234.56" 或 "marginal L = -1234.56"
    matches <- grep("marginal\\s*(likelihood|L)", lines, ignore.case = TRUE, value = TRUE)
    
    if (length(matches) > 0) {
      # 通常最后一行是最终组合结果
      last_match <- tail(matches, 1)
      # 提取负号或正号开头的浮点数 (支持科学计数法)
      num_str <- stringr::str_extract(last_match, "[-+]?[0-9]+\\.?[0-9]*([eE][-+]?[0-9]+)?")
      return(as.numeric(num_str))
    }
    return(NA)
  }, error = function(e) {
    return(NA)
  })
}

# ---------------------------------------------------------
# UI 端
# ---------------------------------------------------------
mod04_BETS_ui <- function(id) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      width = 4,
      h4("BETS Calculator"),
      p("Upload your BEAST Path Sampling / Stepping Stone output files. The Marginal Likelihood (Log MLE) will be automatically extracted.", class = "text-muted"),
      
      hr(),
      # 模型 1：真实时间
      h5(HTML(" Heterochronous Model (<em>M</em><sub>het</sub>)")),
      p("Model with actual sampling dates.", class = "text-muted", style="font-size:0.85em; margin-bottom:5px;"),
      fileInput(ns("het_file"), "Upload Output File (.log/.txt)", accept = c(".log", ".txt", ".out")),
      numericInput(ns("het_mle"), HTML("Log MLE Value (<em>M</em><sub>het</sub>):"), value = NA, step = 0.01),
      
      hr(),
      # 模型 2：无时间信号
      h5(HTML(" Isochronous Model (<em>M</em><sub>iso</sub>)")),
      p("Model with same time points.", class = "text-muted", style="font-size:0.85em; margin-bottom:5px;"),
      fileInput(ns("iso_file"), "Upload Output File (.log/.txt)", accept = c(".log", ".txt", ".out")),
      numericInput(ns("iso_mle"), HTML("Log MLE Value (<em>M</em><sub>iso</sub>):"), value = NA, step = 0.01),
      
      hr(),
      div(class = "d-flex gap-2 mt-2",
          actionButton(ns("calc_bets_btn"), "Evaluate Temporal Signal", class = "btn-primary flex-grow-1"),
          actionButton(ns("reset_bets_btn"), "Reset", class = "btn-outline-secondary")
      )
    ),
    
    mainPanel(
      width = 8,
      uiOutput(ns("bets_result_ui")),
      uiOutput(ns("bets_plot_ui")),
      
      # 参考文献与标准说明卡片
      card(
        class = "mt-4 border-light",
        card_header("Interpretation Scale (Kass & Raftery, 1995)"),
        card_body(
          p("The Bayes Factor (BF) compares the fit of the heterochronous model against the isochronous model. A positive log(BF) indicates support for the presence of temporal signal. Values can be estimated via Path Sampling, Stepping Stone, or Generalized Stepping Stone algorithms."),
          tags$ul(
            tags$li(strong("log(BF) > 5:"), " Very strong evidence for temporal signal."),
            tags$li(strong("3 < log(BF) ≤ 5:"), " Strong evidence."),
            tags$li(strong("1 < log(BF) ≤ 3:"), " Positive evidence."),
            tags$li(strong("0 < log(BF) ≤ 1:"), " Barely worth mentioning."),
            tags$li(strong("log(BF) < 0:"), " Evidence AGAINST temporal signal (Isochronous model is better).")
          )
        )
      )
    )
  )
}

# ---------------------------------------------------------
# Server 端
# ---------------------------------------------------------
mod04_BETS_server <- function(id, global_state) {
  moduleServer(id, function(input, output, session) {
    
    bets_res <- reactiveVal(NULL)
    
    observeEvent(input[["reset_bets_btn"]], {
      shinyjs::reset("het_file")
      shinyjs::reset("iso_file")    
      updateNumericInput(session, "het_mle", value = NA)
      updateNumericInput(session, "iso_mle", value = NA)
      bets_res(NULL)
      showNotification("BETS Calculator reset.", type = "message")
    })
    
    observeEvent(input[["het_file"]], {
      req(input[["het_file"]])
      val <- extract_mle_from_log(input[["het_file"]][["datapath"]])
      if (!is.na(val)) {
        updateNumericInput(session, "het_mle", value = val)
        showNotification("Heterochronous MLE automatically extracted!", type = "message")
      } else {
        showNotification("Could not find Marginal Likelihood in the file. Please enter manually.", type = "warning")
      }
    })
    
    observeEvent(input[["iso_file"]], {
      req(input[["iso_file"]])
      val <- extract_mle_from_log(input[["iso_file"]][["datapath"]])
      if (!is.na(val)) {
        updateNumericInput(session, "iso_mle", value = val)
        showNotification("Isochronous MLE automatically extracted!", type = "message")
      } else {
        showNotification("Could not find Marginal Likelihood in the file. Please enter manually.", type = "warning")
      }
    })
    
    observeEvent(input[["calc_bets_btn"]], {
      if (is.na(input[["het_mle"]]) || is.na(input[["iso_mle"]])) {
        showNotification("Please provide Log MLE values for both models.", type = "error")
        return()
      }
      
      # ln(Bayes Factor) = ln(MLE_heterochronous) - ln(MLE_isochronous)
      log_bf <- input[["het_mle"]] - input[["iso_mle"]]
      
      # 整合两个模型的值和计算结果，方便后续绘图调用
      bets_res(list(
        het = input[["het_mle"]],
        iso = input[["iso_mle"]],
        bf = log_bf
      ))
    })
    
    # 动态渲染图表与导出按钮的 UI
    output$bets_plot_ui <- renderUI({
      req(bets_res())
      card(
        card_header(
          class = "d-flex justify-content-between align-items-center bg-light",
          "Marginal Likelihood Comparison",
          downloadButton(session$ns("download_bets_pdf"), "Export Vector PDF", class = "btn btn-sm btn-outline-secondary")
        ),
        # 增加高度以容纳完整的边框和轴标签
        card_body(plotlyOutput(session$ns("bf_plot"), height = "250px"), class = "p-2")
      )
    })
    
    # 绘图逻辑：哑铃图 (引入与DRT一致的 theme_bw 风格)
    build_bets_plot <- function(res) {
      df <- data.frame(
        Model = factor(c("Isochronous", "Heterochronous"), levels = c("Isochronous", "Heterochronous")),
        MLE = c(res$iso, res$het)
      )
      
      # 动态计算 X 轴延伸范围，避免两端文字被裁切
      min_mle <- min(df$MLE)
      max_mle <- max(df$MLE)
      mle_range <- max_mle - min_mle
      if (mle_range == 0) mle_range <- 1
      
      ggplot(df) +
        # 连接两个数据点的粗线段
        geom_segment(aes(x = res$iso, xend = res$het, y = 1, yend = 1), color = "#BDC3C7", linewidth = 1.0) +
        # 绘制两个模型的点
        geom_point(aes(x = MLE, y = 1, color = Model), size = 3.5) +
        # 在点上方添加模型名称和数值标签 (使用固定的 y 坐标防止偏移)
        geom_text(aes(x = MLE, y = 1.12, label = sprintf("%s\n%.2f", Model, MLE), color = Model), fontface = "bold", size = 3) +
        # 在连接线正下方中心标注 Bayes Factor 的差值
        annotate("text", x = mean(c(res$iso, res$het)), y = 0.88, 
                 label = sprintf("log(BF) = %.2f", res$bf), 
                 color = "#E74C3C", fontface = "bold", size = 3.5) +
        scale_color_manual(values = c("Isochronous" = "#95A5A6", "Heterochronous" = "#2980B9")) +
        scale_y_continuous(limits = c(0.7, 1.3), expand = c(0, 0), breaks = NULL) +
        scale_x_continuous(limits = c(min_mle - mle_range * 0.25, max_mle + mle_range * 0.25)) +
        # 统一使用与 DRT 一致的 theme_bw，附带深色外线框
        theme_bw(base_size = 12) +
        theme(
          axis.title = element_text(face = "bold"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.position = "none", 
          axis.text.x = element_text(size = 10),
          axis.title.x = element_text(margin = margin(t = 12))
        ) +
        labs(x = "Log Marginal Likelihood (Log MLE)", y = "Bayesian Evaluation")
    }
    
    # 渲染 Plotly
    output$bf_plot <- renderPlotly({
      req(bets_res())
      p <- build_bets_plot(bets_res())
      ggplotly(p, tooltip = "none") %>%
        layout(
          hovermode = FALSE,
          # 精确调控图表四周边距，确保外线框完整渲染且不被遮挡
          margin = list(t = 20, b = 40, l = 40, r = 20) 
        ) %>%
        config(displayModeBar = FALSE)
    })
    
    # PDF 导出处理
    output$download_bets_pdf <- downloadHandler(
      filename = function() {
        paste0("EvoTempo_BETS_Dumbbell_", format(Sys.Date(), "%Y%m%d"), ".pdf")
      },
      content = function(file) {
        req(bets_res())
        p <- build_bets_plot(bets_res())
        ggsave(file, plot = p, width = 9, height = 3.5, device = "pdf")
      }
    )
    
    # 渲染结论面板
    output$bets_result_ui <- renderUI({
      req(bets_res())
      bf <- bets_res()$bf  # 提取 BF 用于结论评估
      
      if (bf > 5) {
        status_color <- "bg-primary text-white"
        status_icon <- "check-double"
        conclusion <- "Very Strong Evidence FOR Temporal Signal"
        advice <- "The dataset contains a highly significant temporal signal. It is safe to proceed with tip-dating calibration in BEAST."
      } else if (bf > 3) {
        status_color <- "bg-primary text-white"
        status_icon <- "check"
        conclusion <- "Strong Evidence FOR Temporal Signal"
        advice <- "The dataset contains a significant temporal signal. You can proceed with tip-dating calibration."
      } else if (bf > 1) {
        status_color <- "bg-warning text-dark"
        status_icon <- "circle-exclamation"
        conclusion <- "Positive (Moderate) Evidence FOR Temporal Signal"
        advice <- "Temporal signal is present but not overwhelmingly strong. Consider testing stricter clock models or verifying data quality."
      } else if (bf > 0) {
        status_color <- "bg-secondary text-white"
        status_icon <- "minus"
        conclusion <- "Weak / Inconclusive Evidence"
        advice <- "The temporal signal is barely supported. Tip-dating may lead to wide confidence intervals and unreliable clock rates."
      } else {
        status_color <- "bg-danger text-white"
        status_icon <- "xmark"
        conclusion <- "Evidence AGAINST Temporal Signal"
        advice <- "The isochronous model fits the data better. Tip-dating should NOT be used. Consider using informative clock rate priors or external node calibrations instead."
      }
      
      card(
        class = paste("mb-3", status_color),
        card_header(h4("BETS Conclusion", class="m-0")),
        card_body(
          h2(sprintf("log(Bayes Factor) = %.2f", bf), class="mb-3"),
          h4(icon(status_icon), " ", conclusion),
          p(advice, class="mt-2", style="font-size:1.1em;")
        )
      )
    })
    
  })
}