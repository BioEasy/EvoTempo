## app.R
## EvoTempo: Main Application Entry Point

source("global.R")

# ---------------------------------------------------------
# UI 端：扁平化导航栏布局与全局 UI 定制
# ---------------------------------------------------------
ui <- page_navbar(
  title = tags$span(
    tags$img(
      src = "logo.png",
      width = "30px",
      height = "30px",
      style = "margin-right: 10px;"
    ),
    " EvoTempo"
  ),

  # 全局主色调
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#2E6A8C"
  ),

  # 全局 UI 内容放入 header，避免被识别为导航项
  header = tagList(
    shinyjs::useShinyjs(),

    tags$head(
      tags$style(HTML("
        /* ========================================= */
        /* 1. 顶部导航栏 (Navbar) 专属覆盖            */
        /* ========================================= */

        /* 非激活状态的标签文字 */
        .navbar-nav .nav-link {
          color: rgba(255, 255, 255, 0.8) !important;
        }

        /* 激活状态：白色文字与下划线 */
        .navbar-nav .nav-link.active,
        .navbar-nav .nav-item.show .nav-link {
          color: #FFFFFF !important;
          font-weight: bold !important;
          border-bottom: 2px solid #FFFFFF !important;
        }

        /* 悬停与键盘聚焦状态 */
        .navbar-nav .nav-link:hover,
        .navbar-nav .nav-link:focus {
          color: #FFFFFF !important;
          opacity: 1 !important;
        }

        /* ========================================= */
        /* 2. 页面内 Primary 主题色统一               */
        /* ========================================= */
        .btn-primary {
          background-color: #2E6A8C !important;
          border-color: #2E6A8C !important;
          color: #FFFFFF !important;
        }

        .btn-primary:hover {
          background-color: #049ac0 !important;
          color: #FFFFFF !important;
        }

        .bg-primary {
          background-color: #2E6A8C !important;
          color: #FFFFFF !important;
        }

        .text-primary {
          color: #2E6A8C !important;
        }

        .border-primary {
          border-color: #2E6A8C !important;
        }
      "))
    )
  ),

  # -------------------------------------------------------
  # 主功能模块
  # -------------------------------------------------------
  nav_panel(
    "1. Data & Mantel Test",
    mod01_Mantel_ui("mantel_module")
  ),

  nav_panel(
    "2. Root-to-tip Regression",
    mod02_RTT_ui("rtt_module")
  ),

  nav_panel(
    "3. Date Randomization Test",
    mod03_DRT_ui("drt_module")
  ),

  nav_panel(
    "4. BETS Calculator",
    mod04_BETS_ui("bets_module")
  ),

  nav_spacer(),

  nav_panel(
    "About",
    icon = icon("circle-info"),
    mod05_About_ui("about_module")
  )
)


# ---------------------------------------------------------
# Server 端：模块挂载与状态共享
# ---------------------------------------------------------
server <- function(input, output, session) {

  # 全局状态机
  # 用于在不同模块间传递树对象、采样日期和诊断结果
  global_state <- reactiveValues(
    raw_data = NULL,
    phylo_tree = NULL,
    mantel_p_value = NULL,
    has_confounding = FALSE,
    clusters = NULL,
    extracted_dates = NULL
  )

  # 挂载各模块的 Server 逻辑
  mod01_Mantel_server("mantel_module", global_state)
  mod02_RTT_server("rtt_module", global_state)
  mod03_DRT_server("drt_module", global_state)
  mod04_BETS_server("bets_module", global_state)
  mod05_About_server("about_module")
}


# ---------------------------------------------------------
# 启动应用
# ---------------------------------------------------------
shinyApp(ui, server)