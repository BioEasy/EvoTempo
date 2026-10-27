## R/mod03_DRT.R
## EvoTempo: Module 03 - Date Randomization Test (XML Gen & Log Analysis)

library(shiny)
library(bslib)
library(zip)
library(coda)
library(ggplot2)
library(plotly)

# ---------------------------------------------------------
# UI 
# ---------------------------------------------------------
mod03_DRT_ui <- function(id) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("DRT Workflow"),
      
      radioButtons(ns("drt_tab"), "Select Step:", 
                   choices = c("Step 1: XML Generator" = "gen", 
                               "Step 2: Log Analyzer" = "log"),
                   selected = "gen"),
      hr(),
      
      conditionalPanel(
        condition = sprintf("input['%s'] == 'gen'", ns("drt_tab")),
        h5("XML Generator Settings"),
        p("Upload your BEAST template. The engine will shuffle the sampling dates to generate randomized replicates.", class="text-muted", style="font-size:0.85em;"),
        
        fileInput(ns("base_xml"), "Upload BEAST XML Template", accept = ".xml"),
        numericInput(ns("drt_reps"), "Number of Randomizations", value = 10, min = 5, max = 20),
        textInput(ns("out_prefix"), "File Prefix", value = "DRT_rep"),
        
        hr(),
        checkboxInput(ns("use_cluster_drt"), tags$b("Enable Cluster-based Permutation", style="color:#2E6A8C;"), value = FALSE),
        conditionalPanel(
          condition = sprintf("input['%s'] == true", ns("use_cluster_drt")),
          card(
            class = "bg-light mb-3",
            card_body(
              style = "padding: 10px;",
              p("Required to control for temporal and genetic confounding detected in Mantel Test.", style="font-size:0.8em; margin-bottom: 10px;"),
              radioButtons(ns("cluster_source"), "Cluster Definition:",
                           choices = c("Upload Metadata CSV" = "csv",
                                       "Auto-cluster (Mclust on dates)" = "mclust"),
                           selected = "csv"),
              conditionalPanel(
                condition = sprintf("input['%s'] == 'csv'", ns("cluster_source")),
                fileInput(ns("cluster_csv"), "Upload Cluster Metadata (.csv/.tsv)", accept = c(".csv", ".tsv")),
                uiOutput(ns("cluster_col_selectors"))
              ),
              conditionalPanel(
                condition = sprintf("input['%s'] == 'mclust'", ns("cluster_source")),
                p("Mclust will be fitted to the sampling dates of taxa in the XML. Cluster count is selected by BIC over G = 1..N.",
                  style="font-size:0.8em;"),
                numericInput(ns("mclust_max_g"), "Maximum number of clusters (G):",
                             value = 9, min = 2, max = 30),
                selectInput(ns("mclust_model"), "Mclust model name:",
                            choices = c("E" = "E", "V" = "V"),
                            selected = "E"),
                checkboxInput(ns("mclust_round"), "Round dates to nearest integer before clustering", value = FALSE),
                helpText("TipDatingBeast uses model 'E' on raw dates by default.", style="font-size:0.75em;"),
                helpText("Mclust requires a sufficient number of taxa (n >= 10 recommended).", style="font-size:0.75em;"),
                helpText("If Mclust returns very few clusters, cluster-based DRT may approximate Standard DRT.",
                         style="font-size:0.75em; color:#8a6d3b;")
              )
            )
          )
        ),
        
        div(class = "d-flex gap-2 mt-2",
            downloadButton(ns("generate_zip_btn"), "Generate & Download ZIP", class = "btn-primary flex-grow-1"),
            actionButton(ns("reset_gen_btn"), "Reset", class = "btn-outline-secondary")
        )
      ),
      
      conditionalPanel(
        condition = sprintf("input['%s'] == 'log'", ns("drt_tab")),
        h5("Log Analyzer Settings"),
        p("Evaluate CR1 and CR2 criteria by comparing the original log against randomized logs.", class="text-muted", style="font-size:0.85em;"),
        
        fileInput(ns("orig_log"), "Original Log", accept = ".log"),
        fileInput(ns("rand_logs"), "Randomized Logs (Select all replicates)", multiple = TRUE, accept = ".log"),
        numericInput(ns("burnin"), "Burn-in (%)", value = 10, min = 0, max = 90),
        
        selectizeInput(ns("param_target"), "Rate Parameter (Auto-detected)",
                       choices = c("meanRate", "clock.rate", "clockRate", "ucldMean"),
                       selected = "meanRate", options = list(create = TRUE)),
        
        div(class = "d-flex gap-2 mt-2",
            actionButton(ns("analyze_logs_btn"), "Run Evaluation", class = "btn-primary flex-grow-1"),
            actionButton(ns("reset_log_btn"), "Reset", class = "btn-outline-secondary")
        )
      )
    ),
    
    mainPanel(
      width = 9,
      
      conditionalPanel(
        condition = sprintf("input['%s'] == 'gen'", ns("drt_tab")),
        card(
          class = "border-primary mb-4",
          card_header("XML Generation Engine"),
          card_body(
            h5("How it works:"),
            tags$ul(
              tags$li("The engine automatically parses your BEAST 1.x or BEAST 2.x XML file."),
              tags$li(tags$strong("Standard Mode:"), " Completely randomizes all sampling dates (sampling without replacement)."),
              tags$li(tags$strong("Cluster Mode:"), " Permutes dates in a way that avoids assigning a sequence a date from its own cluster, controlling for temporal and genetic confounding."),
              tags$li(tags$strong("Cluster source options:"), " Upload a metadata CSV (e.g., exported from Module 1), or let the module auto-cluster the sampling dates with Mclust (as in TipDatingBeast)."),
              tags$li("All modified XML files are instantly compressed into a single ZIP archive for easy transfer to your computing cluster.")
            ),
            uiOutput(ns("xml_status_msg"))
          )
        )
      ),
      
      conditionalPanel(
        condition = sprintf("input['%s'] == 'log'", ns("drt_tab")),
        uiOutput(ns("log_results_ui"))
      )
    )
  )
}

# ---------------------------------------------------------
# Server 
# ---------------------------------------------------------
mod03_DRT_server <- function(id, global_state) {
  moduleServer(id, function(input, output, session) {
    
    observeEvent(input[["reset_gen_btn"]], {
      shinyjs::reset("base_xml")
      shinyjs::reset("cluster_csv")
      updateNumericInput(session, "drt_reps", value = 10)
      updateTextInput(session, "out_prefix", value = "DRT_rep")
      updateCheckboxInput(session, "use_cluster_drt", value = FALSE)
      updateRadioButtons(session, "cluster_source", selected = "csv")
      updateNumericInput(session, "mclust_max_g", value = 9)
      updateSelectInput(session, "mclust_model", selected = "E")
      updateCheckboxInput(session, "mclust_round", value = FALSE)
      showNotification("XML Generator settings reset.", type = "message")
    })
    
    observeEvent(input[["reset_log_btn"]], {
      shinyjs::reset("orig_log")
      shinyjs::reset("rand_logs")
      log_eval_results(NULL)
      updateNumericInput(session, "burnin", value = 10)
      showNotification("Log Analyzer results cleared.", type = "message")
    })
    
    cluster_metadata <- reactive({
      if (isTRUE(input[["use_cluster_drt"]]) && identical(input[["cluster_source"]], "csv")) {
        req(input[["cluster_csv"]])
        ext <- tools::file_ext(input[["cluster_csv"]][["name"]])
        sep <- if(ext == "tsv") "\t" else ","
        return(read.csv(input[["cluster_csv"]][["datapath"]], sep = sep, stringsAsFactors = FALSE))
      }
      return(NULL)
    })
    
    output$cluster_col_selectors <- renderUI({
      df <- cluster_metadata()
      if (is.null(df)) return(NULL)
      
      cols <- names(df)
      tagList(
        selectInput(session[["ns"]]("taxon_col"), "Taxon Name Column:", choices = cols, selected = cols[1]),
        selectInput(session[["ns"]]("cluster_col"), "Cluster/Location Column:", choices = cols, selected = if(length(cols)>1) cols[2] else cols[1])
      )
    })
    
    output$generate_zip_btn <- downloadHandler(
      filename = function() {
        paste0("EvoTempo_DRT_XMLs_", format(Sys.Date(), "%Y%m%d"), ".zip")
      },
      content = function(file) {
        req(input[["base_xml"]])
        
        id_notif <- showNotification("Parsing XML and generating permutations. Please wait...", duration = NULL, type = "message")
        
        tryCatch({
          temp_dir <- file.path(tempdir(), paste0("DRT_XML_", as.numeric(Sys.time())))
          dir.create(temp_dir, showWarnings = FALSE)
          
          n_reps <- input[["drt_reps"]]
          prefix <- input[["out_prefix"]]
          use_clust <- isTRUE(input[["use_cluster_drt"]])
          clusters_aligned <- NULL
          
          if (use_clust) {
            xml_lines <- readLines(input[["base_xml"]][["datapath"]], warn = FALSE)
            is_beast2 <- any(grepl("<trait.*traitname=\"date(-forward)?\"", xml_lines))
            xml_taxa <- NULL
            xml_dates <- NULL
            
            if (is_beast2) {
              trait_line <- xml_lines[grep("<trait.*traitname=\"date(-forward)?\"", xml_lines)[1]]
              value_str <- sub(".*value=\"([^\"]+)\".*", "\\1", trait_line)
              pairs <- strsplit(value_str, ",")[[1]]
              pairs <- trimws(pairs)
              xml_taxa  <- sapply(strsplit(pairs, "="), `[`, 1)
              xml_dates <- as.numeric(sapply(strsplit(pairs, "="), `[`, 2))
            } else {
              date_idx <- grep("<date value=", xml_lines)
              xml_taxa  <- character(length(date_idx))
              xml_dates <- numeric(length(date_idx))
              for (i in seq_along(date_idx)) {
                found <- FALSE
                for (j in (date_idx[i]-1):max(1, date_idx[i]-10)) {
                  if (grepl("<taxon id=", xml_lines[j])) {
                    xml_taxa[i] <- sub(".*<taxon id=\"([^\"]+)\".*", "\\1", xml_lines[j])
                    found <- TRUE
                    break
                  }
                }
                if (!found) stop("Failed to automatically extract taxon names from BEAST 1 XML.")
                xml_dates[i] <- as.numeric(sub('.*value="([^"]+)".*', "\\1", xml_lines[date_idx[i]]))
              }
            }

            cluster_source <- input[["cluster_source"]]
            
            if (identical(cluster_source, "csv")) {
              meta_df <- cluster_metadata()
              if (is.null(meta_df)) stop("Cluster metadata is missing. Cannot proceed with Cluster-based DRT.")
              
              t_col <- input[["taxon_col"]]
              c_col <- input[["cluster_col"]]
              
              match_idx <- match(xml_taxa, meta_df[[t_col]])
              if (any(is.na(match_idx))) {
                missing_taxa <- xml_taxa[is.na(match_idx)]
                stop(paste("The following taxa from XML are missing in metadata:",
                           paste(head(missing_taxa, 3), collapse=", "), "..."))
              }
              clusters_aligned <- meta_df[[c_col]][match_idx]
              
            } else {  # cluster_source == "mclust"
              if (!require("mclust", quietly = TRUE)) {
                stop("Package 'mclust' is required for auto-clustering. Please install it: install.packages('mclust').")
              }
              
              dates_for_clust <- xml_dates
              if (isTRUE(input[["mclust_round"]])) {
                dates_for_clust <- round(dates_for_clust)
              }
              
              n_tips <- length(dates_for_clust)
              if (n_tips < 5) stop("Too few taxa (<5) for reliable Mclust clustering.")
              
              max_g <- min(as.integer(input[["mclust_max_g"]]), n_tips)
              if (max_g < 2) stop("Too few taxa for Mclust clustering.")
              
              model_name <- input[["mclust_model"]]
              
              mclust_fit <- tryCatch(
                mclust::Mclust(as.matrix(dates_for_clust), G = 1:max_g, modelNames = model_name),
                error = function(e) stop(paste("Mclust failed:", e$message))
              )
              
              n_found <- mclust_fit$G
              if (is.null(n_found) || n_found < 2) {
                stop("Mclust found only one cluster. Standard DRT is more appropriate for this dataset.")
              }
              message(sprintf("Mclust selected G = %d (model = %s)", n_found, mclust_fit$modelName))
              
              clusters_aligned <- mclust_fit$classification
      
              global_state$mclust_result <- data.frame(
                Taxon_ID = xml_taxa,
                Date     = xml_dates,
                Cluster  = clusters_aligned,
                stringsAsFactors = FALSE
              )
            }
          }
          
          if (exists("generate_randomized_xmls")) {
            generate_randomized_xmls(
              xml_path = input[["base_xml"]][["datapath"]],
              out_dir = temp_dir,
              base_name = prefix,
              nreps = n_reps,
              use_clusters = use_clust,
              clusters = clusters_aligned
            )
          } else {
            stop("XML parsing engine not found. Ensure utils_xml.R is loaded.")
          }
          
          xml_files <- list.files(temp_dir, pattern = "\\.xml$", full.names = TRUE)
          if (length(xml_files) == 0) stop("XML generation failed. No files found.")
          
          old_wd <- getwd()
          setwd(temp_dir)
          on.exit(setwd(old_wd))
          
          zip::zip(zipfile = file, files = basename(xml_files))
          
        }, error = function(e) {
          showNotification(paste("ZIP Export Error:", e$message), type = "error")
        }, finally = {
          removeNotification(id_notif)
        })
      }
    )
    
    output$xml_status_msg <- renderUI({
      if (!is.null(input[["base_xml"]])) {
        mode_txt <- if (isTRUE(input[["use_cluster_drt"]])) {
          if (identical(input[["cluster_source"]], "mclust")) "Cluster-based (Mclust)" else "Cluster-based (CSV)"
        } else {
          "Standard"
        }
        p(class = "text-primary font-weight-bold mt-3",
          paste("Template", input[["base_xml"]][["name"]], "loaded. Ready to generate", input[["drt_reps"]], mode_txt, "replicates!"))
      }
    })
    
    observeEvent(input[["orig_log"]], {
      req(input[["orig_log"]])
      tryCatch({
        header_data <- read.table(input[["orig_log"]][["datapath"]], header = TRUE, comment.char = "#",
                                  nrows = 5, check.names = FALSE, stringsAsFactors = FALSE)
        col_names <- names(header_data)
        
        target_bases <- c("meanClockRate", "meanRate", "clock.rate", "clockRate",
                          "ORCucldMean", "ucedMean", "ucld.mean", "ucldMean")
        found_param <- NULL
        
        for (tb in target_bases) {
          if (tb %in% col_names) {
            found_param <- tb
            break
          }
        }
        
        if (is.null(found_param)) {
          target_patterns <- paste0("(^|\\.|_)", gsub("\\.", "\\\\.", target_bases), "$")
          for (pat in target_patterns) {
            match_idx <- grep(pat, col_names, ignore.case = TRUE)
            if (length(match_idx) > 0) {
              found_param <- col_names[match_idx[1]]
              break
            }
          }
        }
        
        if (!is.null(found_param)) {
          updateSelectizeInput(session, "param_target", choices = col_names, selected = found_param)
          showNotification(paste("Auto-detected rate parameter:", found_param), type = "message")
        } else {
          fallback_sel <- if(length(col_names) > 1) col_names[2] else col_names[1]
          updateSelectizeInput(session, "param_target", choices = col_names, selected = fallback_sel)
          showNotification("Rate parameter not automatically found. Please select manually.", type = "warning")
        }
      }, error = function(e) {})
    })
    
    log_eval_results <- reactiveVal(NULL)
    
    observeEvent(input[["analyze_logs_btn"]], {
    req(input[["orig_log"]], input[["rand_logs"]])
    log_eval_results(NULL)
      id_notif <- showNotification("Parsing MCMC logs and computing 95% HPDs...", duration = NULL, type = "message")
      
    tryCatch({
        burnin_frac <- input[["burnin"]] / 100
        param <- input[["param_target"]]
        
        orig_data <- read.table(input[["orig_log"]][["datapath"]], header = TRUE, comment.char = "#",
                                check.names = FALSE, stringsAsFactors = FALSE)
        if (!param %in% names(orig_data)) {
          stop(sprintf("Parameter '%s' not found in Original Log.", param))
        }
        
        orig_start <- floor(nrow(orig_data) * burnin_frac) + 1
        orig_vals <- orig_data[[param]][orig_start:nrow(orig_data)]
        orig_mean <- mean(orig_vals)
        orig_hpd <- coda::HPDinterval(coda::as.mcmc(orig_vals), prob = 0.95)
        
        rand_stats <- lapply(1:nrow(input[["rand_logs"]]), function(i) {
          r_data <- read.table(input[["rand_logs"]][["datapath"]][i], header = TRUE, comment.char = "#",
                               check.names = FALSE, stringsAsFactors = FALSE)
          if (!param %in% names(r_data)) {
            stop(sprintf("Parameter '%s' not found in randomized log '%s'. All uploaded logs must contain the same parameter.", param, input[["rand_logs"]][["name"]][i]))
          }
          
          r_start <- floor(nrow(r_data) * burnin_frac) + 1
          r_vals <- r_data[[param]][r_start:nrow(r_data)]
          r_hpd <- coda::HPDinterval(coda::as.mcmc(r_vals), prob = 0.95)
          
        data.frame(
            Rep = paste("Randomized", sprintf("%02d", i)),
            Mean = mean(r_vals),
            Lower = r_hpd[1, 1],
            Upper = r_hpd[1, 2],
            Type = "Randomized",
            stringsAsFactors = FALSE
          )
        })
        
        rand_df <- do.call(rbind, rand_stats)
        if (is.null(rand_df) || nrow(rand_df) == 0) {
          stop("Failed to extract parameter from randomized logs.")
        }
        
        lower_vals <- rand_df[["Lower"]]
        upper_vals <- rand_df[["Upper"]]
        
        cr1_pass <- all(orig_mean < lower_vals | orig_mean > upper_vals)
        cr2_pass <- all(orig_hpd[1, 2] < lower_vals | orig_hpd[1, 1] > upper_vals)
        
        orig_df <- data.frame(
          Rep = "Original Dataset",
          Mean = orig_mean,
          Lower = orig_hpd[1, 1],
          Upper = orig_hpd[1, 2],
          Type = "Original",
          stringsAsFactors = FALSE
        )
        
        rand_df <- rand_df[order(rand_df[["Lower"]]), ]
        plot_df <- rbind(orig_df, rand_df)
        
        log_eval_results(list(
          plot_df = plot_df,
          cr1 = cr1_pass,
          cr2 = cr2_pass,
          param = param
        ))
        
      }, error = function(e) {
        showNotification(paste("Analysis Error:", e$message), type = "error")
      }, finally = {
        removeNotification(id_notif)
      })
    })
    
    build_drt_plot <- function(df, param) {
      df[["Rep"]] <- factor(df[["Rep"]], levels = c("Original Dataset", df[["Rep"]][df[["Type"]] == "Randomized"]))
      
      ggplot(df, aes(x = Rep, y = Mean, color = Type,
                     text = paste("Dataset:", Rep,
                                  "<br>Mean:", signif(Mean, 4),
                                  "<br>95% HPD:", signif(Lower, 4), "-", signif(Upper, 4)))) +
        geom_point(size = 3) +
        geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.3, linewidth = 1) +
        geom_hline(yintercept = df[["Lower"]][df[["Type"]] == "Original"], linetype = "dashed", color = "#2E6A8C", alpha=0.5) +
        geom_hline(yintercept = df[["Upper"]][df[["Type"]] == "Original"], linetype = "dashed", color = "#2E6A8C", alpha=0.5) +
        scale_color_manual(name = "Dataset", values = c("Original" = "#2E6A8C", "Randomized" = "#95A5A6")) +
        theme_bw(base_size = 12) +
        labs(x = "Datasets", y = "Substitution rate", color = "Dataset") +
        theme(
          legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          axis.title = element_text(face = "bold")
        )
    }
    
    output$log_results_ui <- renderUI({
      res <- log_eval_results()
      if (is.null(res)) {
        card(
          card_header("Evaluation Results: CR1 & CR2 Criteria"),
          card_body(
            h5("Awaiting Logs..."),
            p("Upload the original log and all randomized replicate logs to evaluate temporal signal strictness.", class = "text-muted")
          )
        )
      } else {
        cr1_color <- ifelse(res[["cr1"]], "success", "danger")
        cr2_color <- ifelse(res[["cr2"]], "success", "warning")
        
        card(
          class = "border-primary",
          card_header("DRT Evaluation Results"),
          card_body(
            layout_columns(
              col_widths = c(6, 6),
              card(class = paste0("bg-", cr1_color, " text-white"),
                   card_body(h5("CR1 (Standard)"),
                             p(ifelse(res[["cr1"]], "PASS: Original mean falls completely outside all randomized 95% HPDs.",
                                      "FAIL: Original mean overlaps with randomized 95% HPDs."), style="font-size:0.9em;"))
              ),
              card(class = paste0("bg-", cr2_color, " text-white"),
                   card_body(h5("CR2 (Stringent)"),
                             p(ifelse(res[["cr2"]], "PASS: Original 95% HPD shows ZERO overlap with any randomized 95% HPD.",
                                      "FAIL: Original 95% HPD overlaps with randomized 95% HPDs."), style="font-size:0.9em;"))
              )
            ),
            div(class = "d-flex justify-content-between align-items-center mt-4 mb-3",
                h5(class="text-primary font-weight-bold m-0", "95% HPD Interval Comparisons"),
                div(class = "d-flex gap-2",
                    downloadButton(session[["ns"]]("download_drt_data"), "Export Plot Data", class = "btn btn-sm btn-outline-secondary"),
                    downloadButton(session[["ns"]]("download_drt_pdf"), "Export Vector PDF", class = "btn btn-sm btn-outline-primary")
                )
            ),
            plotlyOutput(session[["ns"]]("drt_forest_plot"), height = "450px")
          )
        )
      }
    })
    
    output$download_drt_data <- downloadHandler(
      filename = function() {
        paste0("EvoTempo_DRT_Data_", format(Sys.Date(), "%Y%m%d"), ".csv")
      },
      content = function(file) {
        req(log_eval_results())
        write.csv(log_eval_results()$plot_df, file, row.names = FALSE)
      }
    )
    
    output$drt_forest_plot <- renderPlotly({
      req(log_eval_results())
      res <- log_eval_results()
      p <- build_drt_plot(res[["plot_df"]], res[["param"]])
      
      ggplotly(p, tooltip = "text") %>%
        config(displayModeBar = TRUE) %>%
        layout(
          hovermode = "closest",
          legend = list(orientation = "h", x = 0.5, y = 1.15, xanchor = "center", yanchor = "bottom")
        )
    })
    
    output$download_drt_pdf <- downloadHandler(
      filename = function() {
        paste0("EvoTempo_DRT_Plot_", format(Sys.Date(), "%Y%m%d"), ".pdf")
      },
      content = function(file) {
        req(log_eval_results())
        res <- log_eval_results()
        p <- build_drt_plot(res[["plot_df"]], res[["param"]])
        ggsave(file, plot = p, width = 10, height = 6, device = "pdf")
      }
    )
    
  })
}
