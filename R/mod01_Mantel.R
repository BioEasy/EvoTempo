## R/mod01_Mantel.R
## EvoTempo: Module 01 - Data Input, Parsing & Mantel Test

library(shiny)
library(bslib)
library(DT)
library(plotly)

# ---------------------------------------------------------
# UI
# ---------------------------------------------------------
mod01_Mantel_ui <- function(id) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      width = 4,
      h4("1. Data Upload"),
      fileInput(ns("file_upload"), "Sequence / Tree File"),
      selectInput(ns("file_format"), "Format", 
                  choices = c("FASTA" = "fasta", "Newick Tree" = "newick", "NEXUS Tree" = "nexus")),
      
      h4("2. Date Source", class = "mt-4"),
      radioButtons(ns("date_source"), "Where are the dates?", 
                   choices = c("Parse from Taxon ID" = "id", "Upload Metadata (CSV)" = "csv")),
      
      conditionalPanel(
        condition = sprintf("input['%s'] == 'csv'", ns("date_source")),
        fileInput(ns("csv_upload"), "Upload Metadata (.csv)"),
        uiOutput(ns("csv_mapping_ui"))
      ),
      
      h4("3. Parsing Rules", class = "mt-4"),
      accordion(
        id = ns("parse_accordion"),
        accordion_panel(
          "Extraction Method",
          
          conditionalPanel(
            condition = sprintf("input['%s'] == 'csv'", ns("date_source")),
            helpText("Dates will be extracted directly from your CSV. Proceed to Format Transformation if needed.")
          ),
          
          conditionalPanel(
            condition = sprintf("input['%s'] == 'id'", ns("date_source")),
            radioButtons(ns("extract_method"), "Extract date string by:", 
                         choices = c("Delimiter & Order" = "order", "Regular Expression" = "regex")),
            
            conditionalPanel(
              condition = sprintf("input['%s'] == 'order'", ns("extract_method")),
              textInput(ns("delim"), "Delimiter (e.g., _ or |)", value = "_"),
              textInput(ns("order_pos"), "Order (first, last, or index)", value = "last")
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'regex'", ns("extract_method")),
              textInput(ns("regex_pat"), "REGEX Pattern", value = "[0-9]{4}(\\.[0-9]+)?")
            )
          )
        ),
        
        accordion_panel(
          "Format Transformation",
          radioButtons(ns("parse_type"), "Parse as:", 
                       choices = c("Number (Decimal Year)" = "number", "Calendar Date" = "calendar")),
          
          conditionalPanel(
            condition = sprintf("input['%s'] == 'number'", ns("parse_type")),
            checkboxInput(ns("enable_2digit"), "Enable 2-digit year correction", value = FALSE),
            conditionalPanel(
              condition = sprintf("input['%s'] == true", ns("enable_2digit")),
              numericInput(ns("add_val"), "Add to all:", 1900),
              numericInput(ns("threshold"), "...unless less than:", 16),
              numericInput(ns("add_cond"), "...in which case add:", 2000)
            )
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'calendar'", ns("parse_type")),
            textInput(ns("date_fmt"), "Date Format (e.g., yyyy-MM-dd)", value = "yyyy-MM-dd"),
            helpText("Will be automatically converted to decimal years.")
          )
        )
      ),
      
      h4("4. Run Analysis", class = "mt-4"),
      numericInput(ns("nreps"), "Mantel Permutations", value = 1000, min = 100),
      div(class = "d-flex gap-2 mt-2",
          actionButton(ns("run_btn"), "Run Mantel Test", class = "btn-primary flex-grow-1"),
          actionButton(ns("skip_btn"), "Send to RTT (Skip Mantel Test)", class = "btn-primary flex-grow-1")
      ),
      div(class = "d-flex gap-2 mt-2",
          actionButton(ns("reset_mantel_btn"), "Reset", class = "btn-outline-secondary w-100")
      )
    ),
    
    mainPanel(
      width = 8,
      card(
        card_header(
          class = "bg-light d-flex justify-content-between align-items-center",
          strong(" Data Parsing Validation (Double-click cell to edit)"),
          span(class="badge bg-primary", "Auto-updates based on rules")
        ),
        card_body(DTOutput(ns("preview_table")))
      ),
      
      uiOutput(ns("status_ui")),
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header(
            class = "d-flex justify-content-between align-items-center bg-light",
            "Mantel Regression",
            div(
              class = "d-flex gap-2",
              downloadButton(ns("download_mantel_data"), "Export Plot Data", class = "btn btn-sm btn-outline-secondary"),
              downloadButton(ns("download_mantel_pdf"), "Export Vector PDF", class = "btn btn-sm btn-outline-secondary")
            )
          ),
          card_body(plotlyOutput(ns("scatter_plot"), height = "400px"), class = "p-0")
        ),
        uiOutput(ns("tree_ui"))
      )
    )
  )
}

# ---------------------------------------------------------
# Server 
# ---------------------------------------------------------
mod01_Mantel_server <- function(id, global_state) {
  moduleServer(id, function(input, output, session) {
    
    raw_seq_data <- reactiveVal(NULL)
    csv_data <- reactiveVal(NULL)
    parsed_df_rv <- reactiveVal(NULL)
    mantel_res <- reactiveVal(NULL)

    parse_uploaded_file <- function(filepath, fmt) {
      lines <- readLines(filepath, warn = FALSE)
      lines_trim <- trimws(lines)
      lines_non_empty <- lines_trim[nzchar(lines_trim)]
      if (length(lines_non_empty) == 0) {
        stop("The uploaded file is empty.")
      }
      
      first_line <- lines_non_empty[1]
      
      detected_fmt <- fmt
      if (grepl("^#NEXUS", first_line, ignore.case = TRUE)) {
        detected_fmt <- "nexus"
      } else if (startsWith(first_line, ">")) {
        detected_fmt <- "fasta"
      } else if (startsWith(first_line, "(") || any(grepl("\\(.*:.*;", lines_non_empty))) {
        detected_fmt <- "newick"
      }
      
      if (!identical(detected_fmt, fmt)) {
        updateSelectInput(session, "file_format", selected = detected_fmt)
      }
      
      if (detected_fmt %in% c("newick", "nexus")) {
        phy <- NULL
        if (detected_fmt == "nexus") {
          phy <- tryCatch(
            ape::read.nexus(filepath),
            error = function(e) {
              nwk_text <- paste(lines_non_empty, collapse = "")
              ape::read.tree(text = nwk_text)
            }
          )
        } else {
          nwk_text <- paste(lines_non_empty, collapse = "")
          phy <- ape::read.tree(text = nwk_text)
        }
        
        if (inherits(phy, "multiPhylo")) {
          phy <- phy[[1]]
        }
        if (is.null(phy) || !inherits(phy, "phylo")) {
          stop("Failed to parse valid phylogenetic tree from file (check Newick/NEXUS syntax).")
        }
        
        phy$tip.label <- gsub("^['\"]|['\"]$", "", trimws(phy$tip.label))
        
        if (is.null(phy$edge.length)) {
          stop("The uploaded tree does not contain branch lengths.")
        }
        
        g_mat <- ape::cophenetic.phylo(phy)
        g_mat <- g_mat[phy$tip.label, phy$tip.label, drop = FALSE]
        
        return(list(
          names = phy$tip.label,
          phy = phy,
          g = g_mat
        ))
      } else {
        return(read_sequence_data(filepath, detected_fmt, "_"))
      }
    }
    
    observeEvent(input$reset_mantel_btn, {
      shinyjs::reset("file_upload") 
      shinyjs::reset("csv_upload")  
      raw_seq_data(NULL)
      csv_data(NULL)
      parsed_df_rv(NULL)
      mantel_res(NULL)
      showNotification("Mantel Test session reset.", type = "message")
    })
    
    observeEvent(list(input$file_upload, input$file_format), {
      req(input$file_upload)
      tryCatch({
        data <- parse_uploaded_file(input$file_upload$datapath, input$file_format)
        raw_seq_data(data)
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error", duration = 8)
      })
    }, ignoreInit = TRUE)
    
    observeEvent(input$csv_upload, {
      req(input$csv_upload)
      tryCatch({
        df <- read.csv(input$csv_upload$datapath, stringsAsFactors = FALSE, check.names = FALSE)
        csv_data(df)
      }, error = function(e) {
        showNotification("Failed to read CSV. Please ensure it is a valid comma-separated file.", type = "error")
      })
    })
    
    output$csv_mapping_ui <- renderUI({
      req(csv_data())
      cols <- colnames(csv_data())
      ns <- session$ns
      tagList(
        selectInput(ns("csv_taxon_col"), "Select Taxon ID Column", choices = cols),
        selectInput(ns("csv_date_col"), "Select Date Column", choices = cols)
      )
    })
    
    observe({
      req(raw_seq_data())
      labels <- raw_seq_data()$names
      extracted_strs <- rep(NA_character_, length(labels))
      
      if (input$date_source == "id") {
        extracted_strs <- extract_date_strings(
          labels, input$extract_method, input$delim, input$order_pos, input$regex_pat
        )
      } else if (input$date_source == "csv") {
        if (!is.null(csv_data()) && !is.null(input$csv_taxon_col) && !is.null(input$csv_date_col)) {
          df <- csv_data()
          csv_ids <- gsub("^['\"]|['\"]$", "", trimws(as.character(df[[input$csv_taxon_col]])))
          match_idx <- match(labels, csv_ids)
          
          if (all(is.na(match_idx))) {
            match_idx <- sapply(labels, function(lbl) {
              tokens <- trimws(unlist(strsplit(lbl, "[|]")))
              hit <- which(csv_ids %in% c(lbl, tokens))
              if (length(hit) > 0) hit[1] else NA_integer_
            })
          }
          extracted_strs <- as.character(df[[input$csv_date_col]][match_idx])
        }
      }
      
      decimals <- convert_to_decimal_years(
        extracted_strs, input$parse_type, 
        input$add_val, input$threshold, input$add_cond, input$date_fmt,
        enable_2digit = input$enable_2digit
      )
      
      df <- data.frame(
        Taxon_ID = labels,
        Extracted_String = extracted_strs,
        Decimal_Year = round(decimals, 4),
        stringsAsFactors = FALSE
      )
      parsed_df_rv(df)
    })
    
    output$preview_table <- renderDT({
      req(parsed_df_rv())
      datatable(
        parsed_df_rv(),
        colnames = c("Taxon ID", "Extracted String", "Decimal Year"),
        selection = 'none',
        editable = list(target = 'cell', disable = list(columns = c(0, 1))), 
        options = list(pageLength = 8, scrollX = TRUE, dom = 'rtp') 
      ) %>%
        formatStyle(
          'Decimal_Year', 
          backgroundColor = styleEqual(NA, '#ffdddd')
        )
    })
    
    observeEvent(input$preview_table_cell_edit, {
      info <- input$preview_table_cell_edit
      df <- parsed_df_rv()
      df[info$row, info$col + 1] <- as.numeric(info$value)
      parsed_df_rv(df)
    })
    observeEvent(input$skip_btn, {
      req(raw_seq_data(), parsed_df_rv())
      final_dates <- parsed_df_rv()$Decimal_Year
      
      if (any(is.na(final_dates))) {
        showNotification("Error: NA values present in Decimal_Year. Please fix them in the table.", type = "error")
        return()
      }
      global_state$phylo_tree <- raw_seq_data()$phy
      global_state$extracted_dates <- final_dates
      global_state$mantel_p_value <- NULL
      global_state$has_confounding <- NULL
      global_state$clusters <- NULL
      
      showNotification("Data successfully loaded! Please switch directly to the RTT module.", type = "message")
    })
    
    observeEvent(input$run_btn, {
      req(raw_seq_data(), parsed_df_rv())
      final_dates <- parsed_df_rv()$Decimal_Year
      
      if (any(is.na(final_dates))) {
        showNotification("Error: NA values present in Decimal_Year. Please fix them in the table.", type = "error")
        return()
      }
      
      g_mat <- raw_seq_data()$g
      phy_tree <- raw_seq_data()$phy
      n_reps <- input$nreps
      
      id_notif <- showNotification("Running Mantel Test ...", duration = NULL, type = "message")
      
      future::future({
        run_mantel_from_data(phy_tree, g_mat, final_dates, n_reps)
      }, seed = TRUE) %...>% (function(out) {
        cluster_assignment <- NULL
        if (isTRUE(out$is_confounded)) {
          cluster_assignment <- get_clusters(out$phylo_tree, out$dates, do.plot = FALSE)
        }
        
        mantel_res(list(mantel = out, clusters = cluster_assignment))
        
        global_state$phylo_tree <- out$phylo_tree
        global_state$extracted_dates <- out$dates
        global_state$mantel_p_value <- out$p_value
        global_state$has_confounding <- out$is_confounded
        global_state$clusters <- cluster_assignment
        
        if (isTRUE(out$is_confounded)) {
           global_state$metadata <- data.frame(
              Taxon_ID = out$phylo_tree$tip.label,
              Cluster = cluster_assignment,
              stringsAsFactors = FALSE
           )
        }
        
        removeNotification(id_notif)
      }) %...!% (function(err) {
        removeNotification(id_notif)
        showNotification(paste("Error:", err$message), type = "error")
      })
    })
    
    output$status_ui <- renderUI({
      req(mantel_res())
      data <- mantel_res()$mantel
      if(data$is_confounded) {
        card(class = "bg-danger text-white mt-4", card_body(
          div(class = "d-flex justify-content-between align-items-center",
              div(
                h5(icon("triangle-exclamation"), " Significant Confounding Detected"),
                p(sprintf("Mantel r = %.3f, P-value = %.4f", data$r_value, data$p_value), class="mb-0")
              ),
              downloadButton(session$ns("download_cluster_csv"), "Download Cluster Metadata (.csv)", class="btn-light btn-sm text-danger font-weight-bold")
          )
        ))
      } else {
        card(class = "bg-primary text-white mt-4", card_body(
          h5(icon("check-circle"), " No Confounding Detected"),
          p(sprintf("Mantel r = %.3f, P-value = %.4f", data$r_value, data$p_value), class="mb-0")
        ))
      }
    })
    
    output$download_cluster_csv <- downloadHandler(
      filename = function() {
        paste0("EvoTempo_Cluster_Metadata_", format(Sys.Date(), "%Y%m%d"), ".csv")
      },
      content = function(file) {
        req(global_state$metadata)
        write.csv(global_state$metadata, file, row.names = FALSE)
      }
    )
    
    output$download_mantel_data <- downloadHandler(
      filename = function() {
        paste0("EvoTempo_Mantel_Data_", format(Sys.Date(), "%Y%m%d"), ".csv")
      },
      content = function(file) {
        req(mantel_res())
        write.csv(mantel_res()$mantel$plot_data, file, row.names = FALSE)
      }
    )
    
    output$download_mantel_pdf <- downloadHandler(
      filename = function() {
        paste0("EvoTempo_Mantel_Regression_", format(Sys.Date(), "%Y%m%d"), ".pdf")
      },
      content = function(file) {
        req(mantel_res())
        p <- ggplot(mantel_res()$mantel$plot_data, aes(x = Genetic_Distance, y = Temporal_Distance)) +
          geom_point(alpha = 0.3) + 
          geom_smooth(method = "lm", se = FALSE, color = "red") + 
          theme_minimal() +
          labs(x = "Genetic Distance", y = "Temporal Distance")
        
        ggsave(filename = file, plot = p, device = "pdf", width = 6, height = 5)
      }
    )
    
    output$scatter_plot <- renderPlotly({
      req(mantel_res())
      
      p <- ggplot(mantel_res()$mantel$plot_data, aes(x = Genetic_Distance, y = Temporal_Distance)) +
        geom_point(alpha = 0.3) + 
        geom_smooth(method = "lm", se = FALSE, color = "red") + 
        theme_minimal() +
        labs(x = "Genetic Distance", y = "Temporal Distance")
      
      ply <- ggplotly(p, tooltip = "none") %>% layout(hovermode = FALSE)
      ply <- plotly::plotly_build(ply)
      
      if (!is.null(ply$x$data)) {
        ply$x$data <- lapply(ply$x$data, function(trace) {
          trace$hoveron <- NULL
          return(trace)
        })
      }
      
      toWebGL(ply)
    })
    
    output$tree_ui <- renderUI({
      req(mantel_res())
      if (mantel_res()$mantel$is_confounded && !is.null(mantel_res()$clusters)) {
        card(
          card_header(
            class = "d-flex justify-content-between align-items-center bg-light",
            "Phylo-temporal clustering",
            downloadButton(session$ns("download_cluster_pdf"), "Export Vector PDF", class = "btn btn-sm btn-outline-secondary")
          ),
          card_body(plotOutput(session$ns("cluster_tree"), height = "400px"), class = "p-0")
        )
      }
    })
    
    output$download_cluster_pdf <- downloadHandler(
      filename = function() {
        paste0("EvoTempo_Cluster_Tree_", format(Sys.Date(), "%Y%m%d"), ".pdf")
      },
      content = function(file) {
        req(mantel_res()$clusters)
        phy <- mantel_res()$mantel$phylo_tree
        cl <- mantel_res()$clusters
        set.seed(42) 
        cols <- rainbow(max(cl))[sample(max(cl))]
        phy$tip.label <- rep('.', length(phy$tip.label))
        
        pdf(file, width = 8, height = 8)
        plot(phy, tip.color = cols[cl], cex = 4, root.edge = TRUE, no.margin = TRUE)
        dev.off()
      }
    )
    
    output$cluster_tree <- renderPlot({
      req(mantel_res()$clusters)
      phy <- mantel_res()$mantel$phylo_tree
      cl <- mantel_res()$clusters
      set.seed(42) 
      cols <- rainbow(max(cl))[sample(max(cl))]
      phy$tip.label <- rep('.', length(phy$tip.label))
      
      plot(phy, tip.color = cols[cl], cex = 4, root.edge = TRUE, no.margin = TRUE)
    })
  })
}
