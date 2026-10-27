## R/mod02_RTT.R
## EvoTempo: Module 02 - Interactive Root-to-Tip & Regression Grouping

library(shiny)
library(bslib)
library(ggplot2)
library(plotly)
library(ape)
library(ggtree)
library(crosstalk)
library(dplyr)
library(stringr)

# ---------------------------------------------------------
# UI 
# ---------------------------------------------------------
mod02_RTT_ui <- function(id) {
  ns <- NS(id)

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("RTT Regression Settings"),

      selectInput(ns("rooting_method"), "Rooting Strategy",
                  choices = c("Minimize Residual Variance (nve-rms)" = "nve-rms",
                              "Maximize Correlation (r)" = "r")),
      selectInput(ns("grouping_mode"), "Regression Grouping Strategy",
                  choices = c("Global Regression (Strict Clock)" = "global",
                              "Interactive Clades (Local Clocks)" = "local",
                              "Discrete Trait Partitioning" = "metadata")),

      conditionalPanel(
        condition = sprintf("input['%s'] == 'local'", ns("grouping_mode")),
        card(
          class = "bg-primary text-white mb-3",
          card_body(
            class = "p-2",
            h6(" Interactive Clocks"),
            p("Click on internal nodes in the tree to assign separate RTT regression lines to selected clades.",
              style = "font-size:0.9em; margin-bottom: 5px;"),
            actionButton(ns("clear_clocks_btn"), "Clear Local Clocks",
                         class = "btn-light btn-sm w-100")
          )
        )
      ),

      conditionalPanel(
        condition = sprintf("input['%s'] == 'metadata'", ns("grouping_mode")),
        card(
          class = "bg-primary text-white mb-3",
          card_body(
            class = "p-2",
            h6("Upload Trait Data"),
            p("Group isolates based on geographic or host data to compare evolutionary rates.",
              style = "font-size:0.9em; margin-bottom: 5px;"),
            fileInput(ns("trait_csv"), "Metadata (.csv)", accept = ".csv"),
            uiOutput(ns("trait_selector_ui"))
          )
        )
      ),

      checkboxInput(ns("run_permutation"),
                    "Run Permutation Test (calculate P-value)",
                    value = FALSE),

      div(
        class = "d-flex gap-2 mb-3",
        actionButton(ns("calc_rtt_btn"), "Compute RTT",
                     class = "btn-primary flex-grow-1"),
        actionButton(ns("reset_btn"), "Reset",
                     class = "btn-outline-secondary")
      )
    ),

    mainPanel(
      width = 9,
      uiOutput(ns("dependency_warning")),

      uiOutput(ns("rtt_stats_ui")),

      layout_columns(
        col_widths = c(6, 6),

        card(
          card_header(
            class = "d-flex justify-content-between align-items-center bg-light",
            style = "height: 48px;",
            "Phylogenetic Tree",
            downloadButton(ns("download_tree_pdf"), "Export Vector PDF",
                           class = "btn btn-sm btn-outline-secondary")
          ),
          plotlyOutput(ns("tree_plot"), height = "500px")
        ),

        card(
          card_header(
            class = "d-flex justify-content-between align-items-center bg-light",
            style = "height: 48px;",
            "Root-to-tip Regression",
            div(
              class = "d-flex gap-2",
              downloadButton(ns("download_data"), "Export Plot Data",
                             class = "btn btn-sm btn-outline-secondary"),
              downloadButton(ns("download_pdf"), "Export Vector PDF",
                             class = "btn btn-sm btn-outline-secondary")
            )
          ),
          plotlyOutput(ns("rtt_scatter"), height = "500px")
        )
      )
    )
  )
}


# ---------------------------------------------------------
# Server 
# ---------------------------------------------------------
mod02_RTT_server <- function(id, global_state) {
  moduleServer(id, function(input, output, session) {

    rtt_data <- reactiveVal(NULL)
    group_mapping <- reactiveVal(NULL)
    local_clock_counter <- reactiveVal(1)
    selected_clock_nodes <- reactiveVal(setNames(character(0), character(0)))
    raw_trait_data <- reactiveVal(NULL)

    get_group_palette <- function(groups) {
      lvls <- sort(unique(as.character(groups)))
      n <- length(lvls)
      if (n <= 1) {
        return(setNames("#2E6A8C", lvls))
      }
      hues <- seq(15, 375, length = n + 1)
      cols <- grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
      setNames(cols, lvls)
    }

    compute_rtt_core <- function() {
      req(global_state$phylo_tree, global_state$extracted_dates)

      id_notif <- showNotification(
        "Optimizing root and computing distances...",
        duration = NULL
      )

      tryCatch({

        rtt_res <- run_rtt_regression(
          phy = global_state$phylo_tree,
          dates = global_state$extracted_dates,
          stat = input$rooting_method,
          reroot = TRUE
        )

        p_val <- NA

        if (input$run_permutation) {

          showNotification(
            "Running 1000 permutations...",
            type = "message",
            id = id_notif
          )

          nreps <- 1000
          obs_r <- rtt_res$stats$r_value
          nt <- length(global_state$extracted_dates)

          distances <- dist.nodes(
            rtt_res$phylo_rooted
          )[nt + 1, 1:nt]

          dates_vec <- global_state$extracted_dates
          null_r <- numeric(nreps)

          for (i in 1:nreps) {
            shuff_dates <- sample(dates_vec)
            null_lm <- lm(distances ~ shuff_dates)

            null_r[i] <- sqrt(summary(null_lm)$r.squared) *
              sign(coef(null_lm)[2])
          }

          p_val <- sum(abs(null_r) >= abs(obs_r)) / nreps
        }

        rtt_res$stats$p_value <- p_val

        if (is.null(group_mapping())) {
          tips <- rtt_res$phylo_rooted$tip.label

          group_mapping(
            setNames(rep("Global", length(tips)), tips)
          )

          local_clock_counter(1)
        }

        rtt_data(list(
          raw_df = rtt_res$df_plot,
          phylo_rooted = rtt_res$phylo_rooted,
          stats = rtt_res$stats
        ))

      }, error = function(e) {

        showNotification(
          paste("RTT Error:", e$message),
          type = "error"
        )

      }, finally = {

        removeNotification(id_notif)

      })
    }

    output$dependency_warning <- renderUI({
      if (is.null(global_state$phylo_tree)) {
        card(
          class = "bg-warning text-white mb-3",
          card_body(
            h5(" Waiting for Data"),
            "Please upload data in Tab 1 first."
          )
        )
      }
    })

    observeEvent(input$reset_btn, {
      rtt_data(NULL)
      group_mapping(NULL)
      local_clock_counter(1)
      selected_clock_nodes(setNames(character(0), character(0)))
      raw_trait_data(NULL)
    })

    observeEvent(input$clear_clocks_btn, {
      req(global_state$phylo_tree)

      tips <- global_state$phylo_tree$tip.label
      group_mapping(setNames(rep("Global", length(tips)), tips))
      local_clock_counter(1)
      selected_clock_nodes(setNames(character(0), character(0)))
    })

    observeEvent(input$grouping_mode, {
      if (input$grouping_mode == "global" &&
          !is.null(global_state$phylo_tree)) {

        tips <- global_state$phylo_tree$tip.label
        group_mapping(setNames(rep("Global", length(tips)), tips))
        local_clock_counter(1)
        selected_clock_nodes(setNames(character(0), character(0)))
      }
    })


    observeEvent(input$trait_csv, {
      req(input$trait_csv)

      tryCatch({

        file_path <- input$trait_csv$datapath

        first_line <- readLines(file_path, n = 1, warn = FALSE)

        if (length(first_line) == 0 || !nzchar(first_line)) {
          stop("The uploaded metadata file is empty.")
        }

        delimiter_counts <- c(
          comma = lengths(
            gregexpr(",", first_line, fixed = TRUE)
          ),
          semicolon = lengths(
            gregexpr(";", first_line, fixed = TRUE)
          ),
          tab = lengths(
            gregexpr("\t", first_line, fixed = TRUE)
          )
        )

        delimiter <- ","

        if (delimiter_counts["semicolon"] >
            delimiter_counts["comma"] &&
            delimiter_counts["semicolon"] >= delimiter_counts["tab"]) {

          delimiter <- ";"

        } else if (delimiter_counts["tab"] >
                   delimiter_counts["comma"] &&
                   delimiter_counts["tab"] >
                   delimiter_counts["semicolon"]) {

          delimiter <- "\t"
        }

        df <- utils::read.table(
          file_path,
          header = TRUE,
          sep = delimiter,
          stringsAsFactors = FALSE,
          check.names = FALSE,
          comment.char = "",
          quote = "\"",
          strip.white = TRUE
        )

        if (ncol(df) == 1) {

          alternative_delimiters <- c(",", ";", "\t")

          for (candidate_sep in alternative_delimiters) {

            if (identical(candidate_sep, delimiter)) {
              next
            }

            candidate_df <- tryCatch(
              utils::read.table(
                file_path,
                header = TRUE,
                sep = candidate_sep,
                stringsAsFactors = FALSE,
                check.names = FALSE,
                comment.char = "",
                quote = "\"",
                strip.white = TRUE
              ),
              error = function(e) NULL
            )

            if (!is.null(candidate_df) && ncol(candidate_df) > 1) {
              df <- candidate_df
              break
            }
          }
        }

        if (length(names(df)) > 0) {
          names(df)[1] <- sub(
            "^\ufeff",
            "",
            names(df)[1]
          )
        }

        names(df) <- trimws(names(df))

        if (ncol(df) < 2) {
          stop(
            paste(
              "Metadata must contain at least two separate columns.",
              "Please check whether the file uses commas, semicolons,",
              "or tabs as the delimiter."
            )
          )
        }

        raw_trait_data(df)

        showNotification(
          paste(
            "Metadata loaded successfully:",
            nrow(df), "rows and", ncol(df), "columns."
          ),
          type = "message"
        )

      }, error = function(e) {

        showNotification(
          paste("Failed to parse Trait Metadata:", e$message),
          type = "error",
          duration = 8
        )

      })
    })


    output$trait_selector_ui <- renderUI({
      req(raw_trait_data())

      cols <- colnames(raw_trait_data())
      ns <- session$ns

      tagList(
        selectInput(
          ns("trait_taxon_col"),
          "Taxon ID Column",
          choices = cols,
          selected = cols[1]
        ),

        selectInput(
          ns("trait_group_col"),
          "Grouping Trait (e.g., Geography)",
          choices = cols,
          selected = if (length(cols) >= 2) cols[2] else cols[1]
        ),

        actionButton(
          ns("apply_trait_btn"),
          "Apply Trait Grouping",
          class = "btn-light btn-sm w-100"
        )
      )
    })


    observeEvent(input$apply_trait_btn, {
      req(
        raw_trait_data(),
        input$trait_taxon_col,
        input$trait_group_col,
        global_state$phylo_tree
      )

      df <- raw_trait_data()

      if (input$trait_taxon_col == input$trait_group_col) {
        showNotification(
          "Taxon ID Column and Grouping Trait must be different columns.",
          type = "error"
        )
        return()
      }

      if (is.null(rtt_data()) && !is.null(global_state$extracted_dates)) {
        compute_rtt_core()
      }

      phy <- if (!is.null(rtt_data())) {
        rtt_data()$phylo_rooted
      } else {
        global_state$phylo_tree
      }

      tips <- as.character(phy$tip.label)

      metadata_ids <- trimws(as.character(df[[input$trait_taxon_col]]))
      tree_ids <- trimws(tips)

      match_idx <- match(tree_ids, metadata_ids)

      if (all(is.na(match_idx))) {
        non_empty_meta <- which(nzchar(metadata_ids))
        match_idx <- sapply(tree_ids, function(tid) {
          hits <- non_empty_meta[
            startsWith(tid, paste0(metadata_ids[non_empty_meta], "_")) |
              startsWith(tid, paste0(metadata_ids[non_empty_meta], "|")) |
              startsWith(tid, paste0(metadata_ids[non_empty_meta], "/")) |
              startsWith(tid, metadata_ids[non_empty_meta])
          ]
          if (length(hits) > 0) {
            hits[which.max(nchar(metadata_ids[hits]))]
          } else {
            NA_integer_
          }
        })
      }

      n_matched <- sum(!is.na(match_idx))
      if (n_matched == 0) {
        showNotification(
          "No Taxon IDs matched between tree tips and the selected column. Please check 'Taxon ID Column'.",
          type = "error",
          duration = 8
        )
        return()
      }

      mapped_groups <- as.character(
        df[[input$trait_group_col]][match_idx]
      )

      mapped_groups[
        is.na(mapped_groups) | trimws(mapped_groups) == ""
      ] <- "Unknown"

      selected_clock_nodes(setNames(character(0), character(0)))
      group_mapping(setNames(mapped_groups, tips))

      showNotification(
        sprintf(
          "Trait mapping applied (%d/%d taxa matched across %d groups).",
          n_matched, length(tips), length(unique(mapped_groups))
        ),
        type = "message"
      )
    })


    click_observer <- observeEvent(
      plotly::event_data("plotly_click", source = "tree_source"),
      {
        req(input$grouping_mode == "local", global_state$phylo_tree)

        click_data <- plotly::event_data(
          "plotly_click",
          source = "tree_source"
        )
        req(!is.null(click_data))

        node_id <- click_data$customdata

        if (is.null(node_id) && !is.null(click_data$text)) {
          extracted <- stringr::str_extract(
            click_data$text,
            "(?<=Node: )[0-9]+"
          )
          node_id <- as.numeric(extracted)
        }

        if (is.null(node_id) || length(node_id) == 0 || is.na(node_id[1])) {
          return()
        }

        node_id <- as.numeric(node_id[1])

        phy <- if (!is.null(rtt_data())) {
          rtt_data()$phylo_rooted
        } else {
          global_state$phylo_tree
        }

        nt <- length(phy$tip.label)

        if (node_id > nt) {
          clade_tips <- ape::extract.clade(
            phy,
            node = node_id
          )$tip.label

          current_mapping <- group_mapping()

          if (is.null(current_mapping)) {
            current_mapping <- setNames(
              rep("Global", nt),
              phy$tip.label
            )
          }

          new_clock_name <- paste(
            "Local Clock",
            local_clock_counter() + 1
          )

          current_mapping[clade_tips] <- new_clock_name
          group_mapping(current_mapping)
          local_clock_counter(local_clock_counter() + 1)
          cur_nodes <- selected_clock_nodes()
          cur_nodes[as.character(node_id)] <- new_clock_name
          selected_clock_nodes(cur_nodes)

          showNotification(
            paste(
              new_clock_name,
              "assigned to",
              length(clade_tips),
              "sequences. Click 'Compute RTT' to update analysis."
            ),
            type = "message"
          )
        }
      },
      ignoreNULL = TRUE,
      suspended = TRUE
    )


    observeEvent(input$calc_rtt_btn, {
      compute_rtt_core()
    })


    combined_plot_data <- reactive({
      req(rtt_data(), group_mapping())

      df <- rtt_data()$raw_df
      gm <- group_mapping()
      tax_chr <- as.character(df$Taxon)
      grp <- unname(gm[tax_chr])

      if (any(is.na(grp))) {
        gm_trim <- setNames(as.character(gm), trimws(names(gm)))
        grp[is.na(grp)] <- unname(gm_trim[trimws(tax_chr[is.na(grp)])])
      }
      grp[is.na(grp) | !nzchar(trimws(grp))] <- "Unknown"

      df$Group <- as.character(grp)

      df
    })


    group_stats <- reactive({
      req(combined_plot_data())

      df <- combined_plot_data()

      stats_list <- lapply(unique(df$Group), function(grp) {

        sub_df <- df[df$Group == grp & !is.na(df$Date) & !is.na(df$Distance), , drop = FALSE]

        if (nrow(sub_df) > 2 && length(unique(sub_df$Date)) > 1) {

          fit <- tryCatch(lm(Distance ~ Date, data = sub_df), error = function(e) NULL)
          slope <- if (!is.null(fit)) unname(coef(fit)[2]) else NA_real_

          if (!is.null(fit) && !is.na(slope) && is.finite(slope) && abs(slope) > 1e-15) {
            tmrca <- unname(-coef(fit)[1] / slope)
            r_sq <- summary(fit)$r.squared

            label_str <- sprintf(
              "%s: Rate = %.2e | tMRCA = %.0f | R² = %.3f",
              grp, slope, tmrca, r_sq
            )

            return(data.frame(
              Group = grp,
              Rate = slope,
              tMRCA = tmrca,
              R2 = r_sq,
              Label = label_str,
              stringsAsFactors = FALSE
            ))
          }
        }

        data.frame(
          Group = grp,
          Rate = NA_real_,
          tMRCA = NA_real_,
          R2 = NA_real_,
          Label = NA_character_,
          stringsAsFactors = FALSE
        )
      })

      do.call(rbind, stats_list)
    })


    output$rtt_stats_ui <- renderUI({
      req(rtt_data())

      if (!is.null(group_stats()) && nrow(group_stats()) > 1) {

        stats_df <- group_stats()

        cards <- lapply(1:nrow(stats_df), function(i) {

          rate_val <- as.numeric(stats_df$Rate[i])
          tmrca_val <- as.numeric(stats_df$tMRCA[i])
          r2_val <- as.numeric(stats_df$R2[i])

          card(
            class = "border-primary flex-fill",
            card_body(
              class = "p-2",
              h6(
                stats_df$Group[i],
                class = "text-primary m-0 mb-1"
              ),
              div(
                class = "d-flex justify-content-between gap-3",
                style = "font-size:0.95em;",
                span(
                  strong("Rate: "),
                  if (is.na(rate_val)) "N/A" else sprintf("%.4e", rate_val)
                ),
                span(
                  strong("tMRCA: "),
                  if (is.na(tmrca_val)) "N/A" else sprintf("%.2f", tmrca_val)
                ),
                span(
                  strong("R-squared: "),
                  if (is.na(r2_val)) "N/A" else sprintf("%.4f", r2_val)
                )
              )
            )
          )
        })

        div(
          class = "d-flex flex-wrap gap-2 mb-3",
          cards
        )

      } else {

        s <- rtt_data()$stats

        p_val_ui <- if (!is.na(s$p_value)) {
          span(
            strong("p-value: "),
            sprintf("%.4f", s$p_value),
            class = "text-danger"
          )
        } else {
          NULL
        }

        card(
          class = "bg-light border-0 shadow-sm mb-3",
          card_body(
            class = "p-3 d-flex justify-content-around align-items-center flex-wrap",
            span(
              strong("Clock Rate: "),
              sprintf("%.4e", s$clock_rate)
            ),
            span(
              strong("tMRCA: "),
              sprintf("%.2f", s$tmrca)
            ),
            span(
              strong("Correlation (r): "),
              sprintf("%.4f", s$r_value)
            ),
            span(
              strong("R-squared: "),
              sprintf("%.4f", s$r_squared)
            ),
            p_val_ui
          )
        )
      }
    })

    build_scatter_plot <- function(
      df,
      s,
      is_web = FALSE,
      annot_labels = NULL
    ) {

      b_size <- ifelse(is_web, 11, 14)
      t_size <- ifelse(is_web, 12, 14)
      ax_size <- ifelse(is_web, 10, 12)
      lg_title <- ifelse(is_web, 10, 12)
      lg_text <- ifelse(is_web, 9, 11)
      text_size <- ifelse(is_web, 3.5, 4.5)
      valid_lm_groups <- unique(df$Group)[sapply(unique(df$Group), function(g) {
        sub_d <- df$Date[df$Group == g & !is.na(df$Date) & !is.na(df$Distance)]
        length(sub_d) >= 2 && length(unique(sub_d)) >= 2
      })]
      df_smooth <- df[df$Group %in% valid_lm_groups, , drop = FALSE]

      p <- suppressWarnings({

        ggplot(
          df,
          aes(x = Date, y = Distance, color = Group)
        ) +
          geom_point(
            aes(text = Taxon),
            alpha = 0.7,
            size = 2
          ) +
          geom_smooth(
            data = df_smooth,
            method = "lm",
            formula = y ~ x,
            se = FALSE,
            linetype = "solid",
            linewidth = 1
          ) +
          theme_bw(base_size = b_size) +
          theme(
            panel.border = element_rect(
              color = "black",
              fill = NA,
              linewidth = 1
            ),
            panel.grid.major = element_line(
              color = "grey80",
              linetype = "dashed"
            ),
            panel.grid.minor = element_line(
              color = "grey92",
              linetype = "dashed"
            ),
            axis.title = element_text(
              size = t_size,
              face = "bold",
              color = "#2C3E50"
            ),
            axis.text = element_text(
              size = ax_size,
              color = "black"
            ),
            legend.title = element_text(
              size = lg_title,
              face = "bold"
            ),
            legend.text = element_text(
              size = lg_text
            ),
            legend.position = "top",
            legend.justification = "center",
            legend.direction = "horizontal",
            legend.box = "horizontal"
          ) +
          guides(
            color = guide_legend(
              nrow = 1,
              byrow = TRUE,
              override.aes = list(
                alpha = 1,
                linewidth = 1.2,
                size = 2.5
              )
            )
          ) +
          labs(
            x = "Sampling Date",
            y = "Root-to-tip Distance",
            color = "Regression Group"
          )
      })

      if (length(unique(df$Group)) == 1) {
        p <- p + scale_color_manual(
          values = setNames(
            "#2E6A8C",
            unique(as.character(df$Group))
          )
        )

        if (length(p$layers) >= 2) {
          p$layers[[2]]$aes_params$colour <- "#C0392B"
        }

        stats_text <- sprintf(
          "Rate: %.4e\ntMRCA: %.2f\nr: %.4f",
          s$clock_rate,
          s$tmrca,
          s$r_value
        )

        if (!is.na(s$p_value)) {
          stats_text <- paste0(
            stats_text,
            sprintf("\np-val: %.4f", s$p_value)
          )
        }

        if (!is_web) {
          p <- p +
            annotate(
              "label",
              x = Inf,
              y = Inf,
              label = stats_text,
              hjust = 1.05,
              vjust = 1.1,
              size = text_size,
              fontface = "bold",
              color = "#2C3E50",
              fill = "white",
              alpha = 0.9,
              label.size = NA
            )
        }


      } else {
        p <- p + scale_color_manual(values = get_group_palette(df$Group))

        if (!is_web && !is.null(annot_labels)) {

          valid_labels <- na.omit(annot_labels)

          if (length(valid_labels) > 0) {

            combo_text <- paste(
              valid_labels,
              collapse = "\n"
            )

            p <- p +
              annotate(
                "label",
                x = Inf,
                y = Inf,
                label = combo_text,
                hjust = 1.05,
                vjust = 1.1,
                size = text_size,
                fontface = "bold",
                fill = "white",
                alpha = 0.9,
                label.size = NA
              )

          }
        }
      }

      return(p)
    }


    static_scatter_pdf <- reactive({
      req(combined_plot_data(), rtt_data())

      build_scatter_plot(
        combined_plot_data(),
        rtt_data()$stats,
        is_web = FALSE,
        annot_labels = group_stats()$Label
      )
    })


    static_scatter_web <- reactive({
      req(combined_plot_data(), rtt_data())

      build_scatter_plot(
        combined_plot_data(),
        rtt_data()$stats,
        is_web = TRUE,
        annot_labels = group_stats()$Label
      )
    })


    output$rtt_scatter <- renderPlotly({

      suppressWarnings(suppressMessages({

        ggplotly(
          static_scatter_web(),
          tooltip = c("text", "x", "y")
        ) %>%
          layout(
            hovermode = "closest",
            legend = list(
              orientation = "h",
              x = 0.5,
              xanchor = "center",
              y = 1.12,
              yanchor = "bottom",
              title = list(text = "Regression Group")
            ),
            margin = list(t = 75)
          ) %>%
          config(displayModeBar = TRUE)

      }))
    })


    output$download_data <- downloadHandler(
      filename = function() {
        paste0(
          "EvoTempo_RTT_Grouped_",
          format(Sys.Date(), "%Y%m%d"),
          ".csv"
        )
      },
      content = function(file) {
        req(combined_plot_data())
        write.csv(
          combined_plot_data(),
          file,
          row.names = FALSE
        )
      }
    )


    output$download_pdf <- downloadHandler(
      filename = function() {
        paste0(
          "EvoTempo_RTT_",
          format(Sys.Date(), "%Y%m%d"),
          ".pdf"
        )
      },
      content = function(file) {
        req(static_scatter_pdf())

        ggsave(
          filename = file,
          plot = static_scatter_pdf(),
          device = "pdf",
          width = 8,
          height = 6
        )
      }
    )

    build_tree_plot <- reactive({
      req(global_state$phylo_tree)

      phy <- if (!is.null(rtt_data())) {
        rtt_data()$phylo_rooted
      } else {
        global_state$phylo_tree
      }

      options(ignore.negative.edge = TRUE)

      p_tree <- ggtree(phy)

      df <- p_tree$data

      cur_mapping <- group_mapping()
      has_local_clocks <- (
        input$grouping_mode == "local" &&
          !is.null(cur_mapping) &&
          any(cur_mapping != "Global")
      )
      has_trait_groups <- (
        input$grouping_mode == "metadata" &&
          !is.null(cur_mapping) &&
          !(length(unique(cur_mapping)) == 1 && unique(cur_mapping)[1] == "Global")
      )

      node_group <- rep("Global", nrow(df))
      names(node_group) <- as.character(df$node)

      if (has_local_clocks || has_trait_groups) {
        pal <- get_group_palette(cur_mapping)

        for (i in seq_len(nrow(df))) {
          nd <- df$node[i]
          if (df$isTip[i]) {
            lbl <- as.character(df$label[i])
            if (!is.na(lbl) && lbl %in% names(cur_mapping)) {
              node_group[as.character(nd)] <- cur_mapping[[lbl]]
            }
          } else if (has_local_clocks) {
            clade_tips <- tryCatch(
              ape::extract.clade(phy, node = nd)$tip.label,
              error = function(e) character(0)
            )
            if (length(clade_tips) > 0) {
              tip_grps <- unique(cur_mapping[clade_tips])
              if (length(tip_grps) == 1 && !is.na(tip_grps) && tip_grps != "Global") {
                node_group[as.character(nd)] <- tip_grps
              }
            }
          }
        }
      }

      df$clade_group <- node_group[as.character(df$node)]

      df_segs <- df %>%
        filter(node != parent) %>%
        left_join(
          df %>% select(node, px = x, py = y),
          by = c("parent" = "node")
        )
      seg_colors <- rep("#5D6D7E", nrow(df_segs))
      if (has_local_clocks) {
        is_local_seg <- df_segs$clade_group != "Global"
        seg_colors[is_local_seg] <- pal[df_segs$clade_group[is_local_seg]]
      }

      df_tips <- filter(df, isTip == TRUE)
      tip_colors <- rep("#2C3E50", nrow(df_tips))
      if (has_local_clocks) {
        is_local_tip <- df_tips$clade_group != "Global"
        tip_colors[is_local_tip] <- pal[df_tips$clade_group[is_local_tip]]
      } else if (has_trait_groups) {
        tip_colors <- unname(pal[df_tips$clade_group])
        tip_colors[is.na(tip_colors)] <- "#2C3E50"
      }

      p_tree <- p_tree +
        geom_segment(
          data = df_segs,
          aes(
            x = px,
            y = py,
            xend = px,
            yend = y
          ),
          color = seg_colors,
          linewidth = 0.5
        ) +
        geom_segment(
          data = df_segs,
          aes(
            x = px,
            y = y,
            xend = x,
            yend = y
          ),
          color = seg_colors,
          linewidth = 0.5
        ) +
        geom_point(
          data = df_tips,
          aes(
            x = x,
            y = y,
            text = label
          ),
          color = tip_colors,
          size = 1.5
        )

      if (input$grouping_mode == "local") {

        df_internal <- filter(df, isTip == FALSE)
        int_colors <- rep("darkorange", nrow(df_internal))
        int_sizes <- rep(3, nrow(df_internal))
        int_texts <- paste("Node:", df_internal$node)

        if (has_local_clocks) {
          is_local_int <- df_internal$clade_group != "Global"
          int_colors[is_local_int] <- pal[df_internal$clade_group[is_local_int]]
          int_texts[is_local_int] <- paste0(
            "Node: ", df_internal$node[is_local_int],
            " (", df_internal$clade_group[is_local_int], ")"
          )

          # Slightly enlarge the clicked clade root node(s)
          sel_nodes <- names(selected_clock_nodes())
          is_clicked_root <- as.character(df_internal$node) %in% sel_nodes & is_local_int
          int_sizes[is_clicked_root] <- 4.2
        }

        p_tree <- p_tree +
          geom_point(
            data = df_internal,
            aes(
              x = x,
              y = y,
              customdata = node,
              text = int_texts
            ),
            color = int_colors,
            alpha = 0.85,
            size = int_sizes
          )
      }

      p_tree <- p_tree + theme_tree2()

      p_tree
    })


    output$tree_plot <- renderPlotly({
      req(global_state$phylo_tree)

      suppressWarnings(suppressMessages({

        p_tree <- build_tree_plot()

        p <- plotly::ggplotly(
          p_tree,
          tooltip = "text",
          source = "tree_source"
        ) %>%
          plotly::layout(dragmode = "pan") %>%
          plotly::config(displayModeBar = TRUE) %>%
          plotly::event_register("plotly_click")
        session$onFlushed(
          function() {
            click_observer$resume()
          },
          once = TRUE
        )

        p

      }))
    })

    output$download_tree_pdf <- downloadHandler(
      filename = function() {
        paste0(
          "EvoTempo_PhyloTree_",
          format(Sys.Date(), "%Y%m%d"),
          ".pdf"
        )
      },
      content = function(file) {
        req(build_tree_plot())

        ggsave(
          filename = file,
          plot = build_tree_plot(),
          device = "pdf",
          width = 8,
          height = 8
        )
      }
    )

  })
}
