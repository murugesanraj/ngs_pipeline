suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(DT)
  library(plotly)
  library(jsonlite)
  library(yaml)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(htmltools)
})

app_dir <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(app_dir, "R", "utils.R"))) {
  stop("Start the app from its app/ directory or use app/run_app.R")
}
source(file.path(app_dir, "R", "utils.R"), local = TRUE)
cfg <- read_app_config(app_dir)
portal_root <- cfg$paths$portal_root
max_rows <- as.integer(cfg$display$max_table_rows %||% 50000)
max_download_bytes <- as.numeric(cfg$app$max_download_mb %||% 250) * 1024^2
refresh_ms <- as.numeric(cfg$app$refresh_seconds %||% 60) * 1000

theme <- bs_theme(
  version = 5,
  bg = "#FFFFFF",
  fg = cfg$display$primary_color %||% "#111111",
  primary = cfg$display$accent_color %||% "#F1B82D"
)

ui <- page_sidebar(
  title = cfg$app$title %||% "NGS Results Portal",
  theme = theme,
  fillable = TRUE,
  sidebar = sidebar(
    width = 310,
    p(class = "text-muted", paste("Signed in as", current_hpc_user())),
    selectInput("project_id", "Project", choices = character()),
    actionButton("refresh", "Refresh projects", icon = icon("rotate"), class = "btn-primary"),
    hr(),
    uiOutput("project_details"),
    hr(),
    p(class = "small text-muted", "Read-only view. Access follows Hellbender filesystem permissions.")
  ),
  navset_card_tab(
    nav_panel(
      "Overview",
      layout_columns(
        value_box(title = "Samples", value = textOutput("sample_count"), showcase = icon("vials")),
        value_box(title = "Reads after QC", value = textOutput("read_count"), showcase = icon("dna")),
        value_box(title = "Mean Q30", value = textOutput("mean_q30"), showcase = icon("circle-check")),
        value_box(title = "Mean mapping", value = textOutput("mean_mapping"), showcase = icon("chart-column")),
        col_widths = c(3, 3, 3, 3)
      ),
      layout_columns(
        card(card_header("Project summary"), uiOutput("project_summary")),
        card(card_header("PCA / primary figure"), imageOutput("primary_figure", height = "420px")),
        col_widths = c(5, 7)
      )
    ),
    nav_panel(
      "QC",
      layout_columns(
        card(card_header("Post-QC metrics"), plotlyOutput("qc_plot", height = "420px")),
        card(card_header("Sample QC table"), DTOutput("qc_table")),
        col_widths = c(6, 6)
      )
    ),
    nav_panel(
      "Mapping / Assembly",
      card(card_header("Mapping and coverage"), plotlyOutput("mapping_plot", height = "420px")),
      card(card_header("Metrics"), DTOutput("mapping_table"))
    ),
    nav_panel(
      "Expression",
      layout_columns(
        card(
          card_header("Contrast"),
          selectInput("contrast_id", "Comparison", choices = character()),
          imageOutput("volcano", height = "440px")
        ),
        card(card_header("Differential-expression preview"), DTOutput("de_table")),
        col_widths = c(5, 7)
      )
    ),
    nav_panel(
      "Files",
      card(
        card_header("Approved result downloads"),
        selectInput("download_id", "File", choices = character()),
        downloadButton("download_result", "Download selected file", class = "btn-primary"),
        br(), br(),
        DTOutput("download_table")
      )
    ),
    nav_panel(
      "Provenance",
      card(card_header("Run provenance"), verbatimTextOutput("provenance"))
    )
  )
)

server <- function(input, output, session) {
  project_catalog <- reactiveVal(discover_projects(portal_root))

  refresh_projects <- function() {
    catalog <- discover_projects(portal_root)
    project_catalog(catalog)
    choices <- if (nrow(catalog)) setNames(catalog$id, catalog$name) else character()
    current <- input$project_id %||% ""
    selected <- if (nzchar(current) && current %in% catalog$id) {
      current
    } else if (nrow(catalog)) {
      catalog$id[1]
    } else {
      NULL
    }
    updateSelectInput(session, "project_id", choices = choices, selected = selected)
  }

  observeEvent(input$refresh, refresh_projects(), ignoreInit = TRUE)
  observe({
    invalidateLater(refresh_ms, session)
    refresh_projects()
  })

  selected_project <- reactive({
    catalog <- project_catalog()
    req(nrow(catalog), input$project_id)
    row <- catalog[catalog$id == input$project_id, , drop = FALSE]
    validate(need(nrow(row) == 1, "Project is not available"))
    row
  })

  manifest <- reactive({
    row <- selected_project()
    jsonlite::fromJSON(row$manifest_path[[1]], simplifyVector = FALSE)
  })

  samples <- reactive({
    row <- selected_project()
    read_manifest_table(row$project_root[[1]], manifest(), "samples", max_rows)
  })
  qc <- reactive({
    row <- selected_project()
    read_manifest_table(row$project_root[[1]], manifest(), "sample_qc", max_rows)
  })
  mapping <- reactive({
    row <- selected_project()
    read_manifest_table(row$project_root[[1]], manifest(), "mapping_qc", max_rows)
  })

  output$project_details <- renderUI({
    item <- manifest()$project
    tagList(
      strong(item$display_name), br(),
      span(class = "badge text-bg-dark", item$analysis_mode),
      span(class = "badge text-bg-success ms-1", item$status),
      if (nzchar(item$description %||% "")) tagList(br(), small(item$description))
    )
  })

  output$sample_count <- renderText({ nrow(samples()) })
  output$read_count <- renderText({
    data <- qc()
    if (!"reads_after" %in% names(data)) return("N/A")
    format_count(sum(data$reads_after, na.rm = TRUE))
  })
  output$mean_q30 <- renderText({
    data <- qc()
    if (!"q30_percent" %in% names(data) || all(is.na(data$q30_percent))) return("N/A")
    paste0(round(mean(data$q30_percent, na.rm = TRUE), 1), "%")
  })
  output$mean_mapping <- renderText({
    data <- mapping()
    if (!"mapping_percent" %in% names(data) || all(is.na(data$mapping_percent))) return("N/A")
    paste0(round(mean(data$mapping_percent, na.rm = TRUE), 1), "%")
  })

  output$project_summary <- renderUI({
    item <- manifest()
    tags$dl(
      tags$dt("Project ID"), tags$dd(item$project$id),
      tags$dt("Analysis"), tags$dd(item$project$analysis_mode),
      tags$dt("Generated"), tags$dd(item$generated_at),
      tags$dt("Pipeline"), tags$dd(paste(item$pipeline$name, item$pipeline$version)),
      tags$dt("HPC user"), tags$dd(current_hpc_user())
    )
  })

  output$primary_figure <- renderImage({
    item <- manifest()
    relative <- item$figures$pca %||% item$figures$sample_distance %||% ""
    validate(need(nzchar(relative), "No primary figure for this analysis mode"))
    row <- selected_project()
    path <- safe_project_file(row$project_root[[1]], relative, item)
    list(src = path, contentType = "image/png", alt = "Project overview figure")
  }, deleteFile = FALSE)

  output$qc_plot <- renderPlotly({
    data <- qc()
    validate(need(nrow(data), "QC metrics are not available"))
    metrics <- intersect(c("q30_percent", "gc_percent", "duplication_percent"), names(data))
    validate(need(length(metrics), "No plottable QC columns"))
    long <- tidyr::pivot_longer(data, all_of(metrics), names_to = "metric", values_to = "percent")
    plot_ly(long, x = ~sample_id, y = ~percent, color = ~metric, type = "bar") |>
      layout(barmode = "group", yaxis = list(title = "Percent"), xaxis = list(title = "Sample"))
  })
  output$qc_table <- renderDT(datatable(qc(), options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE))

  output$mapping_plot <- renderPlotly({
    data <- mapping()
    validate(need(nrow(data), "Mapping metrics are not available for this mode"))
    metric <- if ("mapping_percent" %in% names(data)) "mapping_percent" else "mean_coverage"
    validate(need(metric %in% names(data), "No plottable mapping or coverage metric"))
    plot_ly(data, x = ~sample_id, y = data[[metric]], type = "bar") |>
      layout(xaxis = list(title = "Sample"), yaxis = list(title = metric))
  })
  output$mapping_table <- renderDT(datatable(mapping(), options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE))

  observe({
    entries <- manifest()$differential_expression %||% list()
    choices <- if (length(entries)) {
      setNames(vapply(entries, `[[`, character(1), "id"), vapply(entries, `[[`, character(1), "label"))
    } else character()
    selected <- if (length(choices)) unname(choices[1]) else NULL
    updateSelectInput(session, "contrast_id", choices = choices, selected = selected)
  })

  selected_contrast <- reactive({
    entries <- manifest()$differential_expression %||% list()
    req(input$contrast_id, length(entries))
    matches <- Filter(function(item) identical(item$id, input$contrast_id), entries)
    validate(need(length(matches) == 1, "Contrast is unavailable"))
    matches[[1]]
  })

  output$volcano <- renderImage({
    entry <- selected_contrast()
    validate(need(nzchar(entry$volcano %||% ""), "Volcano plot is unavailable"))
    row <- selected_project()
    path <- safe_project_file(row$project_root[[1]], entry$volcano, manifest())
    list(src = path, contentType = "image/png", alt = paste(entry$label, "volcano plot"))
  }, deleteFile = FALSE)

  output$de_table <- renderDT({
    entry <- selected_contrast()
    row <- selected_project()
    path <- safe_project_file(row$project_root[[1]], entry$table, manifest())
    table <- readr::read_tsv(path, n_max = max_rows, show_col_types = FALSE, progress = FALSE)
    datatable(table, options = list(pageLength = 20, scrollX = TRUE), rownames = FALSE)
  })

  downloads <- reactive({
    entries <- manifest()$downloads %||% list()
    if (!length(entries)) return(data.frame(id = character(), label = character(), kind = character(), path = character()))
    data.frame(
      id = as.character(seq_along(entries)),
      label = vapply(entries, `[[`, character(1), "label"),
      kind = vapply(entries, `[[`, character(1), "kind"),
      path = vapply(entries, `[[`, character(1), "path"),
      stringsAsFactors = FALSE
    )
  })

  observe({
    data <- downloads()
    choices <- if (nrow(data)) setNames(data$id, data$label) else character()
    selected <- if (length(choices)) unname(choices[1]) else NULL
    updateSelectInput(session, "download_id", choices = choices, selected = selected)
  })
  output$download_table <- renderDT({
    datatable(downloads()[c("label", "kind")], options = list(dom = "t", pageLength = 50), rownames = FALSE)
  })
  output$download_result <- downloadHandler(
    filename = function() {
      data <- downloads()
      req(input$download_id)
      row <- data[data$id == input$download_id, , drop = FALSE]
      basename(row$path[[1]])
    },
    content = function(file) {
      data <- downloads()
      row_data <- data[data$id == input$download_id, , drop = FALSE]
      validate(need(nrow(row_data) == 1, "Download is unavailable"))
      project <- selected_project()
      source <- safe_project_file(project$project_root[[1]], row_data$path[[1]], manifest())
      validate(need(file.info(source)$size <= max_download_bytes, "File exceeds the portal download limit"))
      ok <- file.copy(source, file, overwrite = TRUE)
      if (!ok) stop("Unable to stage download")
    }
  )

  output$provenance <- renderText({
    item <- manifest()
    paste(
      paste("Project:", item$project$id),
      paste("Generated:", item$generated_at),
      paste("Pipeline:", item$pipeline$name, item$pipeline$version),
      paste("Manifest schema:", item$schema_version),
      paste("Portal root:", portal_root),
      paste("Unix user:", current_hpc_user()),
      sep = "\n"
    )
  })
}

shinyApp(ui, server)
