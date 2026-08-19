read_app_config <- function(app_dir) {
  configured <- Sys.getenv("NGS_APP_CONFIG", unset = "")
  candidates <- c(
    configured,
    file.path(app_dir, "config.yml"),
    file.path(app_dir, "..", "config", "app.example.yml")
  )
  candidates <- candidates[nzchar(candidates)]
  config_path <- candidates[file.exists(candidates)][1]
  if (is.na(config_path)) stop("No app configuration file found")
  cfg <- yaml::read_yaml(config_path)
  env_root <- Sys.getenv("NGS_PORTAL_ROOT", unset = "")
  if (nzchar(env_root)) cfg$paths$portal_root <- env_root
  cfg$paths$portal_root <- normalizePath(
    path.expand(cfg$paths$portal_root), mustWork = TRUE
  )
  cfg
}

current_hpc_user <- function() {
  value <- Sys.info()[["user"]]
  if (is.null(value) || !nzchar(value)) value <- Sys.getenv("USER", unset = "unknown")
  value
}

discover_projects <- function(portal_root) {
  entries <- list.files(portal_root, full.names = TRUE, all.files = FALSE, no.. = TRUE)
  entries <- entries[dir.exists(entries)]
  records <- lapply(entries, function(entry) {
    manifest_path <- file.path(entry, "portal", "manifest.json")
    if (!file.exists(manifest_path) || file.access(manifest_path, 4) != 0) return(NULL)
    manifest <- tryCatch(jsonlite::fromJSON(manifest_path, simplifyVector = FALSE), error = function(e) NULL)
    if (is.null(manifest) || is.null(manifest$project$id)) return(NULL)
    project_root <- normalizePath(entry, mustWork = TRUE)
    data.frame(
      id = as.character(manifest$project$id),
      name = as.character(manifest$project$display_name %||% manifest$project$id),
      mode = as.character(manifest$project$analysis_mode %||% "unknown"),
      generated_at = as.character(manifest$generated_at %||% ""),
      project_root = project_root,
      manifest_path = normalizePath(manifest_path, mustWork = TRUE),
      stringsAsFactors = FALSE
    )
  })
  records <- records[!vapply(records, is.null, logical(1))]
  if (!length(records)) {
    return(data.frame(
      id = character(), name = character(), mode = character(),
      generated_at = character(), project_root = character(),
      manifest_path = character(), stringsAsFactors = FALSE
    ))
  }
  result <- do.call(rbind, records)
  result <- result[!duplicated(result$id), , drop = FALSE]
  result[order(tolower(result$name)), , drop = FALSE]
}

`%||%` <- function(left, right) {
  if (
    is.null(left) || length(left) == 0 ||
      (length(left) == 1 && (is.na(left) || !nzchar(as.character(left))))
  ) right else left
}

manifest_paths <- function(manifest) {
  paths <- character()
  append_value <- function(value) {
    if (is.character(value) && length(value) == 1 && nzchar(value)) paths <<- c(paths, value)
  }
  lapply(manifest$tables %||% list(), append_value)
  lapply(manifest$figures %||% list(), append_value)
  for (item in manifest$differential_expression %||% list()) {
    lapply(item[c("table", "volcano", "ma_plot")], append_value)
  }
  for (item in manifest$assemblies %||% list()) {
    lapply(item[c("contigs", "quast_html", "quast_table")], append_value)
  }
  for (item in manifest$downloads %||% list()) append_value(item$path)
  unique(paths)
}

safe_project_file <- function(project_root, relative_path, manifest) {
  if (!is.character(relative_path) || length(relative_path) != 1 || !nzchar(relative_path)) {
    stop("Missing result path")
  }
  if (grepl("^/|(^|/)\.\.(/|$)|[\r\n]", relative_path)) stop("Unsafe result path")
  if (!relative_path %in% manifest_paths(manifest)) stop("Result is not in the manifest allowlist")
  root <- normalizePath(project_root, mustWork = TRUE)
  candidate <- normalizePath(file.path(root, relative_path), mustWork = TRUE)
  prefix <- paste0(root, .Platform$file.sep)
  if (!(identical(candidate, root) || startsWith(candidate, prefix))) {
    stop("Result resolves outside the project directory")
  }
  if (file.access(candidate, 4) != 0) stop("Result is not readable by the current HPC user")
  candidate
}

read_manifest_table <- function(project_root, manifest, key, max_rows = 50000) {
  relative <- manifest$tables[[key]]
  if (is.null(relative) || !nzchar(relative)) return(data.frame())
  path <- safe_project_file(project_root, relative, manifest)
  readr::read_tsv(path, n_max = max_rows, show_col_types = FALSE, progress = FALSE)
}

format_count <- function(value) {
  if (!length(value) || is.na(value)) return("N/A")
  suffix <- if (value >= 1e9) "B" else if (value >= 1e6) "M" else if (value >= 1e3) "K" else ""
  divisor <- if (value >= 1e9) 1e9 else if (value >= 1e6) 1e6 else if (value >= 1e3) 1e3 else 1
  paste0(format(round(value / divisor, 2), trim = TRUE), suffix)
}
