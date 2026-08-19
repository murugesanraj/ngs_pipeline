args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "app/run_app.R"
app_dir <- normalizePath(dirname(script_path), mustWork = TRUE)

host <- Sys.getenv("NGS_PORTAL_HOST", unset = "127.0.0.1")
port <- as.integer(Sys.getenv("NGS_PORTAL_PORT", unset = "0"))
if (is.na(port) || port < 0 || port > 65535) stop("NGS_PORTAL_PORT must be 0-65535")

setwd(app_dir)
shiny::runApp(
  appDir = app_dir,
  host = host,
  port = port,
  launch.browser = interactive(),
  display.mode = "normal"
)

