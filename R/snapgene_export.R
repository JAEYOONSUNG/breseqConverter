# SnapGene's documented CLI requires that the application is not already open:
# https://support.snapgene.com/hc/en-us/articles/10384393330836-Command-Line-Converting-File-Formats
# Never close an existing session or terminate a delayed conversion process.
.breseq_snapgene_is_running <- function() {
  tryCatch({
    if (.Platform$OS.type == "windows") {
      result <- suppressWarnings(system2("tasklist", c("/FO", "CSV", "/NH"),
                                        stdout = TRUE, stderr = FALSE))
      pattern <- '^"SnapGene( Viewer)?[.]exe"'
    } else {
      result <- suppressWarnings(system2("ps", c("-A", "-o", "comm="),
                                        stdout = TRUE, stderr = FALSE))
      pattern <- "(^|/)SnapGene( Viewer)?([.]exe)?$"
    }
    status <- attr(result, "status")
    if (!is.null(status) && status != 0L) return(NA)
    any(grepl(pattern, trimws(result), ignore.case = TRUE))
  }, error = function(e) NA)
}

.breseq_snapgene_run <- function(executable, args, timeout) {
  if (!requireNamespace("processx", quietly = TRUE)) {
    return(list(unavailable = TRUE, pending = FALSE, status = NA_integer_,
                message = "Install the optional processx package for SnapGene conversion."))
  }
  # cleanup=FALSE is intentional: a timeout must not kill SnapGene or a prompt.
  process <- processx::process$new(executable, args, stdout = "|", stderr = "|",
                                  cleanup = FALSE, cleanup_tree = FALSE)
  process$wait(timeout = timeout * 1000)
  if (process$is_alive()) {
    return(list(unavailable = FALSE, pending = TRUE, status = NA_integer_,
                message = paste0("SnapGene did not exit within ", timeout,
                  " seconds; its process was left running and was not terminated.")))
  }
  list(unavailable = FALSE, pending = FALSE, status = process$get_exit_status(),
       message = trimws(paste(process$read_all_output(), process$read_all_error())))
}

.breseq_snapgene_dna_header <- function(file) {
  if (!file.exists(file) || file.info(file)$size <= 19L) return(FALSE)
  connection <- file(file, "rb")
  on.exit(close(connection), add = TRUE)
  header <- readBin(connection, "raw", n = 19L)
  # Cookie packet: type 9, size 14, signature SnapGene, sequence type DNA (1).
  length(header) == 19L &&
    identical(as.integer(header[1:5]), c(9L, 0L, 0L, 0L, 14L)) &&
    identical(header[6:13], charToRaw("SnapGene")) &&
    identical(as.integer(header[14:15]), c(0L, 1L))
}

# Return an explicit conversion status; never fabricate a .dna by renaming text.
# A single edited GenBank replicon is expected.  The caller retains that GenBank
# whether native SnapGene conversion succeeds, is unavailable, or is deferred.
.breseq_convert_genbank_to_dna <- function(input, output, snapgene_path = NULL,
                                          timeout = 60) {
  valid_path <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  if (!valid_path(input) || !file.exists(input) || dir.exists(input))
    stop("input must name an existing GenBank file.")
  if (!valid_path(output) || !grepl("[.]dna$", output, ignore.case = TRUE))
    stop("output must name a SnapGene .dna file.")
  if (!is.numeric(timeout) || length(timeout) != 1L || !is.finite(timeout) || timeout <= 0)
    stop("timeout must be a positive number of seconds.")
  input <- normalizePath(input, mustWork = TRUE)
  lines <- readLines(input, warn = FALSE)
  if (sum(grepl("^LOCUS[[:space:]]", lines)) != 1L ||
      sum(grepl("^ORIGIN([[:space:]]|$)", lines)) != 1L ||
      sum(grepl("^//[[:space:]]*$", lines)) != 1L)
    stop("SnapGene conversion requires one complete GenBank sequence record.")
  result <- function(status, message, success = FALSE, path = NULL, exit_status = NA_integer_) {
    list(success = success, status = status, output_file = path,
         message = message, exit_status = exit_status)
  }
  if (file.exists(output))
    return(result("output_exists", "The output already exists and was not overwritten."))
  if (!is.null(snapgene_path) &&
      (!valid_path(snapgene_path) || !file.exists(snapgene_path) || dir.exists(snapgene_path)))
    return(result("unavailable", "The configured SnapGene executable is unavailable."))
  executable <- tryCatch(suppressMessages(find_snapgene(snapgene_path)), error = function(e) NULL)
  if (is.null(executable) || file.access(executable, mode = 1L) != 0L)
    return(result("unavailable", "A usable SnapGene executable was not found; the GenBank file is retained."))
  running <- .breseq_snapgene_is_running()
  if (is.na(running))
    return(result("process_check_failed", "Could not determine whether SnapGene is running; conversion was not started."))
  if (running)
    return(result("application_running", "SnapGene is already running; conversion was deferred without closing it."))

  directory <- dirname(path.expand(output))
  if (!dir.exists(directory) && !dir.create(directory, recursive = TRUE))
    stop("Could not create output directory: ", directory)
  output <- file.path(normalizePath(directory, mustWork = TRUE), basename(output))
  temporary <- tempfile(".snapgene-conversion-", tmpdir = dirname(output), fileext = ".dna")
  remove_temporary <- TRUE
  on.exit(if (remove_temporary && file.exists(temporary)) unlink(temporary), add = TRUE)
  args <- c("--quiet", "--convert", "SnapGene", "--input", input, "--output", temporary)
  converted <- tryCatch(.breseq_snapgene_run(executable, args, timeout),
    error = function(e) list(unavailable = FALSE, pending = FALSE,
                             status = NA_integer_, message = conditionMessage(e)))
  if (isTRUE(converted$unavailable)) return(result("unavailable", converted$message))
  if (isTRUE(converted$pending)) {
    remove_temporary <- FALSE
    pending <- result("pending", converted$message)
    pending$pending_file <- temporary
    return(pending)
  }
  if (is.na(converted$status) || converted$status != 0L)
    return(result("conversion_failed", converted$message, exit_status = converted$status))
  if (!.breseq_snapgene_dna_header(temporary))
    return(result("invalid_output", "SnapGene did not produce a valid DNA-format header.",
                  exit_status = converted$status))
  if (file.exists(output) || !file.copy(temporary, output, overwrite = FALSE))
    return(result("output_exists", "The destination could not be created without overwriting an existing file."))
  result("converted", "Native SnapGene DNA conversion completed and its binary header was verified.",
         success = TRUE, path = normalizePath(output), exit_status = converted$status)
}
