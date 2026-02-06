#' Find SnapGene Executable Path
#'
#' Automatically detects SnapGene installation across different operating systems.
#' Supports macOS, Windows, and Linux with common installation paths.
#'
#' @param custom_path Optional custom path to SnapGene executable
#' @return Character string with the path to SnapGene executable
#' @export
#' @examples
#' \dontrun{
#' # Auto-detect
#' snapgene_path <- find_snapgene()
#'
#' # Use custom path
#' snapgene_path <- find_snapgene("/custom/path/to/snapgene")
#'
#' # Set via environment variable
#' Sys.setenv(SNAPGENE_PATH = "/path/to/snapgene")
#' snapgene_path <- find_snapgene()
#' }
find_snapgene <- function(custom_path = getOption("breseqConverter.snapgene_path")) {

  # 1. User-specified path takes priority
  if (!is.null(custom_path) && file.exists(custom_path)) {
    message("Using custom SnapGene path: ", custom_path)
    return(normalizePath(custom_path))
  }

  # 2. Check environment variable
  env_path <- Sys.getenv("SNAPGENE_PATH")
  if (nzchar(env_path) && file.exists(env_path)) {
    message("Using SnapGene from SNAPGENE_PATH: ", env_path)
    return(normalizePath(env_path))
  }

  # 3. OS-specific default paths
  os <- Sys.info()["sysname"]

  candidates <- switch(os,
    "Darwin" = c(
      "/Applications/SnapGene.app/Contents/MacOS/SnapGene",
      path.expand("~/Applications/SnapGene.app/Contents/MacOS/SnapGene"),
      "/Applications/SnapGene Viewer.app/Contents/MacOS/SnapGene Viewer",
      path.expand("~/Applications/SnapGene Viewer.app/Contents/MacOS/SnapGene Viewer")
    ),
    "Windows" = {
      prog_files <- Sys.getenv("PROGRAMFILES")
      prog_files_x86 <- Sys.getenv("PROGRAMFILES(X86)")
      local_app <- Sys.getenv("LOCALAPPDATA")
      c(
        file.path(prog_files, "SnapGene", "snapgene.exe"),
        file.path(prog_files_x86, "SnapGene", "snapgene.exe"),
        file.path(local_app, "Programs", "SnapGene", "snapgene.exe"),
        file.path(prog_files, "SnapGene Viewer", "snapgene.exe"),
        file.path(prog_files_x86, "SnapGene Viewer", "snapgene.exe")
      )
    },
    "Linux" = c(
      "/usr/bin/snapgene",
      "/usr/local/bin/snapgene",
      "/opt/snapgene/snapgene",
      "/opt/gslbiotech/snapgene/snapgene.sh",
      path.expand("~/snapgene/snapgene"),
      "/snap/bin/snapgene"
    ),
    character(0)
  )

  # 4. Check PATH
  path_result <- Sys.which("snapgene")
  if (nzchar(path_result)) {
    candidates <- c(unname(path_result), candidates)
  }

  # 5. Return first existing path

  for (path in candidates) {
    if (!is.na(path) && nzchar(path) && file.exists(path)) {
      message("SnapGene found: ", path)
      return(normalizePath(path))
    }
  }

  stop(
    "SnapGene not found. Please either:\n",
    "  1. Set environment variable: Sys.setenv(SNAPGENE_PATH = '/path/to/snapgene')\n",
    "  2. Set R option: options(breseqConverter.snapgene_path = '/path/to/snapgene')\n",
    "  3. Pass custom_path argument directly\n",
    "  4. Install SnapGene in the default location\n\n",
    "Searched locations:\n  ", paste(candidates, collapse = "\n  ")
  )
}

#' Check if SnapGene is available
#'
#' @return Logical indicating if SnapGene is found
#' @export
has_snapgene <- function() {
  tryCatch({
    find_snapgene()
    TRUE
  }, error = function(e) FALSE)
}
