#' Hello world function
#'
#' This function sends hello message from a given name.
#'
#' @param name echo name
#'
#' @return A string
#' @export
hello_name <- function(name) {
  result <- as.character(stringr::str_glue("Hello, {name}"))
  message(result)
  return(result)
}


#' DEV tool: auto style and lint package
#'
#' @return NULL
#' @export
#'
style_and_lint <- function() {
  styler::style_pkg()
  lintr::lint_package()
}

#' Set a default value if an object is null
#'
#' @param lhs An object to set if it's null
#' @param rhs The value to provide if x is null
#'
#' @return rhs if lhs is null, else lhs
#'
#' @author Hadley Wickham
#' @references https://adv-r.hadley.nz/functions.html#missing-arguments
#'
set_if_null <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}

#' Set a default value if an object is NOT null
#'
#' @param lhs An object to set if it's NOT null
#' @param rhs The value to provide if x is NOT null
#'
#' @return lhs if lhs is null, else rhs
#'
#' @author Hadley Wickham
#' @references https://adv-r.hadley.nz/functions.html#missing-arguments
#'
set_if_not_null <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(rhs)
  } else {
    return(lhs)
  }
}

#' Check to ensure a package is installed
#'
#' @param package Name of pacakge to check
#' @param repository Repository that package is available on;
#' choose from 'bioconductor', 'github', or 'cran'
#' @param ... Extra parameters passed to BiocManager::install,
#' remotes::install_github, or install.packages, depending on \code{repository}
#'
#' @importFrom utils menu install.packages
#' @author Paul Hoffman
#' @references https://github.com/satijalab/seurat-wrappers/wiki/Code-Guidelines
check_package <- function(package, repository, ...) {
  if (!requireNamespace(package = basename(path = package), quietly = TRUE)) {
    if (interactive()) {
      message("Package ", package, " is not yet installed")
      message("Install now?")
      choice <- menu(choices = c("yes", "no"))
      if (choice == 1) {
        repository <- match.arg(
          arg = tolower(x = repository),
          choices = c("github", "bioconductor", "cran")
        )
        switch(
          EXPR = repository,
          "github" = remotes::install_github(repo = package, ...),
          "bioconductor" = BiocManager::install(pkgs = package, ...),
          "cran" = install.packages(pkgs = package, ...),
          stop("Unknown repository ", repository, call. = FALSE)
        )
        return(invisible(x = NULL))
      }
    }
    stop("Unable to find package ", package, ", please install", call. = FALSE)
  }
}

#' Check if a matrix is empty
#'
#' Takes a matrix and asks if it's empty (either 0x0 or 1x1 with a value of NA)
#'
#' @param x A matrix
#'
#' @return Whether or not \code{x} is empty
#' @author Paul Hoffman
#' @references https://github.com/satijalab/seurat-wrappers/wiki/Code-Guidelines
#'
is_matrix_empty <- function(x) {
  matrix.dims <- dim(x = x)
  matrix.na <- all(matrix.dims == 1) && all(is.na(x = x))
  return(all(matrix.dims == 0) || matrix.na)
}


#' Convert empty string to "-"
#'
#' @param string string
#'
#' @return
#'
#' @references https://rviews.rstudio.com/2019/08/13/plumber-logging/
convert_empty <- function(str) {
  if (str == "") {
    "-"
  } else {
    str
  }
}
