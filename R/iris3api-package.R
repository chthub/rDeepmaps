#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
#' @importFrom tibble tibble
#' @importFrom glue glue
#' @importFrom S4Vectors isTRUEorFALSE
#' @import Seurat
#' @import Signac
## usethis namespace: end
e1 <- rlang::env(
  jobid = NULL, # job id
  type = NULL, # job type: RNA or Matched?
  species = NULL,
  obj = NULL, # active object
  sub_obj = NULL, # the subset object from selection
  full_obj = NULL, # object with all cells
  meta = NULL, # uploaded metadata
  assay_idx = 1,
  embedding_idx = 2,
  ident_idx = 1,
  new_meta_counter = 1,
  deg = NULL
)

e2<- rlang::env(
  sio = NULL #socketio
)

NULL
