.statgen_state <- new.env(parent = emptyenv())
.statgen_state$verbosity <- "info"

set_verbosity <- function(level) {
  if (!is.character(level) || length(level) != 1L || is.na(level)) {
    stop("level must be a character scalar", call. = FALSE)
  }
  if (!(level %in% c("quiet", "info"))) {
    stop("level must be one of 'quiet' or 'info'", call. = FALSE)
  }
  .statgen_state$verbosity <- level
  invisible(level)
}

get_verbosity <- function() {
  .statgen_state$verbosity
}
