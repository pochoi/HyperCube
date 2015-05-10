#' Hello world
#' 
#' @param name characters string
#' @examples
#' ## Say hello!
#' hello("Amy")
#' @export
hello <-
function(name) {
  print(paste0("Hello ", name, "!"))
}

#' Bye world
#' 
#' @param name characters string
#' @examples
#' ## Say Bye!
#' bye("Po")
#' @export
bye <-
function(name) {
  print(paste0("Good Bye ", name, "~"))
}