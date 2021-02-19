#' Inverse map function 
#' 
#' Computes recombination frequency from map distance
#' 
#' @param x map distance (cM)
#' @param model Either "Haldane" or "Kosambi" 
#' 
#' @return recombination frequency
#' @export

inverse_map_fn <- function(x,model) {
  model <- toupper(model)
  x <- x/100
  if (0 <= r & r < 0.5) {
    if (model=="HALDANE") {r <- (1-exp(-2*x))/2}
    if (model=="KOSAMBI") {r <- (exp(2*x)-exp(-2*x))/(exp(2*x)+exp(-2*x))/2}
  } else {
    r <- NA
  }
  return(r) 
}