#' Map functions
#' 
#' Computes cM map distance from recombination frequency
#' 
#' @param r recombination frequency
#' @param model Either "Haldane" or "Kosambi" 
#' 
#' @return Map distance in cM
#' @export

map_fn <- function(r,model) {
  model <- toupper(model)
  if (length(r) > 1) {
    x <- apply(array(r),1,map_fn,model=model)/100
  } else {
    if (0 <= r & r < 0.5) {
      if (model=="HALDANE") {x <- -1/2*log(1-2*r)}
      if (model=="KOSAMBI") {x <- 1/4*log((1+2*r)/(1-2*r))}
    } else {
      x <- NA
    }
  }
  return(x*100) #return cM
}