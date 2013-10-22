#' Converts bee intertegular distance (IT; in cm) to body mass (gr)
#' 
#' The conversion is done as per Cane et al. 1996¿? formula and the correlation is >0.96.
#'
#' @param IT A vector of IT measurments \code{IT}
#'
#' @return output A vector of bee body masses 
#'
#' @keywords Bee
#'
#' @export
#' 
#' @examples
#' Cane(c(1.2, 2.3, 0.6))

Cane <- function(IT){exp(0.6453 + 2.4691*log(IT))}
