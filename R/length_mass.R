#' Converts length measures to body mass for bees and Carabids
#' 
#' For Bees, the conversion is done as per Cane et al. 1987 (Estimation of bee 
#' size using intertegular span (Apoidea). Journal of the Kansas Entomological 
#' Society 60:145–147). The correlation between IT and body mass is 0.96.
#' For Carabids the conversion is done from Jelaska et al. 2011 (Carabid beetle 
#' diversity and mean individual biomass in beech forests of various ages. 
#' ZooKeys 100:393-405) which in turn is based on Szyszko J (1983) book.
#'
#' @param x A vector of bee intertegular spans (IT) measurments in cm 
#' or carabid beetles body length (mm)\code{}
#'
#' @return A vector of bee body masses is returned in gr for bees and (mg) for Carabids.
#'
#' @keywords Bee, Carabids
#'
#' @export
#' 
#' @examples
#' Cane(c(1.2, 2.3, 0.6))
#' Jelaska(c(1.2, 2.3, 0.6))

Cane <- function(x){exp(0.6453 + 2.4691*log(x))}

Jelaska <- function(x){exp(-8.92804283 + 2.5554921 * log(x))}
