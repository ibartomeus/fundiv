#' @name length_to_mass
#' @aliases Cane
#' @aliases Jelaska
#' 
#' @title Converts length measures to body mass for bees and Carabids
#' 
#' @description For Bees, the conversion is done as per Cane et al. 1987 (Estimation of bee 
#' size using intertegular span (Apoidea). Journal of the Kansas Entomological 
#' Society 60:145–147). The correlation between IT and body mass is 0.96.
#' For Carabids the conversion is done from Jelaska et al. 2011 (Carabid beetle 
#' diversity and mean individual biomass in beech forests of various ages. 
#' ZooKeys 100:393-405) which in turn is based on Szyszko J (1983) book.
#'
#' @param x A vector of bee intertegular spans (IT) measurments in cm 
#' or carabid beetles body length (mm)\code{}
#' @param taxa string character indicating if length measure is from "Bees" or "Carabids"
#'
#' @return A vector of bee body masses is returned in (gr) for bees and for Carabids.
#'
#' @keywords Bee, Carabids
#'
#' @rdname length_to_mass
#' @export
#' 
length_to_mass <- function(x, taxa = c("Bees", "Carabids")){
  if(taxa == "Bees")Out <- Cane(x)
  if(taxa == "Carabids")Out <- Jelaska(x)
  Out
}
#' @examples
#' length_to_mass(c(1.2, 2.3, 0.6), "Bees")
#'
#' @rdname length_to_mass
#' @examples
#' Cane(c(1.2, 2.3, 0.6))
#' @export
Cane <- function(x){exp(0.6453 + 2.4691*log(x))}
#'
#' @rdname length_to_mass
#' @examples
#' Jelaska(c(1.2, 2.3, 0.6))
#' @export
Jelaska <- function(x){exp(-8.92804283 + 2.5554921 * log(x))}
