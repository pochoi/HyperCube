#' foo: A package for computating the notorious bar statistic.
#'
#' The foo package provides three categories of important functions:
#' foo, bar and baz.
#' 
#' @section Foo functions:
#' The foo functions ...
#'
#' @docType package
#' @name foo
NULL
#> NULL

The rat litter data, printed on p. 140 of H. Scheff ́e’s text,
comes from a Ph.D. thesis The Inheritance of Maternal Influences on the Growth of the Rat 
by D. W. Bailey (1953). 
The response measured in the experiment is the (average) weight gain
of an infant rat litter when the infants in the litter are nursed by a rat foster-mother. 
Factor 1 is the genotype of the foster-mother nursing the infants. 
Factor 2 is the genotype of the infant litter.

#' Weigth gain of 61 infant rat litters
#'
#' A dataset containing the (average)  wight gain of an infant rat litter
#' when the infants in the litter are nursed by a rat foster-mother
#' 
#'
#' @format A data frame with 61 rows and 3 variables:
#' \describe{
#'   \item{weight}{the (averge) weight gain of an infant rat litter}
#'   \item{motherGen}{the genotype of the foster-mother nursing the infants}
#'   \item{infantGen}{the genotype of the infant litter}
#' }
#' @source printed on p. 140 of H. Scheff ́e’s text,comes from a Ph.D. thesis
#' The Inheritance of Maternal Influences on the Growth of the Rat by D. W. Bailey (1953).
"litter"

#' Prices of 50,000 round cut diamonds.
#'
#' A dataset containing the prices and other attributes of almost 54,000
#' diamonds.
#'
#' @format A data frame with 53940 rows and 10 variables:
#' \describe{
#'   \item{price}{price, in US dollars}
#'   \item{carat}{weight of the diamond, in carats}
#'   ...
#' }
#' @source \url{http://www.diamondse.info/}
"diamonds"