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
#' A dataset containing records the grape yield harvested in each row of a vinyard
#' in three succes years
#'
#' @format A data frame with 52 rows and 4 variables:
#' \describe{
#'   \item{row}{the vineyard row number}
#'   \item{year1}{reporting the harvest yield in first year}
#'   \item{year2}{reporting the harvest yield in second year}
#'   \item{year3}{reporting the harvest yield in third year}
#' }
#' @source 
"vineyard"


#' Response of 5 different monkey-pairs
#'
#' A dataset containing reports responses to a certain stimulus 
#' that were measured for 5 different monkey-pairs (the subjects) 
#' in 5 different periods under 5 different conditions
#'
#' @format A data frame with 25 rows and 4 variables:
#' \describe{
#'   \item{Cond.n}{the condition}
#'   \item{Monkeys}{the monkey pair}
#'   \item{Period}{the monkey pair}
#'   \item{Response}{responses to a certain stimulus}
#' }
#' @source p. 189 of Scheffes text and reformatted in monkey.RData
"monkey"

Seheult and Tukey (2001) analyzed a three-factor layout in which the response measures
the hardness of dental fillings obtained by 5 Dentists (D) using 8 Gold alloys (G) and
3 Condensation methods (C). The objective of the experiment was to find a dental gold
lling with greater hardness. Condensation, properly carried out, was known to increase
the hardness of a filling. The three condensation techniques used in the experiment were:
  (1) electromalleting, in which blows are delivered mechanically at a steady frequency; (2)
hand malleting, in which a small mallet is used to deliver blows; and (3) hand condensation.
The reported hardness observations are each averages of ten measurements that are not
available. It was reported anecdotally that dentist 5 appeared to be physically tired before
the experiment.

#' The hardness of 120 dental fillings
#'
#' A dataset containing the response measures the hardness of dental filling
#' obtained by 5 Dentists using 8 Gold alloys and 3 Condensation methods.
#' The objective of the experiment was to find a dental gold filling with greater hardness.
#' 
#' @format A data frame with 120 rows and 4 variables:
#' \describe{
#'   \item{y}{response mesures the hardness of dental fillings}
#'   \item{G}{the indice of 8 Gold alloys}
#'   \item{C}{the indice of 3 Condensation methods}
#'   \item{D}{the indice of 5 Dentists}
#' }
#' @source Seheult and Tukey (2001)
"dental"


#' Accelations over time
#'
#' A dataset containing 133 observation of motorsycle
#' acceleration against time in a simulated motorcycle accident.
#' The p = 277 possible observation times constitute the vector t = (1, 2, ... , 277)
#' Accelerations were observed at only q < p of these equally spaced time,
#' sometimes with replication.
#'
#' @format A data frame with 133 rows and 2 variables:
#' \describe{
#'   \item{t}{possible observation times}
#'   \item{accel}{acceleration against time}
#' }
#' @source adapted from Silverman (1985)
"motor"

