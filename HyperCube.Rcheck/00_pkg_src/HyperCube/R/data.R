#' Weigth gain of 61 infant rat litters
#'
#' A dataset containing the (average)  wight gain of an infant rat litter
#' when the infants in the litter are nursed by a rat foster-mother
#' 
#' The rat litter data treated by Scheffe (1959) form an unbalnced two-way layout
#' Each response recorded is the average weight-gain of a rat litter when the infants
#' in the litter are nursed by a rat foster-mother. Factor1 with four levels, is the 
#' genotype of the foster-moather. Factor2 with the same levels, in the genotype of the 
#' infant litter.
#' 
#' The response measured in the experiment is the (average) weight gain
#' of an infant rat litter when the infants in the litter are nursed by a rat foster-mother. 
#' Factor 1 is the genotype of the foster-mother nursing the infants. 
#' Factor 2 is the genotype of the infant litter.
#'
#' @format A data frame with 61 rows and 3 variables:
#' \describe{
#'   \item{weight}{the (averge) weight gain of an infant rat litter}
#'   \item{mother}{the genotype of the foster-mother nursing the infants}
#'   \item{infant}{the genotype of the infant litter}
#' }
#' @source D.W. Bailey, The Inheritance of Maternal Influences on the Growth of the Rat (1953).
"litter"

#' Vineyard
#'
#' A dataset containing records the grape yield harvested in each row of a vinyard
#' in three succeed years
#'
#' @format A data frame with 52 rows and 4 variables:
#' \describe{
#'   \item{row}{the vineyard row number}
#'   \item{year1}{reporting the harvest yield in first year}
#'   \item{year2}{reporting the harvest yield in second year}
#'   \item{year3}{reporting the harvest yield in third year}
#' }
#' @source Simonoff, Jeffrey S. Smoothing methods in statistics. Springer Science & Business Media, 2012.
"vineyard"

#' Response of 5 different monkey-pairs
#'
#' A dataset containing reports responses to a certain stimulus 
#' that were measured for 5 different monkey-pairs (the subjects) 
#' in 5 different periods under 5 different conditions
#'
#' @format A data frame with 25 rows and 4 variables:
#' \describe{
#'   \item{cond}{the condition}
#'   \item{monkeys}{the monkey pair}
#'   \item{period}{the monkey pair}
#'   \item{response}{responses to a certain stimulus}
#' }
#' @source Data from Query no. 113, edited by G.W. Sender, Biometrics, Vol. 11, 1955, p.112
"monkey"


#' The hardness of 120 dental fillings
#'
#' A dataset containing the response measures the hardness of dental filling
#' obtained by 5 Dentists using 8 Gold alloys and 3 Condensation methods.
#' The objective of the experiment was to find a dental gold filling with greater hardness.
#' 
#' 
#' Seheult and Tukey (2001) analyzed a three-factor layout in which the response measures
#' the hardness of dental fillings obtained by 5 Dentists (D) using 8 Gold alloys (G) and
#' 3 Condensation methods (C). The objective of the experiment was to find a dental gold
#' lling with greater hardness. Condensation, properly carried out, was known to increase
#' the hardness of a filling. The three condensation techniques used in the experiment were:
#'   (1) electromalleting, in which blows are delivered mechanically at a steady frequency; (2)
#' hand malleting, in which a small mallet is used to deliver blows; and (3) hand condensation.
#' The reported hardness observations are each averages of ten measurements that are not
#' available. It was reported anecdotally that dentist 5 appeared to be physically tired before
#' the experiment.
#' 
#' @format A data frame with 120 rows and 4 variables:
#' \describe{
#'   \item{y}{response mesures the hardness of dental fillings}
#'   \item{G}{the indice of 8 Gold alloys}
#'   \item{C}{the indice of 3 Condensation methods}
#'   \item{D}{the indice of 5 Dentists}
#' }
#' @source Seheult, A. H. and Tukey, J. W. (2001). Towards robust analysis of variance, Data Analysis from Statistical Foundations (eds. A. K. Mohammed and E. Saleh), 217-244, Nova Science Publishers, New York.
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
#'   \item{t}{The time in milliseconds since impact}
#'   \item{accel}{The recorded head acceleration (in g)}
#' }
#' @source Silverman,  B.W.  (1985)  Some  aspects  of   the   spline smoothing   approach   to  non-parametric  curve  fitting. Journal of the Royal Statistical Society, B, 47, 1-52.
"motor"

#' Canadian earnings
#' 
#' The canadian.earnings data frame has 205 pairs observations on Canadian workers from a 1971 Canadian Census Public Use Tape (Ullah, 1985).
#' 
#' @format A data frame with 205 rows and 2 variables:
#' \describe{
#'   \item{age}{age in years.}
#'   \item{log.income}{logarithm of income.}
#' }
#' @source Ullah, A. (1985). Specification analysis of econometric models. Journal of Quantitative Economics, 2, 187-209
"canadian.earnings"
