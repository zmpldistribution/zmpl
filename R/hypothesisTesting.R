#'Dispersion test
#'
#'@description Tests the null hypothesis of precision fixed in RBS models against the alternative of precision variable.
#'
#'@usage grad.test(modelh0,modelh1)
#'
#'
#' @param modelh0 model under null hypothesis.
#' @param modelh1 model under alternative hypothesis.
#'
#' @return A list with class "htest" containing the following components:
#' @return \code{statistic}	the value of the test statistic.
#' @return \code{parameter}	the degrees of freedom for the test statistic.
#' @return \code{p.value}	the p-value for the test.
#' @return \code{method}	a character string indicating what type of likelihood ratio test was performed.
#' @return \code{data.name} a character string giving the name(s) of the data
#'
#'@author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'@examples
#'
#'##
#'library(alr3)
#'data(landrent)
#'attach(landrent)
#'resp <- I(Y/X1)
#'y1 <-  split(resp, X4)$"1"
#'x21 <-  split(X2, X4)$"1"
#'##Fixed Precision
#'fit0 <- gamlss(y1 ~ x21, family=RBS(mu.link="identity"),method=CG()  )
#'##Varying Precision
#'fit1 <- gamlss(y1 ~ x21,sigma.formula = y1 ~x21, family=RBS(mu.link="identity",sigma.link="sqrt"),method=CG()  )
#'#Precision Test
#'lr.test(fit0,fit1)
#'score.test(fit0,fit1)
#'grad.test(fit0,fit1)
#'wald.test(fit1)
#'@export
#'
grad.test <- function(x,lambda)
{
  DNAME <- deparse(substitute(x))
  METHOD <- "Grandient test"
  n<- length(x)
  piest <- lambda[1]
  thetaest <- lambda[2]
  STATISTIC <- n*piest*piest*( thetaest*thetaest + 3*thetaest +1)
  gl = 1
  PVAL = pchisq(STATISTIC, 1, lower.tail = F)
  names(gl) = "df"
  names(STATISTIC ) = "G"
  RVAL <- list(statistic = STATISTIC, parameter = gl, p.value = PVAL,
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}
