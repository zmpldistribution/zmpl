#'Zero Modified Poisson-Lindley (ZMPL) distribution
#'
#'@description The fuction \code{ZMPL()} defines the ZMPL distribution, a two paramenter
#'distribution. The zero modified Poisson-Lindley distribution is similar
#'to the Poisson-Lindley distribution but allows zeros as y values. The extra parameter
#'models the probabilities at zero. The functions dZMPL, pZMPL, qZMPL and rZMPL define
#'the density, distribution function, quantile function and random generation for
#'the zero modified Poisson-Lindley distribution.
#'plotZMPL can be used to plot the distribution. meanZMPL calculates the expected
#'value of the response for a fitted model.
#'
#'@usage
#'dZMPL(x, mu = 1, sigma = 1, nu=0.1, log = FALSE)
#'pZMPL(q, mu = 1, sigma = 1,  nu=0.1, lower.tail = TRUE, log.p = FALSE)
#'qZMPL(p, mu = 1, sigma = 1,  nu=0.1, lower.tail = TRUE, log.p = FALSE)
#'rZMPL(n, mu = 1, sigma = 1)
#'plotZMPL(mu = .5, sigma = 1,  nu=0.1, from = 0, to = 0.999, n = 101, ...)
#'meanZMPL(obj)
#'
#' @param x,q vector of observations/quantiles
#' @param mu vector of scale parameter values
#' @param sigma vector of shape parameter values
#' @param log, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param from where to start plotting the distribution from
#' @param to up to where to plot the distribution
#' @param ... other graphical parameters for plotting
#'
#'
#' @details  The probability massa function of the is given by
#'
#'\deqn{ f_{Y}(y;\mu,\delta,p) =\frac{[1-p]\sqrt{\delta+1}}{4\,y^{3/2}\,\sqrt{\pi\mu}}\left[y+\frac{\delta\mu}{\delta+1} \right]\exp\left(-\frac{\delta}{4}\left[\frac{y[\delta+1]}{\delta\mu}+\frac{\delta\mu}{y[\delta+1]}-2\right]\right) I_{(0, \infty)}(y)+  pI_{\{0\}}(y).}
#'
#'@return returns a object of the ZMPL distribution.
#'
#'@note For the function ZMPL(), mu is the mean and sigma is the precision parameter and nu is the proportion of zeros of the zero adjusted Birnbaum-Saunders distribution.
#'
#'@references
#'Leiva, V., Santos-Neto, M., Cysneiros, F.J.A., Barros, M. (2015) A methodology for stochastic inventory models based on a zero-adjusted Birnbaum-Saunders distribution.
#'\emph{Applied Stochastic Models in Business and Industry.} \email{10.1002/asmb.2124}.
#'
#'@author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'@examples plotZARBS()
#'dat <- rZARBS(1000); hist(dat)
#'fit <- gamlss(dat~1,family=ZARBS(),method=CG())
#'meanZABS(fit)
#'
#'data(papatoes);
#'fit = gamlss(I(Demand/10000)~1,sigma.formula=~1, nu.formula=~1, family=ZARBS(mu.link="identity",sigma.link = "identity",nu.link = "identity"),method=CG(),data=papatoes)
#'summary(fit)
#'
#'data(oil)
#'fit1 = gamlss(Demand~1,sigma.formula=~1, nu.formula=~1, family=ZARBS(mu.link="identity",sigma.link = "identity",nu.link = "identity"),method=CG(),data=oil)
#'summary(fit1)
#'
#'@export
#'

qZMPL <- function (p, theta = 5, p0 = 0, lower.tail = TRUE, log.p = FALSE)
{
  if (log.p == TRUE)
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE)
    p <- p
  else p <- 1 - p
  if (any(p <= 0) | any(p > 1))
    stop(paste("p must be between 0 and 1", "\n", ""))

  p0plind <- ((theta^2)*(theta+2))/((theta+1)^3)
  pnew0 <- ((p - p0)/(1 - p0))
  pnew <- ifelse(pnew0 > 0, pnew0, 0)
  q <- qpoislind(pnew, theta = theta)

  q
}

#'Zero Modified Poisson-Lindley (ZMPL) distribution
#'
#'@description The fuction \code{ZMPL()} defines the ZMPL distribution, a two paramenter
#'distribution. The zero modified Poisson-Lindley distribution is similar
#'to the Poisson-Lindley distribution but allows zeros as y values. The extra parameter
#'models the probabilities at zero. The functions dZMPL, pZMPL, qZMPL and rZMPL define
#'the density, distribution function, quantile function and random generation for
#'the zero modified Poisson-Lindley distribution.
#'plotZMPL can be used to plot the distribution. meanZMPL calculates the expected
#'value of the response for a fitted model.
#'
#'@usage
#'dZMPL(x, mu = 1, sigma = 1, nu=0.1, log = FALSE)
#'pZMPL(q, mu = 1, sigma = 1,  nu=0.1, lower.tail = TRUE, log.p = FALSE)
#'qZMPL(p, mu = 1, sigma = 1,  nu=0.1, lower.tail = TRUE, log.p = FALSE)
#'rZMPL(n, mu = 1, sigma = 1)
#'plotZMPL(mu = .5, sigma = 1,  nu=0.1, from = 0, to = 0.999, n = 101, ...)
#'meanZMPL(obj)
#'
#' @param x,q vector of observations/quantiles
#' @param mu vector of scale parameter values
#' @param sigma vector of shape parameter values
#' @param log, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param from where to start plotting the distribution from
#' @param to up to where to plot the distribution
#' @param ... other graphical parameters for plotting
#'
#'
#' @details  The probability massa function of the is given by
#'
#'\deqn{ f_{Y}(y;\mu,\delta,p) =\frac{[1-p]\sqrt{\delta+1}}{4\,y^{3/2}\,\sqrt{\pi\mu}}\left[y+\frac{\delta\mu}{\delta+1} \right]\exp\left(-\frac{\delta}{4}\left[\frac{y[\delta+1]}{\delta\mu}+\frac{\delta\mu}{y[\delta+1]}-2\right]\right) I_{(0, \infty)}(y)+  pI_{\{0\}}(y).}
#'
#'@return returns a object of the ZMPL distribution.
#'
#'@note For the function ZMPL(), mu is the mean and sigma is the precision parameter and nu is the proportion of zeros of the zero adjusted Birnbaum-Saunders distribution.
#'
#'@references
#'Leiva, V., Santos-Neto, M., Cysneiros, F.J.A., Barros, M. (2015) A methodology for stochastic inventory models based on a zero-adjusted Birnbaum-Saunders distribution.
#'\emph{Applied Stochastic Models in Business and Industry.} \email{10.1002/asmb.2124}.
#'
#'@author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'@examples plotZARBS()
#'dat <- rZARBS(1000); hist(dat)
#'fit <- gamlss(dat~1,family=ZARBS(),method=CG())
#'meanZABS(fit)
#'
#'data(papatoes);
#'fit = gamlss(I(Demand/10000)~1,sigma.formula=~1, nu.formula=~1, family=ZARBS(mu.link="identity",sigma.link = "identity",nu.link = "identity"),method=CG(),data=papatoes)
#'summary(fit)
#'
#'data(oil)
#'fit1 = gamlss(Demand~1,sigma.formula=~1, nu.formula=~1, family=ZARBS(mu.link="identity",sigma.link = "identity",nu.link = "identity"),method=CG(),data=oil)
#'summary(fit1)
#'
#'@export
#'
pZMPL <- function (q, theta = 5, p0=0, lower.tail = TRUE, log.p = FALSE)
{
  cdf <- rep(0, length(q))
  cdf <- ppoislind(q, theta = theta, lower.tail = TRUE, log.p = FALSE)
  cdf <- p0 + (1 - p0) * cdf

  if (lower.tail == TRUE)
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE)
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}

#'Zero Modified Poisson-Lindley (ZMPL) distribution
#'
#'@description The fuction \code{ZMPL()} defines the ZMPL distribution, a two paramenter
#'distribution. The zero modified Poisson-Lindley distribution is similar
#'to the Poisson-Lindley distribution but allows zeros as y values. The extra parameter
#'models the probabilities at zero. The functions dZMPL, pZMPL, qZMPL and rZMPL define
#'the density, distribution function, quantile function and random generation for
#'the zero modified Poisson-Lindley distribution.
#'plotZMPL can be used to plot the distribution. meanZMPL calculates the expected
#'value of the response for a fitted model.
#'
#'@usage
#'dZMPL(x, mu = 1, sigma = 1, nu=0.1, log = FALSE)
#'pZMPL(q, mu = 1, sigma = 1,  nu=0.1, lower.tail = TRUE, log.p = FALSE)
#'qZMPL(p, mu = 1, sigma = 1,  nu=0.1, lower.tail = TRUE, log.p = FALSE)
#'rZMPL(n, mu = 1, sigma = 1)
#'plotZMPL(mu = .5, sigma = 1,  nu=0.1, from = 0, to = 0.999, n = 101, ...)
#'meanZMPL(obj)
#'
#' @param x,q vector of observations/quantiles
#' @param mu vector of scale parameter values
#' @param sigma vector of shape parameter values
#' @param log, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param from where to start plotting the distribution from
#' @param to up to where to plot the distribution
#' @param ... other graphical parameters for plotting
#'
#'
#' @details  The probability massa function of the is given by
#'
#'\deqn{ f_{Y}(y;\mu,\delta,p) =\frac{[1-p]\sqrt{\delta+1}}{4\,y^{3/2}\,\sqrt{\pi\mu}}\left[y+\frac{\delta\mu}{\delta+1} \right]\exp\left(-\frac{\delta}{4}\left[\frac{y[\delta+1]}{\delta\mu}+\frac{\delta\mu}{y[\delta+1]}-2\right]\right) I_{(0, \infty)}(y)+  pI_{\{0\}}(y).}
#'
#'@return returns a object of the ZMPL distribution.
#'
#'@note For the function ZMPL(), mu is the mean and sigma is the precision parameter and nu is the proportion of zeros of the zero adjusted Birnbaum-Saunders distribution.
#'
#'@references
#'Leiva, V., Santos-Neto, M., Cysneiros, F.J.A., Barros, M. (2015) A methodology for stochastic inventory models based on a zero-adjusted Birnbaum-Saunders distribution.
#'\emph{Applied Stochastic Models in Business and Industry.} \email{10.1002/asmb.2124}.
#'
#'@author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'@examples plotZARBS()
#'dat <- rZARBS(1000); hist(dat)
#'fit <- gamlss(dat~1,family=ZARBS(),method=CG())
#'meanZABS(fit)
#'
#'data(papatoes);
#'fit = gamlss(I(Demand/10000)~1,sigma.formula=~1, nu.formula=~1, family=ZARBS(mu.link="identity",sigma.link = "identity",nu.link = "identity"),method=CG(),data=papatoes)
#'summary(fit)
#'
#'data(oil)
#'fit1 = gamlss(Demand~1,sigma.formula=~1, nu.formula=~1, family=ZARBS(mu.link="identity",sigma.link = "identity",nu.link = "identity"),method=CG(),data=oil)
#'summary(fit1)
#'
#'@export
#'

dZMPL <- function (x, theta,  p0 = 0, log = FALSE)
{
  ans <- NULL
  x1 <- min(x)
  deflat.limit <- - (   (theta^2)*(theta+2)/( (theta^2) +3*theta+1)     )
  if(p0 < deflat.limit) ans <- NaN
  if(p0 > 1) ans <- NaN

  if(p0 == deflat.limit)
  {
    if(x1 == 0) stop(paste("The value of x must be >= 1", "\n", ""))

    if(log == TRUE){
      ans <- log1p(-p0) + dpoislind(x,theta, log = TRUE)
    } else{
      ans <- (1 - p0)*dpoislind(x, theta)
    }
  }
  if(p0 >= deflat.limit )
  {

    if(log==TRUE){
      if(x1 == 0)
      {
        if(length(x) > 1 )
        {
          ans[1] <- log(p0 + (1-p0)*dpoislind(0,theta))
          ans[2:length(x)] <- log((1 - p0)) + dpoislind(x[x!=0],theta,log=TRUE)
        } else{
          ans <- log(p0 + (1-p0)*dpoislind(0,theta))
        }
      }
      else{
        ans <- log((1 - p0)) + dpoislind(x,theta,log=TRUE)
      }
    }
    else{

      if(x1 == 0)
      {
        if(length(x) > 1){
          ans[1] <- p0 + (1-p0)*dpoislind(0,theta)
          ans[2:length(x)] <- (1 - p0)*dpoislind(x[x!=0],theta)
        } else{
          ans <- p0 + (1-p0)*dpoislind(0,theta)
        }

      }
      else{
        ans <- (1 - p0)*dpoislind(x,theta)
      }

    }
  }

  ans
}


#'Zero Modified Poisson-Lindley (ZMPL) distribution
#'
#'@description The fuction \code{ZMPL()} defines the ZMPL distribution, a two paramenter
#'distribution. The zero modified Poisson-Lindley distribution is similar
#'to the Poisson-Lindley distribution but allows zeros as y values. The extra parameter
#'models the probabilities at zero. The functions dZMPL, pZMPL, qZMPL and rZMPL define
#'the density, distribution function, quantile function and random generation for
#'the zero modified Poisson-Lindley distribution.
#'plotZMPL can be used to plot the distribution. meanZMPL calculates the expected
#'value of the response for a fitted model.
#'
#'@usage
#'dZMPL(x, mu = 1, sigma = 1, nu=0.1, log = FALSE)
#'pZMPL(q, mu = 1, sigma = 1,  nu=0.1, lower.tail = TRUE, log.p = FALSE)
#'qZMPL(p, mu = 1, sigma = 1,  nu=0.1, lower.tail = TRUE, log.p = FALSE)
#'rZMPL(n, mu = 1, sigma = 1)
#'plotZMPL(mu = .5, sigma = 1,  nu=0.1, from = 0, to = 0.999, n = 101, ...)
#'meanZMPL(obj)
#'
#' @param x,q vector of observations/quantiles
#' @param mu vector of scale parameter values
#' @param sigma vector of shape parameter values
#' @param log, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param from where to start plotting the distribution from
#' @param to up to where to plot the distribution
#' @param ... other graphical parameters for plotting
#'
#'
#' @details  The probability massa function of the is given by
#'
#'\deqn{ f_{Y}(y;\mu,\delta,p) =\frac{[1-p]\sqrt{\delta+1}}{4\,y^{3/2}\,\sqrt{\pi\mu}}\left[y+\frac{\delta\mu}{\delta+1} \right]\exp\left(-\frac{\delta}{4}\left[\frac{y[\delta+1]}{\delta\mu}+\frac{\delta\mu}{y[\delta+1]}-2\right]\right) I_{(0, \infty)}(y)+  pI_{\{0\}}(y).}
#'
#'@return returns a object of the ZMPL distribution.
#'
#'@note For the function ZMPL(), mu is the mean and sigma is the precision parameter and nu is the proportion of zeros of the zero adjusted Birnbaum-Saunders distribution.
#'
#'@references
#'Leiva, V., Santos-Neto, M., Cysneiros, F.J.A., Barros, M. (2015) A methodology for stochastic inventory models based on a zero-adjusted Birnbaum-Saunders distribution.
#'\emph{Applied Stochastic Models in Business and Industry.} \email{10.1002/asmb.2124}.
#'
#'@author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'@examples plotZARBS()
#'dat <- rZARBS(1000); hist(dat)
#'fit <- gamlss(dat~1,family=ZARBS(),method=CG())
#'meanZABS(fit)
#'
#'data(papatoes);
#'fit = gamlss(I(Demand/10000)~1,sigma.formula=~1, nu.formula=~1, family=ZARBS(mu.link="identity",sigma.link = "identity",nu.link = "identity"),method=CG(),data=papatoes)
#'summary(fit)
#'
#'data(oil)
#'fit1 = gamlss(Demand~1,sigma.formula=~1, nu.formula=~1, family=ZARBS(mu.link="identity",sigma.link = "identity",nu.link = "identity"),method=CG(),data=oil)
#'summary(fit1)
#'
#'@export
#'
rZMPL <- function(n, theta, p0){

  deflat.limit <- - (   ((theta^2)*(theta+2)) /( (theta^2) +3*theta+1)     )

  if(p0 == deflat.limit){
    u <- runif(n)
    p <- (theta*(theta+2))/( (theta^2) +3*theta+1 )
    lambda <- NULL
    x <- NULL
    for(i in 1:n)
    {
      lambda[i] <- ifelse(u[i] <= p, rexp(1,theta), rgamma(1,shape = 2, rate=theta) )
      x[i] <- 1 + rpois(1,lambda[i] )
    }
    x
  }
  else if(deflat.limit < p0 & p0 < 0){
    x<- qZMPL(runif(n),theta=theta,p0=p0)
  }
  else if(0<= p0 & p0 <= 1){
    x <- ifelse(rbinom(n, 1, p0), 0, rpoislind(n, theta))
  }
  else{x <- NaN}
  x
}


