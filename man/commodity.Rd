\name{commodity}
\alias{dataset}
\docType{data}
\title{Data: Monthly growth rates of commodity prices, exchange rate and industrial production index.}
\description{Monthly growth rates of commodity prices, exchange rate and industrial production index from February 1980 until October 2010. Levels of these
series can be downloaded from IMF and Federal Reserve Bank of St. Louis.}
\usage{data("commodity")}
\format{
  A data frame with 441 observations on the following 8 variables.
  \describe{
    \item{\code{X_date_}}{a vector with dates}
    \item{\code{dlnbev}}{a numeric vector}
    \item{\code{dlnind}}{a numeric vector}
    \item{\code{dlnrawm}}{a numeric vector}
    \item{\code{dlnmeta}}{a numeric vector}
    \item{\code{dlnoil}}{a numeric vector}
    \item{\code{dlnipi}}{a numeric vector}
    \item{\code{dlnex}}{a numeric vector}
  }
}
\source{IMF Primary Commodity Prices (http://www.imf.org/external/np/res/commod/index.aspx) and Federal Reserve Bank of St. Louis (https://fred.stlouisfed.org).}
\examples{data(dataset)}
\keyword{datasets}
