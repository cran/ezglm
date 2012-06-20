\name{ezglm}
\alias{ezglm}
\title{selects significant non-additive interaction between two variables using fast GLM implementation.}
\description{selects significant non-additive interaction between two variables using fast GLM implementation.}
\usage{
ezglm(y, x1, x2, thr = 1, family=c("gaussian","binomial"))
}
\arguments{
		\item{y}{response variable, of length \code{n}. This argument should be quantitative for least squares, and a two-level factor for logistic regression.}
		\item{x1}{the first predictor, of length \code{n}.}
		\item{x2}{the second predictor, of length \code{n}.}
		\item{thr}{p-value tolerance. Truncate any p-value to 1 if it is larger than \code{thr}. Defaults value is \code{1}.}
		\item{family}{a character string specifying the model to use, valid options are:
		\itemize{
		\item \code{"gaussian"} least squares regression (regression),
		\item \code{"binomial"} logistic regression (classification).
		}Default is \code{"gaussian"}.}
}

\details{
Motivated by pairwise gene interaction selection in genome-wide association study (GWAS),
this package implements fast and simplified least squares,
and logistic regression for efficiently selecting
the significant non-additive interactions between two variables.
Once a user specifies a reponse variable \code{y} and predictors \code{x1} and \code{x2},
for \code{"gaussian"} a least squares model \code{y ~ x1 + x2 + x1*x2} is fitted;
for \code{"binomial"} a logistic regression \code{logit ~ x1 + x2 + x1*x2} is fitted.
Users can then select significant x1*x2 term using returned p-value of Wald test.
}


\value{
A matrix of coefficients.
}

\author{Yi Yang\cr
Maintainer: Yi Yang  <yiyang@umn.edu>}

\examples{
n = 10000
x1 = rnorm(n)
x2 = rnorm(n)
y1 = sample(c(0,1),n,rep=TRUE)
y2 = rnorm(n)

system.time(m1 <- ezglm(y1, x1, x2, 1, family = "binomial"))
m1

system.time(m2 <- glm(y1~x1+x2+x1*x2, family = binomial))
summary(m2)$coef

system.time(m3 <- ezglm(y2, x1, x2, 1, family = "gaussian"))
m3

system.time(m4 <- glm(y2~x1+x2+x1*x2, family = gaussian))
summary(m4)$coef

}
\keyword{models}
\keyword{regression}
