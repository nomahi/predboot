\name{pred.lasso2}
\alias{pred.lasso2}
\title{Two stage bootstrap confidence intervals for prediction accuracy measures of logistic regression model by lasso estimation}
\usage{
pred.lasso2(model, data, B=1000, alpha=0.05, Cores=detectCores())
}
\arguments{
\item{model}{A model formula for the logistic regression model}

\item{data}{Dataset. The outcome variable must be named as \code{Y} and the ID variable must be named as \code{ID}}

\item{B}{The number of bootstrap resampling (defalt: 1000)}

\item{alpha}{The confidence level for the confidence interval of AUC (1-\code{alpha}; default: 95\%)}

\item{Cores}{The number of cores to be used for parallel computations (if the user don't want to implement parallel computations, please set to 1).}
}
\value{
Two stage bootstrap confidence intervals for prediction accuracy measures of logistic regression model. Note that the point estimates are not provided. These outputs are calculated by \code{pred.lasso}.
\itemize{
\item \code{glm.output}: \code{glmnet} output for the logistic regression.
\item \code{N.obs}: Number of observations.
\item \code{N.boot}: Number of bootstrap resampling.
\item \code{C.Harrell_2BSCI}: Two stage bootstrap confidence interval by the Harrell's bias correction.
\item \code{C.0.632_2BSCI}: Two stage bootstrap confidence interval by the 0.632 method.
\item \code{C.0.632p_2BSCI}: Two stage bootstrap confidence interval by the 0.632+ method.
}
}
\author{Hisashi Noma <noma@ism.ac.jp>}
\references{
Noma, H., Shinozaki, T., Iba, K., Teramukai, S. and Furukawa, T. A. (2020). 
Confidence intervals of prediction accuracy measures for multivariable prediction models based on the bootstrap-based optimism correction methods. 
arXiv:2005.01457.
\url{https://arxiv.org/abs/2005.01457}
}
\description{
Two stage bootstrap confidence intervals for prediction accuracy measures of logistic regression model by lasso estimation. 
}
\examples{
data(exdata)

pred.lasso2(Y ~ A65 + SEX + DIA + HYP + HRT + HIG + SHO + TTR, data=exdata, B=1000)
pred.lasso2(Y ~ A65 + SEX + DIA + HYP + HRT + HIG + SHO + TTR + PMI + HEI + WEI + HTN + SMK1 + SMK2 + LIP + PAN + FAM + ST4, data=exdata, B=1000)
}
