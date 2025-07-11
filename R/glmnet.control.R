#' internal glmnet parameters
#'
#' View and/or change the factory default parameters in glmnet
#'
#' If called with no arguments, \code{glmnet.control()} returns a list with the
#' current settings of these parameters. Any arguments included in the call
#' sets those parameters to the new values, and then silently returns. The
#' values set are persistent for the duration of the R session.
#'
#' @param fdev minimum fractional change in deviance for stopping path; factory
#' default = 1.0e-5
#' @param devmax maximum fraction of explained deviance for stopping path;
#' factory default = 0.999
#' @param eps minimum value of lambda.min.ratio (see glmnet); factory default=
#' 1.0e-6
#' @param big large floating point number; factory default = 9.9e35. Inf in
#' definition of upper.limit is set to big
#' @param mnlam minimum number of path points (lambda values) allowed; factory
#' default = 5
#' @param pmin minimum probability for any class. factory default = 1.0e-9.
#' Note that this implies a pmax of 1-pmin.
#' @param exmx maximum allowed exponent. factory default = 250.0
#' @param prec convergence threshold for multi response bounds adjustment
#' solution. factory default = 1.0e-10
#' @param mxit maximum iterations for multiresponse bounds adjustment solution.
#' factory default = 100
#' @param itrace If 1 then progress bar is displayed when running \code{glmnet}
#' and \code{cv.glmnet}. factory default = 0
#' @param epsnr convergence threshold for \code{glmnet.fit}. factory default =
#' 1.0e-8
#' @param mxitnr maximum iterations for the IRLS loop in \code{glmnet.fit}. factory
#' default = 25
#' @param factory If \code{TRUE}, reset all the parameters to the factory
#' default; default is \code{FALSE}
#' @return A list with named elements as in the argument list
#' @author Jerome Friedman, Kenneth Tay, Trevor Hastie\cr Maintainer: Trevor Hastie
#' \email{hastie@@stanford.edu}
#' @seealso \code{glmnet}
#' @keywords models regression
#' @examples
#'
#' glmnet.control(fdev = 0)  #continue along path even though not much changes
#' glmnet.control()  # view current settings
#' glmnet.control(factory = TRUE)  # reset all the parameters to their default
#'
#' @export glmnet.control
glmnet.control <-
  function (fdev = 1e-05, devmax = 0.999, eps = 1e-06, big = 9.9e+35,
            mnlam = 5, pmin = 1e-09, exmx = 250, prec = 1e-10, mxit = 100,
            itrace = 0, epsnr = 1e-08, mxitnr = 100, factory = FALSE)
{
  list(fdev = fdev, devmax = devmax,
                             eps = eps, big = big, mnlam = mnlam, pmin = pmin,
                             exmx = exmx, prec = prec, mxit = mxit, itrace = itrace,
                             epsnr = epsnr, mxitnr = mxitnr)
}
