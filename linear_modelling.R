#' Interpret glm(er) with a binary response
#' @title Interpret glm(er) with a binary response
#' @param model A model produced by \code{\link[stats]{glm}} or \code{\link[lmer4]{glmer}}, where the response variable is binary.
#' @return Prints one sentence for each predictor variable in the form: "for each increase of 1 unit \code{predictor},
#'     \code{response} increases \code{x} times (p=\code{p value})"
#' @note The response variable must be binary to be accurate, as the function reports the exp of the log odds.
#'
#' @examples
#' interp.glm.bin(model)
#'
#' @export
interp.glm.bin<- function(model){
  for (i in 2:length(summary(model)$coefficients[,1])){
    print(paste0("for each increase of 1 unit ", rownames(summary(model)$coefficients)[i], ", ",
                 strsplit(as.character(summary(model)$call[2])," ~")[[1]][1]," increases ",
                 round(exp(summary(model)$coefficients[i,"Estimate"]),digits = 4),
                 "times (p=",round(summary(model)$coefficients[i,"Pr(>|z|)"],digits=4),")"))
  }}

#' Interpret lm(er)
#' @title Interpret lm(er)
#' @param model A model produced by \code{\link[stats]{lm}} or \code{\link[lmer4]{lmer}}.
#' @return Prints one sentence for each predictor variable in the form: "for each increase of 1 unit (predictor),
#'     response increases \code{x} units (p=(p value))"
#'
#' @examples
#' interp.lm(model)
#'
#' @export
interp.lm <- function(model){
  for (i in 2:length(summary(model)$coefficients[,1])){
    print(paste0("for each increase of 1 unit ", rownames(summary(model)$coefficients)[i],", ",
                 strsplit(as.character(summary(model)$call[2])," ~")[[1]][1]," increases ",
                 round(summary(model)$coefficients[i,"Estimate"],digits = 4),
                 " units (p=",round(summary(model)$coefficients[i,"Pr(>|t|)"],digits=4),")"))
  }}


