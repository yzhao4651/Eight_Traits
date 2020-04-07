### this one help me to get which variables has more missing values and help to set up 
###the condition should larger then 900, or other number 
### in the function BLUP
#out_start <- 4
#out_end <- 25
#for (i in out_start:out_end){
#t <-cbind(colnames(normadata[i]),sum(is.na(normadata[normadata$Year==2017 | normadata$Year==2018,][i])))
#print(t)
#}
    
###function for the ranef of all of traits 
ranefvalue <- function(out_start,out_end,y){
  require(lme4)
  require(Matrix)
  Entry = as.factor(y$Entry)
  Rep = as.factor(y$Rep)
  Year = as.factor(y$Year)
  out_nvar=out_end-out_start+1
  out_variable = colnames(y[out_start:out_end])
  out_beta <- matrix(NA_real_, nrow = length(levels(y$Entry)),
                     ncol = out_nvar,
                     dimnames = list(levels(y$Entry), 
                                     out_variable))

  for (i in out_start:out_end){

      outcome = colnames(y)[i]
      model <- lmer(get(outcome)~ (1|Entry)+ (1|Rep),
                    na.action = na.exclude,
                    data=y)
      beta <- ranef(model)

    out_beta[rownames(beta$Entry)] = beta$Entry[[1]]
  }
 return(out_beta)
}


