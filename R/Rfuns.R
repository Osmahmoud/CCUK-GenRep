require(tools)
require(knitr)
require(rmarkdown)
require(blme)
require(lme4)
require(ggplot2)
require(haven)
require(R2MLwiN)

# require(foreign)
# require(pander)
# require(knitr)
# require(grid)
# require(scales)
# require(gridExtra)
# require(PairedData)

#' Generate individual centre reports - in pdf - for the CCUK analysis using data in a specific format
genReport <- function(output.dir = paste(getwd(), "/Outputs/", sep = ""), Hub = "South West", code.hub = "SWSW", Overwrite = FALSE){
  
  Hub.Name = Hub
  codedHub.Name = code.hub
  output.file = paste("Report_", Hub, ".pdf", sep = "")
  
  ### render documentation: Generate pdf documenting variables based on defined settings
  if (output.file %in% list.files(output.dir) & !(Overwrite)){
    warning("There is already a file with the same name in this location.
The output file name is modified by adding the current date to the end of the name.
\n If you wish to overwrite the existing file, please set the 'Overwrite' to TRUE.",
            call. = FALSE, immediate. = TRUE)
    
    output.file = paste(tools::file_path_sans_ext(output.file),
                        "_", format(Sys.time(), format = "%d%m%y"),
                        ".pdf", sep = "")
  }
  
  render(paste(getwd(), "/Rmd/Master.Rmd", sep = ""), "pdf_document", output_dir = output.dir, output_file = output.file,
         quiet = TRUE)
  cat("Finished - Report for ", Hub, " has been successfully generated ...\n", sep = "")
}

#' VPC estimation - using Simulation method
#' @description The Variance Partition Coefficient of a multilevel model using 'simulation'.
#' This method is described by the Bristol centre for multilevel modelling. You can find it here:
#' http://www.bristol.ac.uk/cmm/software/support/support-faqs/pval.html
#' @param model the fitted model
#' @param X.values matrix of specific values of x by which you wish to estimate the VPC. Each column represents
#' a variable. Each row is a value. Note value for x0 (first column) should always be set as 1.
#' @param m number of simulated variance values.
#' @example
#' fit_1 <- glmer(dento ~ age + sex + (1 | hub), family = binomial("logit"), data = CCUK, nAGQ = 100)
#' X.values = matrix(nrow = 2, byrow = TRUE, data = c(1,0,1,1,0,0))
#' VPC.est(fit_1, X.values, m=5000)
VPC.est = function(model, X.values = matrix(nrow = 2, byrow = TRUE, data = c(1,0,1,1,0,0)), m = 5000){
  # extract sd of random effect (sigma_u)
  H.sigma = sqrt(as.numeric(summary(model)$varcor))
  # extract the estimated parameters
  beta = summary(model)$coefficients[,"Estimate"]
  # Generate 'm' values of level2 residuals from norma(0, (sigma_u))
  L2.Residual = rnorm(m, 0, H.sigma)
  
  # calculate probabilities by adding the XB (i.e., given X.values * beta) to the generated residuals
  XB = apply(sweep(X.values, MARGIN=2, beta, FUN=`*`), 1, sum)
  Prob = c()
  for(i in 1:length(XB)){
    Prob = c(Prob, XB[i] + L2.Residual)
  }
  Prob = exp(Prob) / (1 + exp(Prob))
  Var = Prob*(1-Prob)
  V1 = mean(Var)
  V2 = var(Prob)
  VPC = V2 / (V1 + V2)
  return(list(V1 = V1, V2 = V2, VPC = VPC))
}

#' calculate model statistics
#' @param ML.model multilevel model
#' @param single level model
model.Stat = function(ML.model, single.model){
  VPC = VPC.est(ML.model)
  Pvalue = 1-pchisq(as.numeric(-2*(logLik(single.model)-logLik(ML.model))) ,1)
  
  beta = summary(ML.model)$coefficients[1, "Estimate"]
  H.sigma = sqrt(as.numeric(summary(ML.model)$varcor))
  PI_1 = Trans(beta, 0)
  if (H.sigma >= .01){
    CI.lower_1 = Trans(beta, -1.96*H.sigma)
    CI.upper_1 = Trans(beta, 1.96*H.sigma)
  } else {
    CI.lower_1 = Trans(beta, -1.96*VPC$V2)
    CI.upper_1 = Trans(beta, 1.96*VPC$V2)
  }
  return(list(VPC=VPC, Pvalue=Pvalue, beta=beta, H.sigma=H.sigma, PI_1=PI_1,
              CI.lower_1=CI.lower_1, CI.upper_1=CI.upper_1))
}

#' Transform model parameter estimates from logit scale to proportions
#' @param Beta parameter estimate
#' U a value that can be added to beta (zero: to transform beta; 1.96*sd to transform upper limit of C.I)
Trans <- function(Beta, U){
  odds = exp(Beta+U)
  Pi = odds/(1+odds)
  return(Pi)
}

#' calculate model statistics
#' @param ML.model multilevel model
#' @param single.model single level model
#' @param DF dataframe that contain variables to use for plot
#' @param Response The name of the considered response variable
#' @param Plot logical if true, a plot will be generated.
Uni.Stat = function(ML.model, single.model, DF, Response, Plot=TRUE,
                    codedHub.Name, Hub.Name){
  VPC = VPC.est(ML.model)
  Pvalue = 1-pchisq(as.numeric(-2*(logLik(single.model)-logLik(ML.model))) ,1)
  
  beta = summary(ML.model)$coefficients[1, "Estimate"]
  H.sigma = sqrt(as.numeric(summary(ML.model)$varcor))
  PI_1 = Trans(beta, 0)
  if (H.sigma >= .01){
    CI.lower_1 = Trans(beta, -1.96*H.sigma)
    CI.upper_1 = Trans(beta, 1.96*H.sigma)
  } else {
    CI.lower_1 = Trans(beta, -1.96*VPC$V2)
    CI.upper_1 = Trans(beta, 1.96*VPC$V2)
  }
  U0    <- ranef(ML.model, condVar = TRUE)
  U0se  <- sqrt(attr(U0[[1]], "postVar")[1, , ])
  U0tab <- cbind(U0[[1]], U0se)
  colnames(U0tab)[1] <- "U0"
  U0tab <- cbind(U0tab, Trans(beta, U0tab$U0))
  colnames(U0tab)[3] <- "PIj"
  U0tab <- U0tab[order(U0tab$PIj),]
  U0tab <- cbind(U0tab, c(1:dim(U0tab)[1]))
  colnames(U0tab)[4] <- c("U0_rank")
  #### Use This to report individual differences from mean
  U0tab_1 <- cbind(PIj = U0tab$PIj, lower = Trans(beta, U0tab$U0-1.96*U0tab$U0se), upper = Trans(beta, U0tab$U0+1.96*U0tab$U0se)) -  PI_1
  row.names(U0tab_1) <- row.names(U0tab)
  Hub.Size <- as.matrix(table(DF[['unit']][!(is.na(DF[[Response]]))]))
  U0tab_1 <- cbind(Size = Hub.Size[match(rownames(U0tab_1), rownames(Hub.Size))], U0tab_1)
  
  if(Plot){
    ## plot
    cutoff <- data.frame(x = c(-Inf, Inf), y = c(PI_1,PI_1))
    Ident <- which(rownames(U0tab_1) == codedHub.Name)
    
    P <- ggplot(U0tab, aes(U0_rank, PIj)) + geom_point(size=3, shape=18) +
      scale_x_continuous(breaks=Ident, labels=Hub.Name) +
      geom_segment(aes(x=U0tab$U0_rank, y=Trans(beta, U0tab$U0-1.96*U0tab$U0se), xend=U0tab$U0_rank, yend=Trans(beta, U0tab$U0+1.96*U0tab$U0se)),
                   size=0.6, arrow=arrow(angle=90, ends="both", length = unit(0.08, "inches"))) +
      geom_line(aes(x,y), cutoff, size=0.5, linetype="dashed") + scale_y_continuous(limits = c(0, 1)) + labs(x = "Centre", y = "Proportion")
    
    P <- P + theme(axis.text=element_text(size=14, face="bold"),
                   axis.text.x=element_text(angle=20, vjust=0.5, hjust = 0.5),
                   axis.text.y=element_text(vjust=0.5, hjust = 0.5),
                   plot.margin = unit(c(0.3,0.9,0.3,0.3), "cm"),
                   axis.title=element_text(size = rel(1.8)),
                   axis.title.x=element_text(vjust=-0.5),
                   axis.title.y=element_text(vjust=1.5),
                   panel.background = element_rect(fill = 'white', colour = 'black'),
                   panel.grid.major = element_line(colour = 'grey89'), panel.grid.minor = element_line(colour = 'grey89'))
    return(list(Stat = list(VPC=VPC, Pvalue=Pvalue, beta=beta, U0tab=U0tab_1, PI_1=PI_1, CI.lower_1=CI.lower_1, CI.upper_1=CI.upper_1), Plot=P))
  } else {
    return(list(VPC=VPC, Pvalue=Pvalue, beta=beta, U0tab=U0tab_1, PI_1=PI_1, CI.lower_1=CI.lower_1, CI.upper_1=CI.upper_1))
  }
}

#' Table's caption for latex format
#' @description return a latex format to include table caption
#' @param caption string of your caption
tabCap <- function(caption){
  paste0("\\textit{", caption, "}")
}

#' deterministic decimal places
#' @description specify how many decimal numbers do you want to return
#' @param x the given numer (or a vector of given numbers)
#' @param k number of decimal places required
Specify.dec <- function(x, k = 2){
  format(round(x, k), nsmall=k)
}

#' Show plus sign for positive numbers
#' @description return a text showing the plus '+' sign in front of positive numbers.
Sign.dif <- function(Number){
  if(as.numeric(Number) > 0){
    paste0("+",Number)
  } else
    Number
}

#' Pre-processing data imported by the 'haven' package
#' @description This function is used to present factors by labels, replace NaN by NA (in order to have all missing values coded in the same way)
#' @param df data frame imported from Stata data file by the rio package
#' @export
Post.Haven <- function(df){
  DF = as.data.frame(df)
  DF[is.na(DF)] <- NA
  f = function(x){
    if(haven::is.labelled(x)){Out = haven::as_factor(x)} else {Out = x}
  }
  for(i in 1:dim(DF)[2]){
    DF[,i] <- f(DF[,i])
  }
  return(DF)
}

#' Categorize SDQ measures
#' @param x variable (an SDQ measure) to be categorize
#' @param cut.breaks The cut points for splitting categories
#' @param levels labels of the generated levels
Categ <- function(x, cut.breaks, levels = c("Avg.", "S.raised", "high", "V.high")){
  factor(cut(x, breaks=cut.breaks), labels = levels)
}

#' Extract DIC info from MLwiN Outout
#' @description return a data frame containing information of the DIC of the model fitted using R2MLwiN package
#' @param MLwiN.model The model outcome obtained by the 'R2MLwiN' package (e.g. using \code{sink})
get.DIC <- function(MLwiN.model){
  sink("MLwiN_Out.txt")
  print(summary(MLwiN.model))
  sink()
  MLwiN.Outputs <- readLines("MLwiN_Out.txt")
  unlink("MLwiN_Out.txt")
  LnId <- which(!is.na(stringr::str_match(MLwiN.Outputs, "Bayesian Deviance Information Criterion \\(DIC\\)")))
  List <- strsplit(MLwiN.Outputs[c(LnId+1, LnId+2)], "\\s+")
  Re <- as.numeric(List[[2]]); names(Re) <- List[[1]]
  Re
}

#' Extract fixed effect estimates from MLwiN Outout
#' @description return a data frame containing fixed effect estimates of the model fitted using R2MLwiN package - MIGHT NEED REVISE AS IT MIGHT GET WARNING MESSAGE IF THE CELL (***) IS EMPTY!!!
#' @param MLwiN.model The model outcome obtained by the 'R2MLwiN' package
Get.Fixed.est <- function(MLwiN.model){
  sink("MLwiN_Out.txt")
  print(summary(MLwiN.model))
  sink()
  MLwiN.Outputs <- readLines("MLwiN_Out.txt")
  unlink("MLwiN_Out.txt")
  LnId.start <- which(!is.na(stringr::str_match(MLwiN.Outputs, "The fixed part estimates:")))
  LnId.end <- which(!is.na(stringr::str_match(MLwiN.Outputs, "Signif. codes:\\s*\\.*")))
  # split items by spaces
  List <- strsplit(MLwiN.Outputs[(LnId.start+2):(LnId.end-1)], "\\s+")
  Raw.Re <- data.frame(do.call("rbind", List))
  names(Raw.Re) <- c("", "Coef.", "Std. Err.", "z", "Pr(>|z|)", "", "95% CI.lower", "95% CI.upper", "ESS")
  Raw.Re[,-c(1,6)] <- apply(Raw.Re[,-c(1,6)], 2, function(x){as.numeric(as.character(x))})
  Re <- Raw.Re
  return(Re)
}

#' Extract random part estimates from MLwiN Outout
#' @description return a vector containing random effect estimates of the model fitted using R2MLwiN package
#' @param MLwiN.model The model outcome obtained by the 'R2MLwiN' package
Get.Random.est <- function(MLwiN.model){
  sink("MLwiN_Out.txt")
  print(summary(MLwiN.model))
  sink()
  MLwiN.Outputs <- readLines("MLwiN_Out.txt")
  unlink("MLwiN_Out.txt")
  LnId <- which(!is.na(stringr::str_match(MLwiN.Outputs, "^var_\\w+")))
  List <- strsplit(MLwiN.Outputs[LnId], "\\s+")
  Raw.Re <- data.frame(do.call("rbind", List))
  names(Raw.Re) <- c("", " Coef.", "Std. Err.", "95% CI.lower", "95% CI.upper", "ESS")
  Raw.Re[,-1] <- apply(Raw.Re[,-1], 2, function(x){as.numeric(as.character(x))})
  Re <- Raw.Re
  return(Re)
}

#' Collect statistics for summary table
#' @description Organise all required information for summary table and return them all in a list
#' 
Set.Sumary = function(Model, P, CI.l, CI.u, VPC, Pval, tab.ind){
  n = summary(Model)$devcomp$dims["n"]
  return(list(n=n, P=P, CI.l=CI.l, CI.u=CI.u, VPC=VPC, Pval=Pval, tab=tab.ind))
}

SDQ.Sumary = function(N, P, CI, VPC, Pval, Size, tab.ind, Cat = "Avg.", codedHub.Name){
  n = sum(N)
  P = P[Cat]
  CI.l = CI[Cat, "lower"]
  CI.u = CI[Cat, "upper"]
  Size = Size[codedHub.Name]
  P.ind = tab.ind[Cat,1,codedHub.Name] - P
  l.ind = tab.ind[Cat,2,codedHub.Name] - P
  u.ind = tab.ind[Cat,3,codedHub.Name] - P
  return(list(n=n, P=P, CI.l=CI.l, CI.u=CI.u, VPC=VPC, Pval=Pval, Size=Size, P.ind=P.ind, l.ind=l.ind, u.ind=u.ind))
}

THEME = theme(axis.text=element_text(size=12, face="bold"),
              axis.text.x=element_text(angle=20, vjust=0.5, hjust = 0.5),
              axis.text.y=element_text(vjust=0.5, hjust = 0.5),
              plot.margin = unit(c(0.3,0.9,0.3,0.3), "cm"),
              axis.title=element_text(size = rel(1.4)),
              axis.title.x=element_text(vjust=-0.5),
              axis.title.y=element_text(vjust=1.5),
              panel.background = element_rect(fill = 'white', colour = 'black'),
              panel.grid.major = element_line(colour = 'grey89'), panel.grid.minor = element_line(colour = 'grey89'))