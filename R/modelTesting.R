#' Perform Two-Sample t-tests using `limma` package in R 
#'
#' Fit linear model for each gene given a series of arrays using `limma` package in R.
#' 
#' @param data input a expression data.frame (probes in row; samples in column).
#' @param ctlS input the sample names of the control group as the reference. 
#' @param trtS input the sample names of the treatment group.
#' @param ctlS.name input the name of the control group using in the comparison. 
#' @param trtS.name input the name of the treatment group using in the comparison.
#'  
#' @return a list of two tables. 'dds' contains the original data used to perform the test;
#'  `results` contains the testing statistics. 
#'   
#' @seealso \code{\link{lmFit}} which this function wraps.
#' @export
#' @examples
#' sd <- 0.3*sqrt(4/rchisq(100,df=4))
#' y <- data.frame(matrix(rnorm(100*6,sd=sd),100,6))
#' rownames(y) <- paste("Probe",1:100, sep="_")
#' colnames(y) <- c(paste("A", 1:3, sep=""), paste("B", 1:3, sep=""))
#' y[1:50,4:6] <- y[1:50,4:6] + 2
#' y[sample(seq(1,100), 10), 1] <- NA
#' y[sample(seq(1,100), 10), 3] <- NA
#' y[sample(seq(1,100), 10), 5] <- NA
#' Controls=c("A1","A2","A3")
#' Tumors=c("B1","B2","B3")
#' comTest <- limma_ttest(y, ctlS=Controls, trtS=Tumors)
#' comTest2 <- limma_ttest(y, ctlS=Controls, trtS=Tumors, ctlS.name="Normal", trtS.name="Cancer")
#' names(comTest)
#' comTest[[2]][1:5,]
#' comTest2[[2]][1:5,]
limma_ttest <- function(data, ctlS, trtS, ctlS.name=NULL, trtS.name=NULL){
  dat.ctl <- data[,ctlS]
  dat.trt <- data[,trtS]
  dat.final <- cbind(data[,ctlS], data[,trtS])
  if ((is.null(ctlS.name) & !is.null(trtS.name)) | (!is.null(ctlS.name) & is.null(trtS.name))){
    stop("ctlS.name and trtS.name must coexist!!!")
  }
  if (is.null(ctlS.name) & is.null(trtS.name)){
    dat.target <- c(rep(as.character(substitute(ctlS)), length(ctlS)),
                    rep(as.character(substitute(trtS)), length(trtS)))
    Group <- factor(dat.target, levels=c(substitute(ctlS),substitute(trtS)))
    design <- model.matrix(~Group)
    colnames(design) <- c(substitute(ctlS), 
                          paste(substitute(trtS), "vs", substitute(ctlS), sep=""))
  }
  if (!is.null(ctlS.name) & !is.null(trtS.name)){
    dat.target <- c(rep(as.character(ctlS.name), length(ctlS)),
                    rep(as.character(trtS.name), length(trtS)))
    Group <- factor(dat.target, levels=c(ctlS.name,trtS.name))
    design <- model.matrix(~Group)
    colnames(design) <- c(ctlS.name, 
                          paste(trtS.name, "vs", ctlS.name, sep=""))
  }
  
  fit <- limma::lmFit(dat.final, design)
  fit <- limma::eBayes(fit)
  #print(design)
  #return(fit)
  
  ### save result ###
  ### fold change ###
  group.mean2 <- apply(dat.final,1, function (i) {
    tapply(i, Group, mean, na.rm=TRUE)
  })
  
  group.mean=t(group.mean2)
  colnames(group.mean)
  diff.group <- group.mean[,2] - group.mean[,1]
  foldchange <- sign(diff.group)*2^(abs(diff.group))
  
  tt.res <- matrix(0, nrow(data), 3)
  colnames(tt.res) <- c("p.value", "adj.p.val", "t.stat")
  cN <- colnames(design)[2]
  tt.res[,1] <- fit$p.value[,cN]
  tt.res[,2] <- limma::topTable(fit, coef=cN, adjust="BH", number=nrow(data))[rownames(data),"adj.P.Val"]
  tt.res[,3] <- fit$t[,cN]
  rank.res <- rank(tt.res[,1])
  tt <- data.frame(group.mean, diff.group=diff.group,
                   foldchange=foldchange,  tt.res, Rank=rank.res)
  #return(tt)                   
  #######################
  res.lis <- list()
  res.lis[[1]] <- dat.final
  res.lis[[2]] <- tt
  names(res.lis) <- c("dds", "results")
  return(res.lis)
}





#' Perform One-way ANOVA with contrast using `limma` package in R 
#'
#' Fit linear model for each gene given a series of arrays using `limma` package in R.
#' 
#' @param data input a expression data.frame (probes in row; samples in column).
#' @param fac input the sample grouping info that you want to compare.
#' @param contrast inpt the contrast between any two given groups in 'fac'.
#' 
#' @return a list of two tables. 'dds' contains the original data used to perform the test;
#'  `results` contains the testing statistics. 
#'   
#' @seealso \code{\link{lmFit}} which this function wraps.
#' @export
#' @examples
#' sd <- 0.3*sqrt(4/rchisq(100,df=4))
#' y <- data.frame(matrix(rnorm(100*9,sd=sd),100,9))
#' rownames(y) <- paste("Probe",1:100, sep="_")
#' colnames(y) <- c(paste("A", 1:3, sep=""), paste("B", 1:3, sep=""), paste("C", 1:3, sep=""))
#' y[1:50,4:6] <- y[1:50,4:6] + 2
#' y[51:75,7:9] <- y[51:75,7:9] + 3
#' y[sample(seq(1,100), 10), 1] <- NA
#' y[sample(seq(1,100), 10), 3] <- NA
#' y[sample(seq(1,100), 10), 5] <- NA
#' mFAC <- c(rep("Control", 3), rep("Treat1", 3), rep("Treat2", 3))
#' mContr <- c("Treat1-Control", "Treat2-Control", "Treat2-Treat1")
#' comTest <- limma_anova_contrast(data=y, fac=mFAC, contrast=mContr)
#' names(comTest)
#' names(comTest[["results"]])
#' mContr2 <- c("Treat1-Control", "Treat2-Control")
#' comTest2 <- limma_anova_contrast(data=y, fac=mFAC, contrast=mContr2)
#' names(comTest2)
#' names(comTest2[["results"]])
limma_anova_contrast <- function(data, fac, contrast){
  if (is.factor(fac)) {
    fac <- fac
  } else {
    fac <- factor(fac)
  }
  ### design matrix ###
  design <- model.matrix(~0+fac)
  colnames(design) <- levels(fac)
  
  ### pair-wise comparisons ###
  fit <- limma::lmFit(data, design)
  contrast.matrix <- limma::makeContrasts(contrasts=contrast,
                                          levels=design)
  fit2 <- limma::contrasts.fit(fit, contrast.matrix)
  fit2 <- limma::eBayes(fit2)
  
  ### save result ###
  ### fold change ###
  group.mean2 <- apply(data,1, function (i) {
    tapply(i, fac, mean, na.rm=TRUE)
  })
  
  group.mean=t(group.mean2)
  colnames(group.mean)
  
  fcList <- list()
  for (i in 1:length(contrast)){
    ind <- unlist(strsplit(contrast[i], "-"))
    diff.group <- group.mean[,ind[1]] - group.mean[,ind[2]]
    foldchange <- sign(diff.group)*2^(abs(diff.group))
    tt <- data.frame(group.mean[,ind[1]],
                     group.mean[,ind[2]],
                     diff.group,
                     FoldChange=foldchange,
                     p.value=fit2$p.value[,contrast[i]],
                     t.stat=fit2$t[,contrast[i]],
                     Rank=rank(fit2$p.value[,contrast[i]]))
    colnames(tt) <- c(ind[1], ind[2], "diff.group", "FoldChange", "p.value", "t.stat", "Rank")
    fcList[[i]] <- tt
  }
  over.mod <- limma::topTable(fit2, adjust="BH", number=nrow(data))[rownames(data),]
  colnames(over.mod)[1:length(contrast)] <- paste("diff.", colnames(over.mod)[1:length(contrast)], sep="")
  over.mod <- data.frame(group.mean, over.mod, "Rank"=rank(over.mod$P.Value))
  colnames(over.mod) <- gsub("P.V", "p.v", colnames(over.mod))
  fcList[[length(contrast)+1]] <- over.mod
  names(fcList) <- c(contrast, "overall")

  #######################
  res.lis <- list()
  res.lis[[1]] <- data
  res.lis[[2]] <- fcList
  names(res.lis) <- c("dds", "results")
  return(res.lis)  
}




#' Perform One-way ANOVA with TukeyHSD using `aov` package in R 
#'
#' Fit linear model for each gene given a series of arrays using `aov` package in R.
#' 
#' @param expr input a vector of expression values.
#' @param fac input the sample grouping info that you want to compare.
#' @param include input a subset of data if needed. Default is to use all data.
#' @param showname default is NULL. It will show the name of the method selected,
#'  otherwise the character would be used to assign the rowname of the result.
#'   
#' @return a list of two tables. 'dds' contains the original data used to perform the test;
#'  `results` contains the testing statistics.  
#'  
#' @seealso \code{\link{aov}} and \code{\link{TukeyHSD}} which this function wraps.
#' @export
#' @examples
#' sd <- 0.3*sqrt(4/rchisq(100,df=4))
#' y <- data.frame(matrix(rnorm(100*12,sd=sd),100,12))
#' rownames(y) <- paste("Probe",1:100, sep="_")
#' colnames(y) <- c(paste("A", 1:3, sep=""), paste("B", 1:3, sep=""),
#'  paste("C", 1:3, sep=""), paste("D", 1:3, sep=""))
#' y[1:50,4:6] <- y[1:50,4:6] + 2
#' y[51:75,7:9] <- y[51:75,7:9] + 3
#' y[sample(seq(1,100), 10), 1] <- NA
#' y[sample(seq(1,100), 10), 3] <- NA
#' y[sample(seq(1,100), 10), 5] <- NA
#' mFAC <- c(rep("Control", 3), rep("Treat1", 3), rep("Treat2", 3) , rep("Treat3", 3))
#' aovTukeyHSD_Expr(expr=y[1,], fac=mFAC, showname="expr")
#' aovTukeyHSD_Expr(expr=y[1,], fac=mFAC, include=c("Control", "Treat1", "Treat2"), showname="expr")
aovTukeyHSD_Expr <- function(expr, fac, include=NULL, showname=NULL){
  if (is.factor(fac)) {
    fac <- fac
  } else {
    fac <- factor(fac)
  }
  
  if (is.null(include)==FALSE){
    expr <- as.numeric(expr[which(fac %in% include)])
    fac <- factor(fac[which(fac %in% include)], levels=include)
  } 
  aovFit <- try(aov(as.numeric(expr)~fac, na.action=na.omit), silent=TRUE)
  if(class(aovFit)[1]!='try-error' & length(unique(na.omit(as.numeric(expr))))>2) {
    aovP <- summary(aovFit)[[1]][1, 5] # ANOVA P
    names(aovP) <- 'p.value'
    tukey <- TukeyHSD(aovFit)[[1]]
    dif <- tukey[, 'diff']
    names(dif) <- paste('diff', names(dif), sep=".")
    FC <- sign(dif)*2^(abs(dif))
    names(FC) <- paste('FC', names(FC), sep=".")
    padj <- tukey[, 'p adj']
    names(padj) <- paste('padj', names(padj), sep=".")
  } else {
    aovP <- NA
    dif <- rep(NA, length=length(levels(fac)))
    padj <- rep(NA, length=length(levels(fac)))
    FC <- rep(NA, length=length(levels(fac)))
  }
  res <- t(data.frame(c(aovP, padj, dif, FC)))
  iName <- ifelse(is.null(showname), substitute(expr), showname)
  rownames(res) <- as.character(iName)
  return (res)
}





#' Perform One-way ANOVA with TukeyHSD using `aov` package in R 
#'
#' Fit linear model for each gene given a series of arrays using `aov` package in R.
#' 
#' @param mat input a matrix of expression values.
#' @param fac input the sample grouping info that you want to compare.
#' @param include input a subset of data if needed. Default is to use all data.
#' 
#' @return a list of two tables. 'dds' contains the original data used to perform the test;
#'  `results` contains the testing statistics.  
#'  
#' @seealso \code{\link{aov}} and \code{\link{TukeyHSD}} which this function wraps.
#' @export
#' @examples
#' sd <- 0.3*sqrt(4/rchisq(100,df=4))
#' y <- data.frame(matrix(rnorm(100*12,sd=sd),100,12))
#' rownames(y) <- paste("Probe",1:100, sep="_")
#' colnames(y) <- c(paste("A", 1:3, sep=""), paste("B", 1:3, sep=""),
#'  paste("C", 1:3, sep=""), paste("D", 1:3, sep=""))
#' y[1:50,4:6] <- y[1:50,4:6] + 2
#' y[51:75,7:9] <- y[51:75,7:9] + 3
#' y[sample(seq(1,100), 10), 1] <- NA
#' y[sample(seq(1,100), 10), 3] <- NA
#' y[sample(seq(1,100), 10), 5] <- NA
#' mFAC <- c(rep("Control", 3), rep("Treat1", 3), rep("Treat2", 3) , rep("Treat3", 3))
#' comTest <- aovTukeyHSD_VecMat(mat=y, fac=mFAC)
#' comTest2 <- aovTukeyHSD_VecMat(mat=y, fac=mFAC, include=c("Control", "Treat1", "Treat2"))
aovTukeyHSD_VecMat <- function(mat, fac, include=NULL){
  if(ncol(mat)!=length(fac)) stop(sprintf("mat ncol: %d\nfac length: %d\n", ncol(mat), length(fac)))
  
  if (is.null(include)==FALSE){
    mat <- mat[,which(fac %in% include)]
    fac <- factor(fac[which(fac %in% include)], levels=include)
  }
  
  res <- NULL
  for (i1 in 1:nrow(mat)){
    res <- rbind(res, aovTukeyHSD_Expr(expr=mat[i1,],
                                       fac=fac,
                                       showname=rownames(mat)[i1]))
  }
  rank.res <- rank(as.numeric(res[,"p.value"]))
  res.out <- data.frame(res, Rank=rank.res)
  return(list(dds=mat, results=res.out))
} 





















#' Perform Correlation tests using `Hmisc` package in R 
#'
#' rcorr Computes a matrix of Pearson's r or Spearman's rho rank correlation
#'  coefficients for all possible pairs of columns of a matrix. Missing values
#'  are deleted in pairs rather than deleting all rows of x having any missing
#'  variables. Ranks are computed using efficient algorithms, using midranks for
#'  ties. 
#'  
#' @param x a numeric matrix with at least 5 rows and at least 2 columns
#'  (if y is absent). For print, x is an object produced by rcorr.
#' @param y a numeric vector or matrix which will be concatenated to x. If y is
#'  omitted for rcorr, x must be a matrix.
#' @param method specifies the type of correlations to compute. Spearman correlations
#'  are the Pearson linear correlations computed on the ranks of non-missing elements,
#'  using midranks for ties.
#' @param showname default is NULL. It will show the name of the method selected,
#'  otherwise the character would be used to assign the rowname of the result.
#'   
#' @return the correlation coefficient and p-value. 
#'  
#' @seealso \code{\link{rcorr}} which this function wraps.
#' @export
#' @examples
#' rcorr_XY(x=mtcars[,1], y=mtcars[,2], method='pearson')
#' rcorr_XY(x=mtcars[,1], y=mtcars[,2], method='spearman', showname="IC50")
rcorr_XY <- function(x, y, method=c("pearson", "spearman"), showname=NULL){
  df <- data.frame(x, y)
  df1 <- df[complete.cases(df), ]
  method <- match.arg(method)
  if (method=="pearson") {
    tp <- Hmisc::rcorr(x=df1$x, y=df1$y, type="pearson")
    res <- data.frame(r.coef=tp$r[1, 2], pearson.pval=tp$P[1, 2])
  }
  if (method=="spearman") {
    tp <- Hmisc::rcorr(x=df1$x, y=df1$y, type="spearman")
    res <- data.frame(rho.coef=tp$r[1, 2], spearman.pval=tp$P[1, 2])
  }
  iName <- ifelse(is.null(showname), as.character(method), showname)
  rownames(res) <- iName
  return(res)
}


#' Perform Correlation tests between each row in a matrix and a vector using `Hmisc` package in R 
#'
#' rcorr Computes a matrix of Pearson's r or Spearman's rho rank correlation
#'  coefficients for all possible pairs of columns of a matrix. Missing values
#'  are deleted in pairs rather than deleting all rows of x having any missing
#'  variables. Ranks are computed using efficient algorithms, using midranks for
#'  ties. 
#'  
#' @param vec a numeric vector which will be concatenated to the row of 'mat'. 
#' @param mat a numeric matrix with at least 5 rows and at least 2 columns.
#'  For print, x is an object produced by rcorr.
#' @param method specifies the type of correlations to compute. Spearman correlations
#'  are the Pearson linear correlations computed on the ranks of non-missing elements,
#'  using midranks for ties. 
#'  
#' @return the correlation coefficient and p-value.
#'   
#' @seealso \code{\link{rcorr}} which this function wraps.
#' @export
#' @examples
#' mat <- matrix(rnorm(200), ncol=50) 
#' rownames(mat) <- c("A", "B","C", "D")
#' vec <- rnorm(50)
#' rcorr_VecMat(vec=vec, mat=mat, method='pearson')
#' rcorr_VecMat(vec=vec, mat=mat, method='spearman')
rcorr_VecMat <- function(vec, mat, method=c("pearson", "spearman")){
  if(ncol(mat)!=length(vec)) stop(sprintf("mat ncol: %d\nvec length: %d\n", ncol(mat), length(vec)))
  
  res <- t(mapply(rcorr_XY, data.frame(t(mat)),
                  data.frame(vec), showname=rownames(mat), MoreArgs=list(method=method)))
  rank.res <- rank(as.numeric(res[,2]))
  res.out <- data.frame(as.numeric(res[,1]), as.numeric(res[,2]),
                        adj.p.val=p.adjust(as.numeric(res[,2]), "BH"),
                        Rank=rank.res)
  colnames(res.out) <- c(colnames(res), "adj.p.val", "Rank")
  return(res.out)
}





#' Fit Proportional Hazards Regression Model in R  
#'
#' Fits a Cox proportional hazards regression model.  Time dependent 
#' variables, time dependent strata, multiple events per subject, and 
#' other extensions are incorporated using the counting process 
#' formulation of Andersen and Gill.
#'
#' @param expr input a vector of expression values. 
#' @param time for right censored data, this is the follow up time. For 
#'  interval data, the first argument is the starting time for the interval.
#' @param event The status indicator, normally 0=alive, 1=dead. Other choices
#'  are 'TRUE'/'FALSE' ('TRUE' = death) or 1/2 (2=death). 
#'  For interval censored data, the status indicator is 0=right 
#'  censored, 1=event at 'time', 2=left censored, 3=interval 
#'  censored.  Although unusual, the event indicator can be 
#'  omitted, in which case all subjects are assumed to have an event.
#' @param showname default is NULL. It will show the name of the input covariate,
#'  otherwise the character would be used to assign the rowname of the result.
#' @param LRT p-value is return from likelihood ratio test as default. Otherwise is
#'  logrank test version (wald test). 
#'   
#' @return the correlation coefficient and p-value.
#'   
#' @seealso \code{\link{coxph}} and \code{\link{Surv}} which this function wraps.
#' @export
#' @examples
#' hmohiv<-read.table("http://www.ats.ucla.edu/stat/r/examples/asa/hmohiv.csv", sep=",", header = TRUE) 
#' attach(hmohiv)
#' coxph_Expr(expr=age, time=time, event=censor)
#' detach(hmohiv)
coxph_Expr <- function(expr, time, event, showname=NULL, LRT=TRUE){
  fit <- try(coxph(Surv(time, event) ~ expr, na.action=na.omit), silent=TRUE)
  # in case model fails, i.e. no one dies, returns NA
  if(class(fit)=='try-error') {
    res <- data.frame(coef_exp=NA, coef_exp_P=NA, ciL=NA, ciU=NA, p.value=NA)
    iName <- ifelse(is.null(showname), substitute(expr), showname)
    rownames(res) <- iName
    return (res)
  }
  smry <- summary(fit)
  if(LRT){
    pval <- smry$logtest[3]  # likelihood ratio test: more powerful at small sample size
  } else {
    pval <- smry$sctest[3]  # to make it consistent with coef_exp_P; also this is the logrank test version: wald test, performance not good; asymtotically, logrank test is the LRT for cox PH comparing two groups
  }
  cc <- smry$coefficients
  ci <- smry$conf.int
  ciL <- ci[, 'lower .95']
  ciU <- ci[, 'upper .95']
  coef_exp <- cc[, 2]
  coef_exp_P <- cc[, 5]
  res <- data.frame(coef_exp=coef_exp, coef_exp_P=coef_exp_P, ciL=ciL, ciU=ciU, p.value=pval)
  iName <- ifelse(is.null(showname), substitute(expr), showname)
  rownames(res) <- iName
  return (res)
} 



#' Fit Proportional Hazards Regression Model for an expression matrix in R   
#'
#' Fits a Cox proportional hazards regression model.  Time dependent 
#' variables, time dependent strata, multiple events per subject, and 
#' other extensions are incorporated using the counting process 
#' formulation of Andersen and Gill.
#'
#' @param mat input a expression data.frame (probes in row; samples in column). 
#' @param time for right censored data, this is the follow up time. For 
#'  interval data, the first argument is the starting time for the interval.
#' @param event The status indicator, normally 0=alive, 1=dead. Other choices
#'  are 'TRUE'/'FALSE' ('TRUE' = death) or 1/2 (2=death). 
#'  For interval censored data, the status indicator is 0=right 
#'  censored, 1=event at 'time', 2=left censored, 3=interval 
#'  censored.  Although unusual, the event indicator can be 
#'  omitted, in which case all subjects are assumed to have an event.
#' @param LRT p-value is return from likelihood ratio test as default. Otherwise is
#'  logrank test version (wald test).
#'  
#' @return the correlation coefficient and p-value.
#'   
#' @seealso \code{\link{coxph}} and \code{\link{Surv}} which this function wraps.
#' @export
#' @examples
#' mat <- matrix(rnorm(400), ncol=100) 
#' rownames(mat) <- c("A", "B","C", "D")
#' hmohiv<-read.table("http://www.ats.ucla.edu/stat/r/examples/asa/hmohiv.csv", sep=",", header = TRUE) 
#' coxph_VecMat(mat=mat, time=hmohiv$time, event=hmohiv$censor)
coxph_VecMat <- function(mat, time, event, LRT=TRUE){
  if(ncol(mat)!=length(time)) stop(sprintf("mat ncol: %d\ntime length: %d\n", ncol(mat), length(time)))
  
  res <- NULL
  for (i1 in 1:nrow(mat)){
    res <- rbind(res, coxph_Expr(expr=mat[i1,],
                                 time=time, event=event,
                                 showname=rownames(mat)[i1],
                                 LRT=LRT))
  }
  rank.res <- rank(as.numeric(res[,"p.value"]))
  res.out <- data.frame(res,
                        adj.p.val=p.adjust(as.numeric(res[,"p.value"]), "BH"),
                        Rank=rank.res)
  return(res.out)
}  

  











#' Perform Fisher's exact test between each row in a matrix and a vector using `stats` package in R 
#'
#' Performs Fisher's exact test for testing the null of independence of rows and columns
#'  in a contingency table with fixed marginals.
#'   
#' @param group a vector of factor which will be concatenated to the row of 'data'. 
#' @param data a matrix with at least 2 outcomes (0/1; Mut/WT, etc).
#' @param rowlevel specifies the order of factor use in a contingency table.
#' @param ... other parameters in fisher.test function. 
#' 
#' @return the contingency table and p-value. 
#'  
#' @seealso \code{\link{fisher.test}} which this function wraps.
#' @export
#' @examples
#' m0 <- matrix(0, 5, 30)
#' dat <- apply(m0, c(1,2), function(x) sample(c(0,1),1))
#' rownames(dat) <- paste("R", 1:5, sep="")
#' colnames(dat) <- paste("C", 1:30, sep="")
#' si.gp <- factor(c(rep("E", 15), rep("M", 15)), levels=c("E", "M"))
#' stats_fisher.test(dat, si.gp, rowlevel=c("0", "1"))
stats_fisher.test <- function(data, group, rowlevel, ...){
  if(is.factor(group)!=TRUE) stop("group must be a factor!!!")
  
  fe.res <- matrix(0, nrow(data), sum(length(levels(group))*length(rowlevel), 1))
  for (i1 in 1:nrow(data)){
    a <- factor(data[i1,], levels=rowlevel)
    fit <- table(a, group)
    res <- stats::fisher.test(fit, ...)
    fe.res[i1, 1:(ncol(fe.res)-1)] <- as.vector(fit)
    fe.res[i1, ncol(fe.res)] <- res$p.value
  }
  fe.res <- data.frame(fe.res)
  
  my.col <- NULL
  for (x in 1:length(levels(group))){
    for (y in 1:length(rowlevel)){
      my.col <- c(my.col, paste(levels(group)[x], rowlevel[y], sep="_"))
    }
  }
  
  colnames(fe.res) <- c(my.col,
                        "p.value")
  rownames(fe.res) <- rownames(data)
  fe.res <- data.frame(fe.res,
                       adj.p.val=p.adjust(fe.res[,"p.value"], "BH"),
                       Rank=rank(fe.res[,"p.value"]))
  
  #######################
  res.lis <- list()
  res.lis[[1]] <- data
  res.lis[[2]] <- fe.res
  names(res.lis) <- c("dds", "results")
  return(res.lis)
}





#' Perform Pearson's Chi-squared Test between each row in a matrix and a vector using `stats` package in R 
#'
#''chisq.test' performs chi-squared contingency table tests and goodness-of-fit tests.
#' 
#' @param group a vector of factor which will be concatenated to the row of 'data'. 
#' @param data a matrix with at least 2 outcomes (0/1; Mut/WT, etc).
#' @param rowlevel specifies the order of factor use in a contingency table.
#' @param ... other parameters in chisq.test function. 
#' 
#' @return the contingency table and p-value.  
#' 
#' @seealso \code{\link{chisq.test}} which this function wraps.
#' @export
#' @examples
#' m0 <- matrix(0, 5, 30)
#' dat <- apply(m0, c(1,2), function(x) sample(c(0,1),1))
#' rownames(dat) <- paste("R", 1:5, sep="")
#' colnames(dat) <- paste("C", 1:30, sep="")
#' si.gp <- factor(c(rep("E", 15), rep("M", 15)), levels=c("E", "M"))
#' stats_chisq.test(dat, si.gp, rowlevel=c("0", "1"))
stats_chisq.test <- function(data, group, rowlevel, ...){
  if(is.factor(group)!=TRUE) stop("group must be a factor!!!")
  
  fe.res <- matrix(0, nrow(data), sum(length(levels(group))*length(rowlevel), 1))
  for (i1 in 1:nrow(data)){
    a <- factor(data[i1,], levels=rowlevel)
    fit <- table(a, group)
    res <- stats::chisq.test(fit, ...)
    fe.res[i1, 1:(ncol(fe.res)-1)] <- as.vector(fit)
    fe.res[i1, ncol(fe.res)] <- res$p.value
  }
  fe.res <- data.frame(fe.res)
  
  my.col <- NULL
  for (x in 1:length(levels(group))){
    for (y in 1:length(rowlevel)){
      my.col <- c(my.col, paste(levels(group)[x], rowlevel[y], sep="_"))
    }
  }
  
  colnames(fe.res) <- c(my.col,
                        "p.value")
  rownames(fe.res) <- rownames(data)
  fe.res <- data.frame(fe.res,
                       adj.p.val=p.adjust(fe.res[,"p.value"], "BH"),
                       Rank=rank(fe.res[,"p.value"]))
  
  #######################
  res.lis <- list()
  res.lis[[1]] <- data
  res.lis[[2]] <- fe.res
  names(res.lis) <- c("dds", "results")
  return(res.lis)
}






#' Performs Wilcoxon tests between each row in a matrix and a vector using `stats` package in R 
#'
#''Wilcox.test' performs one- and two-sample Wilcoxon tests on vectors of data; the latter
#'  is also known as 'Mann-Whitney' test.
#'   
#' @param group a vector of factor which will be concatenated to the row of 'data'. 
#' @param data a matrix with at least 2 outcomes (0/1; Mut/WT, etc).
#' @param gpMat specifies whether group info is in thr matrix. Default is FALSE.
#' @param ... other parameters in wilcox.test function. 
#' 
#' @return the contingency table and p-value. 
#'  
#' @seealso \code{\link{wilcox.test}} which this function wraps.
#' @export
#' @examples
#' dat <- matrix(rnorm(10000*15*2), 5, 30)
#' rownames(dat) <- paste("R", 1:5, sep="")
#' colnames(dat) <- paste("C", 1:30, sep="")
#' si.gp <- factor(c(rep("E", 15), rep("M", 15)), levels=c("E", "M"))
#' stats_wilcox.test(dat, si.gp, gpMat=FALSE)
#' m0 <- matrix(0, 5, 30)
#' dat <- apply(m0, c(1,2), function(x) sample(c(0,1),1))
#' dat[1,] <- 0
#' dat[2,1] <- NA
#' rownames(dat) <- paste("R", 1:5, sep="")
#' colnames(dat) <- paste("C", 1:30, sep="")
#' si.gp <- rchisq(30,2)
#' stats_wilcox.test(dat, si.gp, gpMat=TRUE) 
stats_wilcox.test <- function(data, group, gpMat=FALSE, ...){
  if (gpMat==FALSE){
    if(is.factor(group)!=TRUE) stop("group must be a factor!!!")
  
    fe.res <- matrix(0, nrow(data), 2)
    for (i1 in 1:nrow(data)){
      res <- stats::wilcox.test(data[i1,] ~ group, na.action=na.omit, ...)

      fe.res[i1, 1] <- res$statistic
      fe.res[i1, 2] <- res$p.value
    }
    fe.res <- data.frame(fe.res)  
    colnames(fe.res) <- c("W.statistic", "p.value")
    rownames(fe.res) <- rownames(data)
    fe.res <- data.frame(fe.res,
                         adj.p.val=p.adjust(fe.res[,"p.value"], "BH"),
                         Rank=rank(fe.res[,"p.value"]))
  }
  if (gpMat==TRUE){
    
    fe.res <- matrix(0, nrow(data), 2)
    for (i1 in 1:nrow(data)){
      if (length(levels(factor(na.omit(data[i1,]))))==2){
        res <- stats::wilcox.test(group ~ data[i1,], na.action=na.omit, ...)
      
        fe.res[i1, 1] <- res$statistic
        fe.res[i1, 2] <- res$p.value
      }
      if (length(levels(factor(na.omit(data[i1,]))))!=2){
        
        fe.res[i1, 1] <- NA
        fe.res[i1, 2] <- NA
      }
    }
         
      fe.res <- data.frame(fe.res)  
      colnames(fe.res) <- c("W.statistic", "p.value")
      rownames(fe.res) <- rownames(data)
      fe.res <- data.frame(fe.res,
                           adj.p.val=p.adjust(fe.res[,"p.value"], "BH"),
                           Rank=rank(fe.res[,"p.value"]))
  }
  
  #######################
  res.lis <- list()
  res.lis[[1]] <- data
  res.lis[[2]] <- fe.res
  names(res.lis) <- c("dds", "results")
  return(res.lis)
}



