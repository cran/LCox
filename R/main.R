#' @title Hypothesis Testing for longitudinal gene profiles and the survival outcome.
#' @description \code{LCox} is a function that performs hypothesis tests of the associations
#' between longitudinal gene profiles and the survival outcome.
#' @param data A data frame contains longitudinal gene expression data. This data
#' frame must contain a column of \code{ID} to identify each patient and a column of
#' \code{years} to indicate the follow-up time for each data point.
#' @param data.id A data frame contains the survival outcome and important covariates.
#' This data frame must contain a column of \code{ID} to identify each patient,
#' a column of \code{fstat} to indicate the censoring time, and a column of \code{ftime}
#' to indicate the survival time.
#' @param geneID A vector of integers indicates the column numbers in \code{data}
#' for the genes of interest.
#' @param varID A vector of integers indicates the column numbers in \code{data.id}
#' for the important confounding covariates that need to be included.
#' @param PLOT A logical value indicates whether a graph showing the fitted lines is desired.
#' If TRUE, a figure "fitted.pdf" will be saved to the current working directory.
#' @param optns A list of options control parameters for the FPCA model.
#' @return returns a matrix with one column being the p values and the other
#' column being the number of eigenfunctions (K).
#' @references LCox: A tool for selecting genes related to survival outcomes
#' using longitudinal gene expression data. Jiehuan Sun, Jose D. Herazo-Maya,
#' Jane-Ling Wang, Naftali Kaminski, and Hongyu Zhao.
#' @export
#' @examples
#' data.list = simudata()
#' data = data.list$data
#' data.id = data.list$data.id
#' res = LCox(data = data, data.id = data.id, geneID = 3:4)
#' res = LCox(data = data, data.id = data.id, geneID = 3:4, varID = 4)
#' res = LCox(data = data, data.id = data.id, geneID = 3:4, PLOT=TRUE)
LCox <- function(data=NULL, data.id=NULL, geneID=3:4, varID = NULL, PLOT = FALSE,
                     optns=list(dataType='Sparse', FVEthreshold=0.95, methodBwMu='CV')){

    Tmax = max(data$years)

    #### separate into two groups ####
    iiid = (1:nrow(data.id))[data.id$ftime < Tmax]

    data2 = data[is.element(data$ID,iiid),]
    data.id2 = data.id[iiid,]
    if(length(iiid)>0){
        data.id = data.id[-iiid,]
    }

    data = data[!is.element(data$ID,iiid),]
    data.id = rbind(data.id,data.id2)

    genename <- colnames(data)[geneID]
    iter = length(genename)


    TT <- lapply(unique(data$ID),function(x){
        data$years[data$ID==x]
    })

    TT2 <- lapply(unique(data2$ID),function(x){
        data2$years[data2$ID==x]
    })

    res <- sapply(1:iter,function(i){

        XY = lapply(unique(data$ID),function(x){
            data[data$ID==x,genename[i]]
        })
        XY2 = lapply(unique(data2$ID),function(x){
            data2[data2$ID==x,genename[i]]
        })

        FPC = FPCA(XY,TT,optns=optns)

        if(length(TT2) > 0){

            toGrid = sort(unique(data2$years))
            mu = ConvertSupport(fromGrid=FPC$workGrid, toGrid = toGrid, mu = FPC$mu)
            cov = ConvertSupport(fromGrid=FPC$workGrid, toGrid = toGrid, Cov = FPC$fittedCov)
            phi = ConvertSupport(fromGrid=FPC$workGrid, toGrid = toGrid, phi = FPC$phi)
            #phi = as.matrix(phi)
            cov = cov + diag(FPC$rho,length(toGrid),length(toGrid))
            scores = sapply(1:length(TT2),function(ii){
                TID = sapply(TT2[[ii]],function(x){which(abs(toGrid-x)<1e-8 ) })
                if(is.matrix(phi)){
                    diag(FPC$lambda,length(FPC$lambda))%*%(t(phi[TID,,drop=FALSE])%*%solve(cov[TID,TID,drop=FALSE])%*%(XY2[[ii]]-mu[TID]))
                }else{
                    diag(FPC$lambda,length(FPC$lambda))%*%(phi)%*%solve(cov[TID,TID,drop=FALSE])%*%(XY2[[ii]]-mu[TID])
                }

            })

            if(is.matrix(scores)){
                scores = t(scores)
            }else{scores = as.matrix(scores)}


            gene=FPC$xiEst
            gene = rbind(gene,scores)
        }else{
            gene=FPC$xiEst
        }

        if(PLOT){
            plotFF(data = data, data2 = data2, gene=gene, FPC=FPC)
        }

        #print(dim(gene))
        #print(Tmax)
        #print(dim(data.id2))
        #print(FPC$sigma2)

        if(is.null(varID)){
            # fit0 <- coxph(Surv(ftime, fstat) ~ 1, data = data.id)
            fit1 = coxph(Surv(ftime, fstat) ~ gene, data = data.id)
            xx = anova(fit1)
            c(xx["Pr(>|Chi|)"][2,1],ncol(gene))
        }else{
            names = colnames(data.id)[varID]
            form0 = as.formula(paste("Surv(ftime, fstat) ~",paste(names,collapse="+")))
            form1 = as.formula(paste("Surv(ftime, fstat) ~ gene + ",paste(names,collapse="+")))
            fit0 = coxph(form0, data = data.id)
            fit1 = coxph(form1, data = data.id)
            xx = anova(fit0,fit1)
            c(xx["P(>|Chi|)"][2,1],ncol(gene))
        }

    })

    res = t(res)
    colnames(res) = c("Pvalues", "K")
    res
}

#' @title Plot Function
#' @description \code{plotFF} generates the fitted lines on the data using FPCA.
#' @keywords internal
#' @param data Longitudinal gene expression data for type one patients.
#' @param data2 Longitudinal gene expression data for type two patients.
#' @param gene Functional scores for the patients.
#' @param FPC A FPCA object returned by the FPCA function.
plotFF <- function(data = NULL, data2 = NULL, gene=NULL, FPC=NULL){

    tmp = rbind(data,data2)
    fittedline = FPC$phi %*% t(gene) + FPC$mu

    yli = c(min(tmp[,3],min(fittedline)),max(tmp[,3],max(fittedline)))
    xli = c(min(tmp[,2]),max(tmp[,2]))

    pdf(file="fitted.pdf", height = 4, width = 8)
    layout(matrix(1:2,ncol=2))
    plot(xli[1],yli[1],xlim = xli,ylim = yli,col = "white",
         xlab = "years", ylab="Orignal gene profiles",bty="L")
    for(id in unique(tmp$ID)){
        lines(tmp$years[tmp$ID==id],tmp[tmp$ID==id,3])
    }
    # lines(FPC$workGrid, FPC$mu, col="red")

    plot(xli[1],yli[1],xlim = xli,ylim = yli,col = "white",
         xlab = "years", ylab="Smoothed gene profiles",bty="L")
    for(i in 1:ncol(fittedline)){
        lines(FPC$workGrid, fittedline[,i])
    }
    dev.off()
}

#' @title Function to Simulate Testing Data
#' @description \code{simudata} is a function that generates a simulated data.
#' @param n Total number of patients.
#' @param nf The average number of follow-up visits per patient.
#' @return returns a list with following objects.
#' \item{data}{Longitudinal gene expression data.}
#' \item{data.id}{Survival data with important covariates.}
#' @references LCox: A tool for selecting genes related to survival outcomes
#' using longitudinal gene expression data. Jiehuan Sun, Jose D. Herazo-Maya,
#' Jane-Ling Wang, Naftali Kaminski, and Hongyu Zhao.
#' @export
#' @examples
#' data.list = simudata()
simudata <- function(n = 50, nf =5){

    nt = sample(2:nf,n,replace=TRUE)
    ID = rep(1:n,nt)
    lam = 3
    theta = 10

    AB1 = cbind(runif(n,0,4),runif(n,0,8))
    AB2 = cbind(runif(n,0,4),runif(n,0,8))

    Y = rep(0,n)
    for(i in 1:n){
        Atrue = AB1[i,1]
        Btrue = AB1[i,2]

        logr = log(runif(1))
        ff = function(t){(lam/theta^lam)*t^(lam-1)*exp(Atrue/8 + Btrue*t/8)}
        f = function(x){ integrate(ff,0,x)$value + logr }
        Y[i] = uniroot(f,c(0,50))$root
    }

    C = rweibull(n,3,8.5)
    U = pmin(Y,C)
    Delta = as.numeric(Y<C)
    data.id = data.frame(ID=1:n,fstat=Delta,ftime=U)

    S = sapply(1:n,function(i){c(0,sort(runif(nf-1,0,1))) })
    data = NULL
    years = NULL

    for(i in 1:n){
        tmp1 = AB1[i,1]+AB1[i,2]*S[,i] + rnorm(nf,0,1)
        tmp2 = AB2[i,1]+AB2[i,2]*S[,i] + rnorm(nf,0,1)

        idx = sort(c(1,sample(2:nf,nt[i]-1)))
        years = c(years,S[idx,i])
        data = rbind(data,cbind(tmp1[idx],tmp2[idx]) )
    }

    SS = U[ID]
    ID = ID[SS > years]
    data = data[SS > years,]
    years = years[SS > years]
    data = data.frame(ID=ID,years=years,data=data)
    colnames(data)[1:2 +2] = paste("Gene",1:2,sep="")

    genename = colnames(data)[3:ncol(data)]
    iter = length(genename)
    TT = lapply(unique(data$ID),function(x){
        data$years[data$ID==x]
    })

    data.id$age = runif(nrow(data.id), min=40, max = 80)

    list(data = data, data.id = data.id)
}
