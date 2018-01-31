
#' model optimization for fitting exponential decay models to normalized data
#'
#' The modOptimization function finds the estimates of model parameters by maximum likelihood, for a single gene on a specified list of models, and saves a tab delimited text file of the results named,"[geneID]_results.txt".
#' The function does the following for each gene:
#' (1) it calculates log likelihood for each point in a 2 dimensional grid of evenly spaced alpha and beta values within the alpha and beta bounds specified using the null model (in which all treatment alphas are equivalent and all betas are equivalent).
#' (2) it calculates log likelihood for each point in a 1 dimensional range of evenly spaced alpha values within the alpha bounds using the single exponential null model (in which all treatment alphas are equivalent).
#' (3) For each of the grid points with the highest log likelihood from steps (1) and (2) 25 starting parameter value sets that are normally distributed around these points are generated.
#' (4) Parameter values are optimized for maximum likelihood using each of these 50 starting parameter sets using pre-compiled C++ functions loaded from dynamically linked libraries stored in the package on all models specified in the models argument.
#' (5) evaluates parameter estimates of all 50 optimizations based on the reported maximum liklihood upon convergence. Only parameter estimates that converged on the same and highest maximum likelihood are returned.
#' (6) returns the optimized parameter estimates, with model selection statistics.
#'
#' @param gene geneID from \code{data} to be modeled
#' @param data decay data data.frame with columns including
#' @param models vector spceifying which models to run optimization on (e.g., c("mod1", "mod239"))
#' @param alpha.bounds vector of length 2 with lower and upper bounds for alpha
#' @param beta.bounds vector of length 2 with lower and upper bounds for beta
#' @param group grouping matrix for alphas or betas
#' @param mod data.frame specifying alpha and beta group pairs for each model
#' @param file.only logical; should output only be written to file (TRUE) or also return a data.frame of the results (FALSE)
#' @param path specify folder for output to be written
#'
#' @export
#'
#' @return returns (if \code{file.only = FALSE}) and writes to \code{path} a data frame of model optimization results for \code{models} one row for each for
#' \code{gene} using values for it found in \code{data}, the columns of the data frame are:
#' geneID, mod (model), model estimates [alpha_treatment1, ..., alpha_treatmentn, beta_treatment1, ..., beta_treatmentn, sigma2],	logLik (maximum log likelihood),
#' nPar (number of parameters in the model), nStarts (number of parameter starting value sets (of 50) that converged on a maximum likelihood peak), J (number of
#' parameter starting value sets that converged on the highest - within 1e-4 - maximum likelihood of all parameter starting value sets), range.LL (range of maximum
#' likelihoods values reached by algorithm convergence from all parameter starting value sets),	nUnique.LL (number of unique maximum likelihoods values reached by
#' algorithm convergence from all parameter starting value sets), C.alpha (sum of all coefficients of variation for each column of alpha estimates),
#' C.beta (sum of all coefficients of variation for each column of beta estimates), C.tot (C.alpha+C.beta),	AICc (calculated from the single highest maximum
#' likelihood of all parameter starting value sets),	AICc_est (calculated from the log likelihood value computed by using the mean of each parameter from all
#' optimizations that converged on the highest maximum likelihood of all starting parameter value sets.)
#'
#' @examples
#' wd = getwd()
#' setwd(paste0(find.package("RNAdecay"), "/src"))
#' if(file.exists(TMB::dynlib("general_Exp_2sse"))) {
#' unlink(c(TMB::dynlib("general_Exp_2sse"),"general_Exp_2sse.o"))}
#' if(file.exists(TMB::dynlib("general_dExp_2sse"))) {
#' unlink(c(TMB::dynlib("general_dExp_2sse"),"general_dExp_2sse.o"))}
#' TMB::compile("general_Exp_2sse.cpp")
#' TMB::compile("general_dExp_2sse.cpp")
#' dyn.load(TMB::dynlib("general_Exp_2sse"))
#' dyn.load(TMB::dynlib("general_dExp_2sse"))
#' setwd(wd); rm(wd)
#'
#' modOptimization(gene = "Gene_BooFu",
#'                 data = data.frame(geneID=rep("Gene_BooFu",30),
#'                             treatment=c(rep("WT",15),rep("mut",15)),
#'                             t.decay=rep(c(0,7.5,15,30,60),6),
#'                             rep=rep(paste0("rep",c(rep(1,5),rep(2,5),rep(3,5))),2),
#'                             value= c(0.9173587, 0.4798672, 0.3327807, 0.1990708, 0.1656554,
#'                                      0.9407511, 0.7062988, 0.3450886, 0.3176824, 0.2749946,
#'                                      1.1026497, 0.6156978, 0.4563346, 0.2865779, 0.1680075,
#'                                      0.8679866, 0.6798788, 0.2683555, 0.5120951, 0.2593122,
#'                                      1.1348219, 0.8535835, 0.6423996, 0.5308946, 0.4592902,
#'                                      1.1104068, 0.5966838, 0.3949790, 0.3742632, 0.2613560)),
#'                 alpha.bounds = c(1e-4,0.75),
#'                 beta.bounds = c(1e-3,0.075),
#'                 models = "mod1",
#'                 group = t(matrix(c(1,2,1,1,NA,NA),nrow=2,
#'                           dimnames=list(c("treat1","treat2"),c("mod1","mod2","mod3")))),
#'                 mod = as.data.frame(t(matrix(c(1,1,1,2,1,3,2,1,2,2,2,3),nrow=2,
#'                         dimnames=list(c("a","b"),paste0("mod",1:6))))),
#'                 file.only = FALSE)
#'


modOptimization = function(gene, data, alpha.bounds, beta.bounds, models, group, mod,file.only=TRUE, path="modeling_results") {
  if (!file.exists(path)) {dir.create(path)}

  genoSet = 1:(length(unique(data$rep)) * length(unique(data$t.decay)))
  nTreat=length(unique(data$treatment))
  try(if(sum(grepl(paste0("_",nTreat,"sse"),names(getLoadedDLLs())))!=2) stop("high performance objective functions need to be loaded as compilied dynamically linked libraries (.dll or .so) before modeling. You can check which are loaded using 'getLoadedDLLs()'", call. = FALSE))
  nSet=length(genoSet)*nTreat

  try(if(nTreat>4) stop("modOptimization can only handle up to 4 treatments."))

  gene = as.character(gene)

  # pull out the data the specific gene data
  gdata=data[data$geneID==as.character(gene),]
  t <- gdata$t.decay
  m <- gdata$value
  const = constraint_fun_list_maker(mods = mod,groups = group)
  eps = 1e-4

  #### Determine starting points. ####

  ###Define grid for evaluating model 240 while determining starting points.
  A <- seq(0,alpha.bounds[2],by=1e-3)
  ###Define grid for evaluating model 239 while determining starting points.
  a <- seq(alpha.bounds[1],alpha.bounds[2],by=1e-3)
  b <- seq(beta.bounds[1],beta.bounds[2],by=1e-3)

  ### Pick 25 points near the "peak" of Model null double exponential model.
  ### Find the SSE surface for the double exponential model...
  sse.nulldExp <- sapply(b,function(x){
    sapply(a,FUN=sse.nulldExp.fun,b=x,m=m,t=t)
  })

  ###Find the location of the minimum.
  loc239 <- which(sse.nulldExp==min(sse.nulldExp), arr.ind=TRUE)

  ###Pick 25 points near the "peak" of Model 240.
  ###Find the SSE surface for model 240...
  sse.nullExp <- sapply(A,FUN=sse.nullExp.fun,m=m,t=t)
  ###Find the location of the minimum.
  loc240 <- which(sse.nullExp==min(sse.nullExp))

  ###Create the matrix of starting points.
  ###Pick the alpha values.
  aX0 <- matrix(c(stats::rnorm(25*nTreat,mean=a[loc239[1,"row"]],sd=0.01), stats::rnorm(25*nTreat,mean=A[loc240],sd=0.01)), ncol=nTreat)
  ###Fix the alpha values outside of [1e-4,0.75]
  aX0 <- ifelse(aX0<alpha.bounds[1],alpha.bounds[1],aX0)
  aX0 <- ifelse(aX0>alpha.bounds[2],alpha.bounds[2],aX0)

  ###Pick the beta values.
  bX0 <- matrix(c(stats::rnorm(25*nTreat,mean=b[loc239[1,"col"]],sd=0.01), stats::runif(25*nTreat,beta.bounds[1],beta.bounds[2])), ncol=nTreat)
  ###Fix the beta values outside of [1e-3,0.075]
  bX0 <- ifelse(bX0<beta.bounds[1],beta.bounds[1],bX0)
  bX0 <- ifelse(bX0>beta.bounds[2],beta.bounds[2],bX0)

  ###Combine alpha and beta values into one matrix.
  X0 <- cbind(aX0,bX0)
  colnames(X0) <- c(paste0("alpha.int",1:nTreat),paste0("beta.int",1:nTreat))

  #set default parameters
  par.default <- c(rep(0,nTreat),rep(1,nTreat))
  names(par.default) <- c(paste0("a",1:nTreat),paste0("b",1:nTreat))


  # OPTIMIZATION OF THE DOUBLE EXPONENTIAL MODELS

  #### Create the objective function in R. ### This has to be done for each gene! ####
  obj.dExp = TMB::MakeADFun(data = list(t = t, m = m),parameters = par.default,silent=TRUE,
                            DLL=  paste0("general_dExp_",nTreat,"sse"))

  results4ab <- sapply(models[models %in% rownames(mod)[mod$b != max(mod$b)]],
                       function(x,eps,nSet,md,grp,alpha.bounds,beta.bounds){
                         fits <- apply(X0,1,function(starts,y) {
                           fit <- nloptr::slsqp(x0=as.numeric(starts[1:ncol(X0)]),fn=obj.dExp$fn,gr=obj.dExp$gr,heq=const[[x]],
                                       lower = c(rep(alpha.bounds[1],nTreat),rep(beta.bounds[1],nTreat)), upper = c(rep(alpha.bounds[2],nTreat),rep(beta.bounds[2],nTreat)))
                           unlist(c(geneID=gene,mod=x,as.numeric(fit$par),fit$value,fit$convergence))
                         },y=x)

                         fits = data.frame(t(fits))
                         fits[,-c(1:2)] = sapply(fits[,-c(1,2)],function(x)as.numeric(as.character(x)))
                         fits = as.data.frame(fits)
                         colnames(fits) =  c("geneID","mod",c(paste0("alpha_",as.character(unique(gdata$treatment))),paste0("beta_",as.character(unique(gdata$treatment)))),"SSE","conv.code")

                         fits=fits[fits$conv.code == 4,]
                         fits$sigma2 = sig.fun(x=fits$SSE,n=nSet)
                         fits$logLik = LL.fun(x=fits$SSE,y=fits$sigma2,n=nSet)

                         max.LL= max(fits$logLik)
                         range.LL = max.LL - min(fits$logLik)
                         n.LL = length(unique(round(fits$logLik,4)))
                         tmp = fits[fits$logLik > (max.LL - eps),]
                         C.alpha <- C_eps.fun(tmp[,grep("alpha",colnames(fits))])
                         C.beta <- C_eps.fun(tmp[,grep("beta",colnames(fits))])
                         C.tot <- C.alpha+C.beta
                         par.est= colMeans(tmp[,c(grep("alpha",colnames(fits)),grep("beta",colnames(fits)))])
                         sigma2=mean(tmp$sigma2)
                         nPar=fun_nPar(x,mod=md,group=grp)
                         AICc = fun_aicc(max.LL,nPar,nSet)
                         AICc_est = fun_aicc(LL.fun(x = obj.dExp$fn(par.est),y = sig.fun(x = obj.dExp$fn(par.est), n = nSet), n = nSet), p = nPar, n = nSet)
                         fit = c(as.character(fits$geneID[1]),as.character(fits$mod[1]),par.est,sigma2,max.LL,nPar,nrow(fits),nrow(tmp),range.LL,n.LL,C.alpha,C.beta,C.tot,AICc,AICc_est)
                         names(fit) = c("geneID","mod",c(paste0("alpha_",as.character(unique(gdata$treatment))),paste0("beta_",as.character(unique(gdata$treatment)))),
                                        "sigma2","logLik","nPar","nStarts","J","range.LL","nUnique.LL","C.alpha","C.beta","C.tot","AICc","AICc_est")
                         return(fit)
                       },eps=eps,nSet=nSet,grp=group,md=mod,alpha.bounds=alpha.bounds,beta.bounds=beta.bounds)
  results4ab = as.data.frame(t(results4ab))
  results4ab[,-c(1,2)] = sapply(results4ab[,-c(1,2)],function(x) as.numeric(as.character(x)))

  # OPTIMIZATION OF THE SINGLE EXPONENTIAL MODELS

  X0 <- aX0
  colnames(X0) <- c(paste0("alpha.int",1:nTreat))

  obj.Exp = TMB::MakeADFun(data = list(t = t, m = m),parameters = par.default[1:nTreat],silent=TRUE,
                           DLL = paste0("general_Exp_",nTreat,"sse"))

  results4a <- sapply(models[models %in% rownames(mod)[mod$b == max(mod$b)]],
                      function(x,eps,nSet,md,grp,alpha.bounds){

                        fits <- apply(X0,1,function(starts,y) {
                          fit <- nloptr::slsqp(x0=as.numeric(starts[1:ncol(X0)]),fn=obj.Exp$fn,gr=obj.Exp$gr,heq=const[[x]],
                                      lower = rep(alpha.bounds[1],nTreat), upper = rep(alpha.bounds[2],nTreat))
                          unlist(c(geneID=gene,mod=x,as.numeric(fit$par),rep(0,nTreat),fit$value,fit$convergence))
                        },y=x)

                        fits = data.frame(t(fits))
                        fits[,-c(1:2)] = sapply(fits[,-c(1,2)],function(x)as.numeric(as.character(x)))
                        fits = as.data.frame(fits)
                        colnames(fits) =  c("geneID","mod",c(paste0("alpha_",as.character(unique(gdata$treatment))),paste0("beta_",as.character(unique(gdata$treatment)))),"SSE","conv.code")

                        fits=fits[fits$conv.code == 4,]
                        fits$sigma2 = sig.fun(x=fits$SSE,n=nSet)
                        fits$logLik = LL.fun(x=fits$SSE,y=fits$sigma2,n=nSet)

                        max.LL= max(fits$logLik)
                        range.LL = max.LL - min(fits$logLik)
                        n.LL = length(unique(round(fits$logLik,4)))
                        tmp = fits[fits$logLik > (max.LL - eps),]
                        C.alpha <- C_eps.fun(tmp[,grep("alpha",colnames(fits))])
                        C.beta <- C_eps.fun(tmp[,grep("beta",colnames(fits))])
                        C.tot <- C.alpha+C.beta
                        par.est= colMeans(tmp[,c(grep("alpha",colnames(fits)),grep("beta",colnames(fits)))])
                        sigma2=mean(tmp$sigma2)
                        nPar=fun_nPar(x,mod=md,group=grp)
                        AICc = fun_aicc(max.LL,nPar,nSet)
                        AICc_est = fun_aicc(LL.fun(x = obj.Exp$fn(par.est[1:nTreat]),y = sig.fun(x = obj.Exp$fn(par.est[1:nTreat]), n = nSet), n = nSet), p = nPar, n = nSet)
                        fit = c(as.character(fits$geneID[1]),as.character(fits$mod[1]),par.est,sigma2,max.LL,nPar,nrow(fits),nrow(tmp),range.LL,n.LL,C.alpha,C.beta,C.tot,AICc,AICc_est)
                        names(fit) = c("geneID","mod",c(paste0("alpha_",as.character(unique(gdata$treatment))),paste0("beta_",as.character(unique(gdata$treatment)))),
                                       "sigma2","logLik","nPar","nStarts","J","range.LL","nUnique.LL","C.alpha","C.beta","C.tot","AICc","AICc_est")
                        return(fit)
                      },eps=eps,nSet=nSet,grp=group,md=mod,alpha.bounds=alpha.bounds)
  results4a = as.data.frame(t(results4a))
  results4a[,-c(1,2)] = sapply(results4a[,-c(1,2)],function(x) as.numeric(as.character(x)))

  results=rbind(results4a,results4ab)
  results=results[order(as.numeric(gsub("mod","",results$mod))),]

  utils::write.table(results, paste0(path,"/",gene,"_results.txt"),sep="\t")
  cat(gene,"done     \n")#; utils::timestamp()

return( if(file.only) invisible(NULL) else results)

}

