pred.lasso <- function(model,data,B=1000,alpha=0.05,Cores=detectCores()){

  print("This computation will take a lot of time. Please don't stop the computation at least several hours to obtain the outputs.")

　　N <- dim(data)[1]

　　boot <- function(model, data){
  
  　　N <- dim(data)[1]
    
  　　bs.i <- sample(1:N, N, replace=TRUE)
 　　 train.i <- data[ bs.i,]
 　　 test.i  <- data[-bs.i,]
 
    X.b <- model.matrix(model,data=train.i)
    X.o <- model.matrix(model,data=data)
    X.t <- model.matrix(model,data=test.i)
 
    gm1.i <- cv.glmnet(x=X.b[,-1], y=train.i$Y, family="binomial", alpha=1, foldid=train.i$ID, grouped=FALSE)
 
    coef.i <- coef(gm1.i, s="lambda.min")
  
    prob.b <- as.vector( 1/(1+exp(-X.b%*%coef.i)) )
    prob.o <- as.vector( 1/(1+exp(-X.o%*%coef.i)) )
    prob.t <- as.vector( 1/(1+exp(-X.t%*%coef.i)) )
  
    AUC.b <- roc(train.i$Y ~ prob.b, levels=c(0,1), direction="<")$auc
    AUC.o <- roc(data$Y ~ prob.o  , levels=c(0,1), direction="<")$auc
    AUC.t <- roc(test.i$Y ~ prob.t , levels=c(0,1), direction="<")$auc

 　　 return(c(AUC.b, AUC.o, AUC.t))

　　}
  
  X.train <- model.matrix(model,data=data)
  
  gm1 <- cv.glmnet(x=X.train[,-1], y=data$Y, family="binomial", alpha=1, foldid=data$ID, grouped=FALSE)
  coef <- coef(gm1, s="lambda.min")
   

  #calculation apparent AUC & test AUC
  prob <- as.vector( 1/(1+exp(-X.train%*%coef)) )
  ROC.app <- roc(data$Y ~ prob, levels=c(0,1), direction="<")
  AUC.app <- ROC.app$auc
  delong1 <- ci.auc(ROC.app,conf.level=1-alpha)[1]
  delong2 <- ci.auc(ROC.app,conf.level=1-alpha)[3]
  
  cl <- makeSOCKcluster(Cores)
  registerDoSNOW(cl)

  block <- ceiling(B/Cores)
  block0 <- c(1, block*(1:(Cores-1)) + 1)
  block1 <- c(block*(1:(Cores-1)), B)

  boot.res <- foreach(b = 1:Cores, .combine = rbind, .packages=c("MASS","pROC","glmnet")) %dopar% {

   R1 <- NULL
   for(iter in block0[b]:block1[b]) R1 <- rbind(R1, boot(model, data))
   R1

  }

  #bootstrap SD and 95%CI of apparent C
  boot.app.C <- boot.res[,1]
	
  # bias corrected AUC estimate (ordinary bootstrap)
  AUC.boot <- AUC.app - mean(boot.res[,1] - boot.res[,2])	

  AUC.boot.CL1 <- quantile(boot.app.C, 0.5*alpha)
  AUC.boot.CL2 <- quantile(boot.app.C, 1-0.5*alpha)
    
  # bias corrected AUC estimate (bootstrap .632)
  AUC.loocv <- mean(boot.res[,3])	
  AUC.632 <- 0.368*AUC.app + 0.632*AUC.loocv
  AUC.loocv.SD <- sd(boot.res[,3])
    
  # bias corrected AUC estimate (bootstrap .632+)
  if (AUC.loocv<=0.5){
  	R <- 1
  } else if (AUC.app > AUC.loocv){
 	R <- (AUC.app - AUC.loocv)/(AUC.app - 0.5)
  } else {
	R <- 0
  }
    
  w <- 0.632/(1-0.368*R)
  AUC.632p <- (1-w)*AUC.app + w*max(AUC.loocv,0.5) 
  
  delta1 <- AUC.boot - AUC.app
  delta2 <- AUC.632 - AUC.app
  delta3 <- AUC.632p - AUC.app
  
  R <- list(
  　glm.output=gm1,
  　N.obs=N,
   N.boot=B,
   C.Apparent=as.numeric(AUC.app),
   C.DeLongCI=c(delong1,delong2),
   C.Apparent_BootstrapCI=c(AUC.boot.CL1,AUC.boot.CL2),
   C.Harrell=AUC.boot,
   C.Harrell_LSCI=(c(AUC.boot.CL1,AUC.boot.CL2) + delta1),
   C.0.632=AUC.632,
   C.0.632_LSCI=(c(AUC.boot.CL1,AUC.boot.CL2) + delta2),
   C.0.632p=AUC.632p,
   C.0.632p_LSCI=(c(AUC.boot.CL1,AUC.boot.CL2) + delta3)
   )
   
  stopCluster(cl)
  
  return(R)

} 
