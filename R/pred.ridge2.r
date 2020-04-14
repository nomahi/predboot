pred.ridge2 <- function(model,data,B=1000,alpha=0.05,Cores=detectCores()){

  print("This computation will take a lot of time. Please don't stop the computation at least several hours to obtain the outputs.")

  opts <- list(progress = function(x) print(paste0("The ",x,"th bootstrap is completed.")))

　　N <- dim(data)[1]

  cl <- makeSOCKcluster(Cores)
  registerDoSNOW(cl)

　　boot <- function(model, data){
  
  　　N <- dim(data)[1]
    
  　　bs.i <- sample(1:N, N, replace=TRUE)
 　　 train.i <- data[ bs.i,]
 　　 test.i  <- data[-bs.i,]
 
    X.b <- model.matrix(model,data=train.i)
    X.o <- model.matrix(model,data=data)
    X.t <- model.matrix(model,data=test.i)
 
    gm1.i <- cv.glmnet(x=X.b[,-1], y=train.i$Y, family="binomial", alpha=0, foldid=train.i$ID, grouped=FALSE)
 
    coef.i <- coef(gm1.i, s="lambda.min")
  
    prob.b <- as.vector( 1/(1+exp(-X.b%*%coef.i)) )
    prob.o <- as.vector( 1/(1+exp(-X.o%*%coef.i)) )
    prob.t <- as.vector( 1/(1+exp(-X.t%*%coef.i)) )
  
    AUC.b <- roc(train.i$Y ~ prob.b, levels=c(0,1), direction="<")$auc
    AUC.o <- roc(data$Y ~ prob.o  , levels=c(0,1), direction="<")$auc
    AUC.t <- roc(test.i$Y ~ prob.t , levels=c(0,1), direction="<")$auc

 　　 return(c(AUC.b, AUC.o, AUC.t))

　　}
       	 
  R1 <- foreach(b = 1:B, .combine = rbind, .packages=c("MASS","pROC","glmnet"), .options.snow = opts) %dopar% {

   re.i <- sample(1:N,N,replace=TRUE)
   train <- data[re.i,]
	
   X.train <- model.matrix(model,data=train)
   
   train$ID <- 1:N
  
   gm1 <- cv.glmnet(x=X.train[,-1], y=train$Y, family="binomial", alpha=0, foldid=train$ID, grouped=FALSE)
   coef <- coef(gm1, s="lambda.min")

   prob <- as.vector( 1/(1+exp(-X.train%*%coef)) )
   ROC.app <- roc(train$Y ~ prob, levels=c(0,1), direction="<")
   AUC.app <- ROC.app$auc

   boot.res <- NULL
   for(iter in 1:B) boot.res <- rbind(boot.res, boot(model, data=train))

   # bias corrected AUC estimate (ordinary bootstrap)
   AUC.boot <- AUC.app - mean(boot.res[,1] - boot.res[,2])	

   # bias corrected AUC estimate (bootstrap .632)
   AUC.loocv <- mean(boot.res[,3])	
   AUC.632 <- 0.368*AUC.app + 0.632*AUC.loocv
    
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

   c(AUC.boot,AUC.632,AUC.632p)
  
  }

  Q1 <- quantile(R1[,1],c(0.5*alpha,1-0.5*alpha))
  Q2 <- quantile(R1[,2],c(0.5*alpha,1-0.5*alpha))
  Q3 <- quantile(R1[,3],c(0.5*alpha,1-0.5*alpha))

  R <- list(
  　N.obs=N,
   N.boot=B,
   C.Harrell_2BSCI=Q1,
   C.0.632_2BSCI=Q2,
   C.0.632p_2BSCI=Q3
   )
   
  stopCluster(cl)
  
  return(R)

} 

