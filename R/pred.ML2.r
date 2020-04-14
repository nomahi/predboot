pred.ML2 <- function(model,data,B=1000,alpha=0.05,Cores=detectCores()){

  print("This computation will take a lot of time. Please don't stop the computation at least 10 to 20 minutes to obtain the outputs.")

  opts <- list(progress = function(x) print(paste0("The ",x,"th bootstrap is completed.")))

　　N <- dim(data)[1]

  cl <- makeSOCKcluster(Cores)
  registerDoSNOW(cl)

　　boot <- function(model, data){
  
  　　N <- dim(data)[1]
    
  　　bs.i <- sample(1:N, N, replace=TRUE)
 　　 train.i <- data[ bs.i,]
 　　 test.i  <- data[-bs.i,]
 
 　　 gm1.i <- glm(model, data=train.i, family=binomial)
  
 　　 prob.b <- predict(gm1.i, type=c("response"))
 　　 prob.o <- predict(gm1.i, newdata=data,type = c("response"))
 　　 prob.t <- predict(gm1.i, newdata=test.i,type = c("response"))
  
 　　 AUC.b <- roc(train.i$Y ~ prob.b, levels=c(0,1), direction="<")$auc
 　　 AUC.o <- roc(data$Y ~ prob.o, levels=c(0,1), direction="<")$auc
 　　 AUC.t <- roc(test.i$Y ~ prob.t , levels=c(0,1), direction="<")$auc

 　　 return(c(AUC.b, AUC.o, AUC.t))

　　}
       	 
  R1 <- foreach(b = 1:B, .combine = rbind, .packages=c("MASS","pROC"), .options.snow = opts) %dopar% {

   re.i <- sample(1:N,N,replace=TRUE)
   train <- data[re.i,]
	
   gm1 <- glm(model, data=train, family=binomial)

   prob <- predict(gm1,type=c("response"))
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

