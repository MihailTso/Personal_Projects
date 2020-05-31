#Black Litterman model

#Download Sp500 stocks
library(rvest)
library(TTR)
library(plyr)
library(quantmod)
library(alphavantager)

setwd('/Users/mihailtsonev/Desktop/Research_03_14/')
library(readr)
sp500 <- read.csv("sp500_symbols.csv")

data=as.data.frame(sp500)


data[,'year']=as.character(substr(data$date,1,4))
data[,'month']=as.character(substr(data$date,5,6))
data[,'month']=sprintf("%02d", as.numeric(data$month))
data$date=paste(data$year,data$month,'01',sep = '-')

data=data[order(data$date),]
data_temp=as.matrix(data[,2:507])

data_temp=apply(data_temp, 2, as.numeric)
returns=diff(data_temp)/data_temp[-nrow(data_temp),]
rownames(returns)=data[1:nrow(data)-1,1]

for(i in 1:ncol(returns)){
  returns[is.na(returns[,i]), i] <- mean(returns[,i], na.rm = TRUE)
}

returns=returns[, colSums(is.na(returns)) != nrow(returns)]
rm(data,data_temp,sp500)


all_stocks=stockSymbols(exchange = c("AMEX", "NASDAQ", "NYSE"))
information=all_stocks[which(all_stocks$Symbol %in% colnames(returns)),c('Name',"Symbol","Industry","Sector")]


######Sector Downloading########
### This piece aggregates sector return and also create logical matrix pointing stock to sector
sectors_data=c()
names=c()
allocation=matrix(0,nrow = ncol(returns),ncol = length(unique(information$Sector)))
colnames(allocation)=unique(information$Sector)
rownames(allocation)=colnames(returns)

for (p in unique(information$Sector)){
  print(p)
  if(!is.na(p)){
    names=rbind(p,names)
    temp=returns[,which(colnames(returns) %in% information[which(information$Sector==p),"Symbol"])]
    allocation[which(colnames(returns) %in% information[which(information$Sector==p),"Symbol"]),which(colnames(allocation)==p)]=1
    
    if (any(is.na(temp))){
      temp=temp[, colSums(is.na(temp)) != nrow(temp)]
    }
    if(!is.null(ncol(temp))){
      temp=rowMeans(temp)
    }
    sectors_data=cbind(sectors_data,temp)
    #rm(temp)
  }
}


colnames(sectors_data)=names
date=row.names(sectors_data)
#sectors_data=cbind(sectors_data,date)


####Stocks which are in allocation table
stocks=rownames(allocation)[which(rowSums(allocation)==1)]
allocation=allocation[which(rowSums(allocation)==1),]
allocation=allocation[,which(!is.na(colnames(allocation)))]

####Stock min and max invest
max_invest=0.2
min_invest=0


#########Sector min and max invest
#colnames(allocation)
#[1] "Transportation"
#[2]"Technology"
#[3]"Health Care"
#[4]"Miscellaneous"
#[5]"Consumer Services"
#[6] NA
#[7]"Finance"
#[8]"Consumer Durables"
#[9]"Consumer Non-Durables"
#[10]"Public Utilities"     
#[11] "Energy"
#[12]"Capital Goods"
#[13]"Basic Industries" 
#             #1    #2  #3    #4  #5  #6  #7  #8  #9  #10   #11 #12   #13
sectors_max=c(0.05,0.25,0.20,0.05,0.05,0,.05,.05,0.05,0.05,0.05,0.1,0.05)
sectors_min=c(0,0,0,0,0,0,0,0,0,0,0,0,0)


################# BLACK LITTERMAN################
#View1
view1=matrix(c(0),ncol=12)
colnames(view1)=as.vector(names)
colnames(view1)
view1[which(colnames(view1)=='Finance')]=-1
view1[which(colnames(view1) %in% c('Health Care'))]=1
#view2
view2=matrix(c(0),ncol=12)
colnames(view2)=as.vector(names)
colnames(view2)
view2[which(colnames(view1)=='Consumer Non-Durables')]=1

P=rbind(view1,view2)

#######Expected Return
Q <- c(0.10,0.1)

######Temporary Uncertainty of views
Omega <- matrix(c(0.015^2,0,0,0.015^2),2,2)


get_sd_from_quantile_score <- function(m, p, x) {
  get_quantile_score <- function(y) { 
    qnorm(mean=m,sd=y,p=p)
  }
  f <- function(y) { (get_quantile_score(y) - x)^2 }
  opt <- optim(par = 1, fn = f, method = 'CG')
  return(opt$par)
}

get_sd_from_quantile_score(0.10,0.025,0.07)
0.015^2

#########Black Litterman adjustment of Expectation
var2 <- t(P) %*% ginv(Omega) %*% P
var12 <- ginv(ginv(cov(sectors_data))+var2)
var3 <- (t(P) %*% ginv(Omega) %*% Q) + (ginv(cov(sectors_data)) %*% colMeans(sectors_data)) 


###Updated means including Black Litterman
mhat <- var12 %*% var3
row.names(mhat)=colnames(sectors_data)
t(mhat)
colMeans(sectors_data)


###updated variance including Black Litterman
Blk_var=ginv(ginv(cov(sectors_data))+var2)
Blk_var
cov(sectors_data)

###############RUN OPTIMAL PORTFOLIO ON SECTOS AFTER BLACK LITTERMAN ADJUSTMENT################################

covariance=Blk_var
return=t(mhat)

max_invest=0.4
min_invest=0.01

library(quadprog)

results=c()

library(Matrix)
pd_D_mat <- nearPD(covariance)

p=0

for(i in seq (0, max(return), length.out = 100)){
  
  Dmat <- as.matrix(pd_D_mat$mat)
  dvec <- matrix(return, nrow=length(return), ncol=1)
  A.Equality <- matrix(rep(1,length(return)), ncol=1)
  
  
  ######################
  #A.Equality - all must equal to 1;
  #dvec - this is given return;
  #rep(min_invest, length(return)) - this is for minimum Stock investment
  #rep(-(max_invest), length(return)) - maximum stock investment
  ########################
  
  
  ####NO SECTOR ALLOCATION
  #Amat <- cbind(A.Equality, dvec, diag(length(return)), -diag(length(return)))
  #bvec <- c(1, i, rep(min_invest, length(return)), rep(-(max_invest), length(return)))
  
  #### SECTOR ALLOCATION
  Amat <- cbind(A.Equality, dvec, diag(length(return)), -diag(length(return)))
  bvec <- c(1, i, rep(min_invest, length(return)),rep(-max_invest, length(return)) )
  
  qp=tryCatch(solve.QP(Dmat, dvec, Amat, bvec, meq=2),error=function(e) NA)
  
  p=p+1
  if(is.na(qp)){
    print(paste('Trial',p,'return',round(i,3),'did not work.'))
  }else{print(paste('Trial',p,'return',round(i,3),'worked.'))}
  
  
  
  if(!is.na(qp)){
    result=qp$solution
    result[which(result<1e-3)]=0
    final_return=(return%*%result)
    final_variance=(result%*%covariance%*%result)
    
    results=cbind(results,c(final_return,final_variance,list(result)))
  }
}


print(results)
rownames(results)=c('return','volatility','weights')

xmin=min(unlist(results[2,]))
xmax=max(unlist(results[2,]))

ymin=min(unlist(results[1,]))
ymax=max(unlist(results[1,]))

plot(unlist(results[2,]),unlist(results[1,]),main = "Optimal Portfolios Sectors", xlim = c(xmin,xmax),ylim=c(ymin,ymax),xlab="Volatility", ylab="Return")  
par(new=TRUE)

sharpes=unlist(results[1,])/unlist(results[2,])
optimal=which(sharpes==max(sharpes))
par(new=TRUE)
plot(results[2,optimal],results[1,optimal],xlim = c(xmin,xmax),ylim=c(ymin,ymax),xlab="Volatility", ylab="Return",pch = 19, col='green')

optimal_sector_weights=unlist(results[3,optimal])

names(optimal_sector_weights)=colnames(sectors_data)
optimal_sector_weights=as.data.frame(t(optimal_sector_weights))
optimal_sector_weights=optimal_sector_weights[colnames(allocation)]
####################################################################################################################################
#AFTER WE OPTIMIZE


################################################################################################################################################


max_invest=.3
min_invest=0


portfolio_data=returns[,which(colnames(returns) %in% stocks)]
covariance=cov(portfolio_data)
return=colMeans(portfolio_data)

library(quadprog)

results=c()

library(Matrix)
pd_D_mat <- nearPD(covariance)

p=0
for(i in seq (0, max(return), length.out = 100)){

  Dmat <- as.matrix(pd_D_mat$mat)
  dvec <- matrix(return, nrow=length(return), ncol=1)
  A.Equality <- matrix(rep(1,length(return)), ncol=1)
  
  
  ######################
  #A.Equality - all must equal to 1;
  #dvec - this is given return;
  #rep(min_invest, length(return)) - this is for minimum Stock investment
  #rep(-(max_invest), length(return)) - maximum stock investment
  ########################
  
  
  ####NO SECTOR ALLOCATION
  #Amat <- cbind(A.Equality, dvec, diag(length(return)), -diag(length(return)))
  #bvec <- c(1, i, rep(min_invest, length(return)), rep(-(max_invest), length(return)))
  
  #### SECTOR ALLOCATION
  Amat <- cbind(A.Equality, dvec, diag(length(return)), -diag(length(return)),allocation)
  bvec <- c(1, i, rep(min_invest, length(return)), rep(-(max_invest), length(return)),optimal_sector_weights)
  
  qp=tryCatch(solve.QP(Dmat, dvec, Amat, bvec, meq=2),error=function(e) NA)
  
  p=p+1
  if(is.na(qp)){
    print(paste('Trial',p,'return',round(i,3),'did not work.'))
  }else{print(paste('Trial',p,'return',round(i,3),'worked.'))}
  
  
  
  if(!is.na(qp)){
    result=qp$solution
    result[which(result<1e-3)]=0
    final_return=(return%*%result)
    final_variance=(result%*%covariance%*%result)
    
    results=cbind(results,c(final_return,final_variance,list(result)))
  }
}

print(results)
rownames(results)=c('return','volatility','weights')

xmin=min(unlist(results[2,]))
xmax=max(unlist(results[2,]))

ymin=min(unlist(results[1,]))
ymax=max(unlist(results[1,]))

plot(unlist(results[2,]),unlist(results[1,]), xlim = c(xmin,xmax),ylim=c(ymin,ymax),xlab="Volatility", ylab="Return")  
par(new=TRUE)
############Min Variance#######################
min_variance=which(results[2,]==min(unlist(results[2,])))
plot(results[2,min_variance],results[1,min_variance],xlim = c(xmin,xmax),ylim=c(ymin,ymax),xlab="Volatility", ylab="Return",pch = 19, col='blue')
##############Optimal Portfolio#############
sharpes=unlist(results[1,])/unlist(results[2,])
optimal=which(sharpes==max(sharpes))
par(new=TRUE)
plot(results[2,optimal],results[1,optimal],xlim = c(xmin,xmax),ylim=c(ymin,ymax),xlab="Volatility", ylab="Return",pch = 19, col='green')


#Pick desired close return
request <- as.numeric(readline(prompt="Enter desired return: "))
picked=which(abs(unlist(results[1,])-request)==min(abs(unlist(results[1,])-request)))

par(new=TRUE)
plot(results[2,picked],results[1,picked],xlim = c(xmin,xmax),ylim=c(ymin,ymax),xlab="Volatility", ylab="Return",pch = 19, col='red')

final_weights=data.frame(weight=unlist(results[3,picked]))
rownames(final_weights)=(names(return))

#########CHECK ALLOCATION WITHIN SECTORS###############

check_allocation=matrix(0,ncol=ncol(allocation))
colnames(check_allocation)=colnames(allocation)

weights=unlist(results[3,optimal])
for (i in colnames(allocation)){
  check_allocation[which(colnames(check_allocation)==i)]=sum(weights[which(allocation[,which(colnames(allocation)==i)]==1)])
}



optimal_sector_weights
check_allocation



