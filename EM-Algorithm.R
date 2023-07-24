# ======================================================================== #
# 本程序实现了输出分布为P维正态分布的隐马尔可夫模型的EM算法，并利用Rubin在
# 其1998年的文章中提出的参数扩展方法来实现了对EM算法的加速，最后得到了估计
# 性能依旧良好，但是运行速度加快的PXEM算法
# ======================================================================== #
#library(mvtnorm)#加载一些会使用到的包
#library("mclust")
#library(MASS)

# ======================================================================== #
# 随机生成状态数为K的初始分布pi
# ======================================================================== #
initProbabilityVector=function(K) {
  probabilityVector=rep(NA, times=K)
  probabilityVectorBuilt=FALSE
  while(TRUE) {#直到ture输出
    probabilityVector[1:(K-1)] <- runif(n=K-1, min=0, max=2/K)#0-k/2的k-1个均匀分布数据
    sum <- sum(probabilityVector[1:(K-1)])#求和
    if(sum<1) {#出现一个ture的判断，如果是false重新赋值
      probabilityVector[K] <- 1-sum
      return(probabilityVector)
    }
  }
}

# ======================================================================== #
#随机生成初始转移概率矩阵
# ======================================================================== #
initA=function(K) {
  A=matrix(rep(NA, times=K^2), nrow=K)
  for(i in 1:K) {
    A[i,]=initProbabilityVector(K=K)
  }
  return(A)
}

# ======================================================================== #
#将P维数据X进行聚类，分为k类，得到初始迭代的均值向量（含k个PxP的向量的列表）
#与协方差矩阵（含k个PxP的矩阵的列表）
# ======================================================================== #
initMGauss <- function(X, nbS){
  dimObs <- ncol(X)
  clust <- mclust::Mclust(X , G = nbS)
  initMu <- t(apply(clust$parameters$mean,2,cbind))
  cluster <- clust$classification
  initSigma <- matrix(0, ncol = nbS*dimObs, nrow =dimObs)
  ystate <- list()
  for(k in 1:nbS){
    ystate[[k]] = X[cluster==k,]
    initSigma[,((k-1)*dimObs+1):(k*dimObs)] <- cov(ystate[[k]])
  }
  initSigma <- t(initSigma)
  # order for state use mu
  ordMu <- order(initMu[,1])
  # # reorder for mu
  muAllOuter <- initMu[ordMu,]
  # reorder for sigma
  sigmaAllOuter=initSigma[sapply(ordMu,function(k){((k-1)*dimObs+1):(k*dimObs)}),]
  muAllOuter <- t(muAllOuter)
  sigmaAllOuter <- t(sigmaAllOuter)
  
  
  #把上面得到的矩阵形式的参数转为list形式，因为我们是把均值向量和协方差矩阵用list进行存储
  mu <- list(NA)
  for(i in 1:ncol(muAllOuter)){
    mu[[i]] <- muAllOuter[,i]
  }
  
  sigma <- list(NA)
  for(i in 1:ncol(muAllOuter)){
    if(i == 1){
      sigma[[i]] <- sigmaAllOuter[,1:nrow(sigmaAllOuter)]
    }else{
      sigma[[i]] <- sigmaAllOuter[,((i-1)*nrow(sigmaAllOuter)+1):(i*nrow(sigmaAllOuter))]
    }
    
  }
  return(list(muAllOuter = mu,sigmaAllOuter = sigma))
}

# ======================================================================== #
#观测概率矩阵
# ======================================================================== #
EmisPr <- function(X,nbK, mu, sigma)
{ 
  nbS  = nbK#状态数
  nbT=nrow(X)#观测数
  emisPr = matrix(NA,nbK,nbT)#观测概率矩阵
  for(k in 1:nbS){
    emisPr[k,]  = dmvnorm(X, mu[[k]], sigma[[k]])
  }
  emisPr
}
# ======================================================================== #
#compute LogSum(Q) = LogSumExp(logQ)，用于解决数据下溢
# ======================================================================== #
LogSumExp <- function(logQ)
{# Q is a vector
  if(any(logQ==Inf)){
    stop("positive inifinite value in logQ \n")
  }else{
    idInf = which(logQ == -Inf)
    if(length(idInf)>0) logQ = logQ[-idInf]
  }
  if(length(logQ)==0){
    LSE = -Inf
  }else{
    maxlogQ  = max(logQ)
    LSE = maxlogQ + log(sum(exp(logQ - maxlogQ)))
  }
  LSE
}

# ======================================================================== #
# return a forward probability matrix in log scale，前向算法
# ======================================================================== #
Forward <- function(logInit, logTrans, logEmiss)
{
  nbObser = ncol(logEmiss); nbState = nrow(logEmiss)
  logF    = matrix(NA, nbState, nbObser)
  
  logF[,1] = logInit + logEmiss[,1]
  for(tt in 2:nbObser)
  {
    logF[,tt] = sapply(1:nbState, function(k){
      logEmiss[k,tt] + LogSumExp(logF[,tt-1] + logTrans[,k])
    })
  }  
  return(logF)
}

# ======================================================================== #
# return a backward probability matrix in log scale，后向算法
# ======================================================================== #
Backward <- function(logTrans, logEmiss)
{
  nbObser = ncol(logEmiss); nbState = nrow(logEmiss)
  logB    = matrix(NA, nbState, nbObser)
  
  logB[,nbObser] = 0
  for(tt in (nbObser-1):1)
  {
    logB[,tt] = sapply(1:nbState, function(k){
      LogSumExp(logB[,tt+1] + logTrans[k,] + logEmiss[,tt+1])
    })
  }
  return(logB)
} 

# ======================================================================== #
# 后验概率tau
# ======================================================================== #
postPr <- function(logF, logB){
  lognorm <- LogSumExp(logF[,ncol(logF)])
  logPosPr = logF + logB - lognorm
  exp(logPosPr)
}

# ======================================================================== #
# 计算给定模型和观测，在时刻t处于状态i且在时刻t+1处于状态j的概率
# forwardP是前向概率，A是状态转移概率矩阵，B是观测矩阵，backwardP是后向概率
# ======================================================================== #
getkeci <- function(forwardP,A,B,backwardP){
  logA = log(A)
  logB = log(B)
  nState <- nrow(forwardP)#状态数
  nObser <- ncol(forwardP)-1#观测数-1,公式里是到T-1，不是到T，要注意
  keci <- array(NA,c(nState,nState,nObser))
  #先算keci的分子部分
  for(tt in 1:nObser){
    for (ii in 1:nState) {
      for(jj in 1:nState){
        keci[ii,jj,tt] <- forwardP[ii,tt]+logA[ii,jj]+logB[jj,tt+1]+backwardP[jj,tt+1]
      }
      
    }
    
  }
  logNorm = LogSumExp(forwardP[,ncol(forwardP)])
  #分子/分母，注意取了对数，所以相减
  for(tt in 1:nObser){
    for (ii in 1:nState) {
      for(jj in 1:nState){
        keci[ii,jj,tt] <- keci[ii,jj,tt] - logNorm
      }
    }
  }
  keci <- exp(keci)#取e，转换回去
  return(keci)
}

# ======================================================================== #
# 维特比算法
# ======================================================================== #
Viterbi <- function(initPr, transPr,emisPr, mu, sigma)
{
  nbT = ncol(emisPr)#观测数
  nbS = nrow(emisPr)#状态数
  logF = matrix(0, ncol=nbT,nrow=nbS)
  logF[,1] = log(initPr) + log(emisPr[,1])
  
  for (t in 2:nbT) {
    for (k in 1:nbS) {
      logF[k,t] = log(emisPr[k,t]) +  max(logF[,t-1] + log(transPr[,k]))
    }
  }
  
  path = rep(NA, nbT)
  path[nbT] = which.max(logF[,nbT])
  
  for (t in (nbT-1):1) {
    path[t] = which.max(logF[,t] + log(transPr[,path[t+1]]))
  }
  return(path)
}

# ======================================================================== #
# EM算法
# ======================================================================== #
EM <- function(y, nbK,iterMax=1000,epsilon=1e-3){
  #P维高斯分布，输入的y是一个P列矩阵，行数是样本总数
  #nbK是待估计高斯分布的数量
  dimy <- ncol(y)#维数P
  N = nrow(y)#观测数，即样本数
  chushi=initMGauss(y,nbK)#先对数据进行一步聚类
  mu = chushi$muAllOuter #初始迭代的均值向量
  sigma = chushi$sigmaAllOuter #初始迭代的协方差矩阵
  pi=c(rep(1/nbK,nbK))#初始迭代的初始分布，取为均匀分布
  a <- 0.01/(nbK-1)
  A = matrix(c(rep(a,nbK)),nbK,nbK,byrow = T)#初始迭代的转移概率矩阵，取对角线上为0.99，其余为均匀
  diag(A) = c(rep(0.99,nbK))
  for(t in 1:iterMax){
    #旧参数，用于最后的收敛判断
    osumP1 <- matrix(0,dimy,nbK)
    for(z in 1:nbK){
      osumP1[,z] <- mu[[z]]
    }
    for(z in 1:nbK){
      osumP1 <- cbind(osumP1,sigma[[z]])
    }
    osumP2 <- cbind(pi,A)
    # E 步
    B=EmisPr(y,nbK,mu,sigma)#观测概率矩阵
    logB = log(B)#取对数避免数据下溢
    logA = log(A)
    logpi = log(pi)
    forwardP <- Forward(logpi,logA,logB)
    backwardP <- Backward(logA,logB)
    tau = t(postPr(forwardP,backwardP))#后验概率矩阵
    # M 步
    for(k in 1:nbK){
      a <- c(rep(0,dimy))
      b <- matrix(0,dimy,dimy)
      for(i in 1:N){
        a <- a + tau[i,k]*y[i,]
      }
      mu[[k]] = a/sum(tau[,k])
      muk <- mu[[k]]
      for(i in 1:N){
        b <- b + tau[i,k]*(t(t(y[i,]-muk))%*%t(y[i,]-muk))
      }
      sigma[[k]] = b/sum(tau[,k])
    }#更新均值向量与协方差矩阵
    
    keci <- getkeci(forwardP,A,B,backwardP)#计算给定模型和观测，在时刻t处于状态i且在时刻t+1处于状态j的概率
    
    for (ii in 1:nbK) {
      for(jj in 1:nbK){
        fenzi = c(rep(NA,N-2))
        for(tt in 2:N-1){
          fenzi[tt] = keci[ii,jj,tt]
        }
        A[ii,jj] <- LogSumExp(log(fenzi)) - LogSumExp(log(tau[1:(N-2),ii]))
      }
    }#更新状态转移概率矩阵
    A <- exp(A)#注意要把取log转换回去
    pi <- colMeans(tau)#初始分布这里一直有问题，困扰了很久，如果是按公式取为tau[1,]，就会出现某一个状态是1，其余为0，明显不合理
    #取列均值的话可以得到比较正常的数字，但是经过检验和真实的值差的比较大
    #pi <- tau[1,]
    
    # 迭代停止
    #更新后的参数
    csumP1 <- matrix(0,dimy,nbK)
    for(z in 1:nbK){
      csumP1[,z] <- mu[[z]]
    }
    for(z in 1:nbK){
      csumP1 <- cbind(csumP1,sigma[[z]])
    }
    csumP2 <- cbind(pi,A)
    sum1 = sum(abs(osumP1 - csumP1))
    sum2 = sum(abs(osumP2 - csumP2))
    if(sum1<epsilon & sum2<epsilon){break}
  }
  
  parameter <- list(mu = mu, sigma = sigma, Init = pi, Trans = A)
  
  return(list(parameter=parameter,iter = t, tau=tau,data = y))
}


# ======================================================================== #
# 生成初始训练数据y
# ======================================================================== #
nState = 3#设定状态数
pi=initProbabilityVector(nState);pi#生成数据的初始分布
A <- initA(nState);A#生成数据的状态转移概率矩阵
mu_set <- list(c(-3,-2,-3),c(-2,-1,0),c(1,0,1))#生成数据的nState个总体的均值向量，以列表形式存储，列表的每个位置是个dim行1列的向量，dim是数据维数
sigma_set <- list(matrix(c(1,0,0,0,1,0,0,0,1),c(3,3),byrow = T),
                  matrix(c(0.4,0,0,0,0.4,0,0,0,0.4),c(3,3),byrow = T),
                  matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1),c(3,3),byrow = T))#生成数据的nState个总体的协方差矩阵，以列表形式存储，列表的每个位置是个dim*dim的矩阵
N = 1000#观测数
y <- matrix(NA,N,3)#空矩阵，用来记录生成数据，列数是数据的维数
oripath <- rep(NA,N)#用来记录生成数据的马氏链
for(t in 1:N){
  alpha_set <- pi
  #1.抽子总体
  for(i in 2:t){
    alpha_set <- alpha_set%*%A
  }
  sp <- sample(1:3, 1, alpha_set, replace=TRUE)
  oripath[t] <- sp
  meansp <- mu_set[[sp]]
  sigmasp <- sigma_set[[sp]]
  #2.生成随机数
  y[t,] <- mvrnorm(1, meansp, sigmasp)
}


# ======================================================================== #
# EM算法测试
# ======================================================================== #
system.time(res <- EM(y,nState))#对上面我们自己生成的数据y进行EM算法，会一并给出运行时间
res$parameter#估计出来的参数
res$iter#收敛迭代次数

data <- load("/Users/tonnytang/Desktop/simObs.rdata")#读入测试数据
simObs_1 <- as.matrix(simObs_1)#我们是把初始数据当成矩阵输入到EM算法中，所以要把数据转成矩阵
simObs_2 <- as.matrix(simObs_2)


system.time(res1 <- EM(y = simObs_1,nbK = 3))#对老师给的两组测试数据都做EM算法,nbK是马氏链状态数
res1$parameter#估计出来的参数
res1$iter#收敛迭代次数

system.time(res2 <- EM(y = simObs_2,nbK = 3))
res2$parameter
res2$iter
# ======================================================================== #
# 维特比算法测试模型状态分类正确率
# ======================================================================== #
curres <- res1#当前对哪个数据集进行状态分类正确率判断，可以取res,res2,PXres1,PXres2,这里PXres1和2都是下面利用加速后的EM算法PXEM得到的结果
initPr <- (curres$parameter)$Init#对当前数据集估计出来的初始分布
transPr <- (curres$parameter)$Trans#对当前数据集估计出来的状态转移概率矩阵
mu <- (curres$parameter)$mu#对当前数据集估计出来的均值向量
sigma <- (curres$parameter)$sigma#对当前数据集估计出来的协方差矩阵
emisPr <- EmisPr(curres$data,length(initPr),mu,sigma)#构建维特比算法所需要的参数
V <- Viterbi(initPr,transPr,emisPr,mu,sigma)#维特比算法

load("/Users/tonnytang/Desktop/simdata_paras.rdata")#读入当前训练集的真实参数
path <- simdata_1$simHMM[,2] - V#这里也可以取simdata_2，把生成原始数据的马氏链的状态跳转情况记为path
#path <- oripath - V #这个是对我们自己生产的测试数据来做判断
condition = 0#用于记录状态被正确判断的观测数
for(i in 1:length(path)){
  if(path[i]==0){
    condition = condition + 1
  }
};condition

# ======================================================================== #
# 以下利用参数扩展思想来加速EM算法
# 利用参数扩展的思想修改前向算法
# ======================================================================== #
PXForward <- function(logInit, logTrans, logEmiss,alpha)
{
  nbObser = ncol(logEmiss); nbState = nrow(logEmiss)
  logF    = matrix(NA, nbState, nbObser)
  
  logF[,1] = logInit + logEmiss[,1]
  for(tt in 2:nbObser)
  {
    logF[,tt] = sapply(1:nbState, function(k){
      logEmiss[k,tt] + LogSumExp(logF[,tt-1] + logTrans[,k])
    })
  }  
  return(logF + log(alpha))
}

# ======================================================================== #
# PXEM算法，利用参数扩展思想加速EM算法
# ======================================================================== #
PXEM <- function(y, nbK,iterMax=1000,epsilon=1e-3){
  #P维高斯分布，输入的y是一个P列矩阵，行数是样本总数
  #nbK是待估计高斯分布的数量
  dimy <- ncol(y)#维数P
  N = nrow(y)#观测数，即样本数
  chushi=initMGauss(y,nbK)#先对数据进行一步聚类
  mu = chushi$muAllOuter #初始迭代的均值向量
  sigma = chushi$sigmaAllOuter #初始迭代的协方差矩阵
  pi=c(rep(1/nbK,nbK))#初始迭代的初始分布，取为均匀分布
  a <- 0.01/(nbK-1)
  A = matrix(c(rep(a,nbK)),nbK,nbK,byrow = T)#初始迭代的转移概率矩阵，取对角线上为0.99，其余为均匀
  diag(A) = c(rep(0.99,nbK))
  alpha = 1
  for(t in 1:iterMax){
    #旧参数，用于最后的收敛判断
    osumP1 <- matrix(0,dimy,nbK)
    for(z in 1:nbK){
      osumP1[,z] <- mu[[z]]
    }
    for(z in 1:nbK){
      osumP1 <- cbind(osumP1,sigma[[z]])
    }
    osumP2 <- cbind(pi,A)
    # E 步
    B=EmisPr(y,nbK,mu,sigma)#观测概率矩阵
    logB = log(B)#取对数避免数据下溢
    logA = log(A)
    logpi = log(pi)
    forwardP <- PXForward(logpi,logA,logB,alpha)
    backwardP <- Backward(logA,logB)
    tau = t(postPr(forwardP,backwardP))#后验概率矩阵
    # M 步
    for(k in 1:nbK){
      a <- c(rep(0,dimy))
      b <- matrix(0,dimy,dimy)
      for(i in 1:N){
        a <- a + tau[i,k]*y[i,]
      }
      mu[[k]] = a/sum(tau[,k])
      muk <- mu[[k]]
      for(i in 1:N){
        b <- b + tau[i,k]*(t(t(y[i,]-muk))%*%t(y[i,]-muk))
      }
      sigma[[k]] = b/sum(tau[,k])
    }#更新均值向量与协方差矩阵
    
    keci <- getkeci(forwardP,A,B,backwardP)#计算给定模型和观测，在时刻t处于状态i且在时刻t+1处于状态j的概率
    
    for (ii in 1:nbK) {
      for(jj in 1:nbK){
        fenzi = c(rep(NA,N-2))
        for(tt in 2:N-1){
          fenzi[tt] = keci[ii,jj,tt]
        }
        A[ii,jj] <- LogSumExp(log(fenzi)) - LogSumExp(log(tau[1:(N-2),ii]))
      }
    }#更新状态转移概率矩阵
    A <- exp(A)#注意要把取log转换回去
    alpha <- sum(A[1,])
    A <- A/alpha
    pi <- colMeans(tau)#初始分布这里一直有问题，困扰了很久，如果是按公式取为tau[1,]，就会出现某一个状态是1，其余为0，明显不合理
    #取列均值的话可以得到比较正常的数字，但是经过检验和真实的值差的比较大
    #pi <- tau[1,]
    
    # 迭代停止
    #更新后的参数
    csumP1 <- matrix(0,dimy,nbK)
    for(z in 1:nbK){
      csumP1[,z] <- mu[[z]]
    }
    for(z in 1:nbK){
      csumP1 <- cbind(csumP1,sigma[[z]])
    }
    csumP2 <- cbind(pi,A)
    sum1 = sum(abs(osumP1 - csumP1))
    sum2 = sum(abs(osumP2 - csumP2))
    if(sum1<epsilon & sum2<epsilon){break}
  }
  
  parameter <- list(mu = mu, sigma = sigma, Init = pi, Trans = A)
  
  return(list(parameter=parameter,iter = t, tau=tau,data = y))
}

# ======================================================================== #
# PXEM算法测试
# ======================================================================== #
system.time(PXres1 <- PXEM(y = simObs_1,nbK = 3))
PXres1$parameter#加速后的估计参数
PXres1$iter#加速后的迭代次数

system.time(PXres2 <- PXEM(y = simObs_2,nbK = 3))
PXres2$parameter#加速后的估计参数
PXres2$iter#加速后的迭代次数
#可以将PXres1和PXres2拿到上面的维特比算法中进行检验