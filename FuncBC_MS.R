######################################################################################
# Purpose:       Copula-based Functional Bayes Classification on Multiple Sclerosis  #
#                                                                                    #
#                                                                                    #
# Author:        Wentian Huang                                                       #
# Contact:       wh365@cornell.edu                                                   #
#                                                                                    #
# Last updated:  2020-07-24                                                          #
#                                                                                    #
# Acknowledgement: the MRI/DTI data were collected at Johns Hopkins University       #
#                  and the Kennedy-Krieger Institute                                 #
######################################################################################


################################## Load Pkgs ###################################


Packages <- c("pls", "fda", "fda.usc", "refund", "Matrix", "MASS", "parallel", "copula", "pracma", 
              "KernSmooth", "caret")

lapply(Packages, library, character.only = TRUE)






################################################################################
#                                                                              
#                CLASSIFIERS on Fully Observed Functional Data
#
# Note: all functional data used here should be pre-smoothed 
# Arguments: 
#   - C_train: score matrix of training functions on Jtotal basis. Jtotal > CV candidate J values
#   - Y_train: response vector of training
#   - c_test: score vector / matrix of test functions on Jtotal basis
#   - J: cut-off dimention value
#   - h0, h1: length J bandwidth vectors for KDE of group 0, 1; scores on PC
#   - hh0, hh1: bandwidth vectors for KDE of group 0, 1; scores on PLS
#   - copula_type: type of parameteric copulas to use: Gaussian or t
################################################################################

########################## Kernel Density Est. w. Gaussian###########################
myKDE <- function(xs,t,h){
  kernelValues <- rep(0,length(xs))
  for(i in 1:length(xs)){
    transformed = (t - xs[i]) / h
    kernelValues[i] <- dnorm(transformed, mean = 0, sd = 1) / h
  }
  return(sum(kernelValues) / length(xs))
}

########################### General Bayes Classifier ###########################

BC<-function(C_train, Y_train, c_test, J, h0, h1){
  C_train = C_train[, 1:J]
  c_test = c_test[, 1:J]
  
  m1 = as.matrix(C_train[Y_train == 1,]) #training data in group 1
  m0 = as.matrix(C_train[Y_train == 0,]) #training data in group 0
  
  vv.R = rep(4,nrow(c_test)) #vector to store estimated group label for each test case. 
  #4 is used as original entry value since it would be easier to debug if returned value is not 0 or 1
  
  for(i in 1:nrow(x.tt)){
    test.k1 = rep(4,J) # vector to store individual densities
    test.k0 = rep(4,J)
    for(j in 1:J){ #kernel density estimation on each basis for test case i in c_test
      f1x = myKDE(as.matrix(m1)[,j], c_test[i,j], h1[j])
      f0x = myKDE(as.matrix(m0)[,j], c_test[i,j], h0[j])
      test.k1[j] = f1x
      test.k0[j] = f0x
    }
    
    #below is for when densities estimated above for each basis j are so small that the KDE function returns 0.
    #To avoid the classifier returning NA values after calculating density ratios f1j/f0j: f1j as j-th score at group 1
    
    #when for observation i, there exists both f1j=0 and f0j=0 for some j, observation i is labeled -1, i.e. always 
    #counted as misclassified case;
    
    #when for observation i, only group 1 has f1j=0 for some j, it is classified to group 0;
    
    #when for observation i, only group 0 has f0j=0 for some j, it is classified to group 1;
    
    #otherwise density ratio is estimated to decide the group label.
    
    if(sum(test.k1 == 0) > 0 & sum( test.k0 == 0) > 0){
      vv.R[i] = -1
    }else if(sum(test.k1 == 0) >0 & sum(test.k0 == 0) == 0){
      vv.R[i] = 0
    }else if(sum(test.k1 == 0) == 0 & sum(test.k0 == 0) > 0){
      vv.R[i] = 1
    }else{
      ratio.ind = sum(log(test.k1 / test.k0)) + log(sum(Y_train) / (length(Y_train) - sum(Y_train)))
      vv.R[i] = as.numeric(ratio.ind > 0)
    }
  }
  
  return(vv.R) #return estimated class label
}




#################### Bayes Classifier using Copulas ####################

BC.copula<-function(C_train_PC, C_train_PLS, Y_train, c_test_pc, c_test_pls, J, h0, h1, hh0, hh1, score_type, copula_type){
  if(score_type == "PC"){
    C_train = C_train_PC
    c_test = c_test_pc
  }else{
    C_train = C_train_PLS
    c_test = c_test_pls
  }
  C_train = C_train[, 1:J]
  c_test = c_test[, 1:J]
  
  m1 = as.matrix(C_train[Y.t == 1,]) #group 1 training cases
  m0 = as.matrix(C_train[Y.t == 0,]) #group 0 training cases
  
  U1 = matrix(rep(0, nrow(as.matrix(m1)) * J), ncol = J) #matrix to store estimated CDF of corresponding score
  U0 = matrix(rep(0, nrow(as.matrix(m0)) * J), ncol = J)
  
  for(j in 1:J){
    U1[,j] = match(m1[,j], sort(m1[,j])) / (nrow(m1) + 1) #estimating CDF for copula fitting
    U0[,j] = match(m0[,j], sort(m0[,j])) / (nrow(m0) + 1)
  }
  
  #similarly, generating estimated CDF of scores from 1 to J for test cases
  U1.test = matrix(rep(0,nrow(c_test) * J), ncol = J)
  U0.test = matrix(rep(0,nrow(c_test) * J), ncol = J)
  
  for(i in 1:nrow(c_test)){
    u1t = apply(t(m1), 2 ,'<=', as.vector(c_test[i,]))
    U1.test[i,] = (rowSums(u1t) + 1) / (nrow(m1) + 2)
    
    
    u0t = apply(t(m0), 2, '<=', as.vector(c_test[i,]))
    U0.test[i,] = (rowSums(u0t) + 1) / (nrow(m0) + 2)
  }
  
  if(copula_type == "Gaussian"){
    
    #generating Gaussian copula densities for group 1 and 0
    fit1.cop <- fitCopula(normalCopula(dim = J, dispstr = "un"), U1, method = "itau")
    cop1 = normalCopula(param = unname(coef(fit1.cop)), dim = J, dispstr = "un")
    
    fit0.cop <- fitCopula(normalCopula(dim = J, dispstr = "un"), U0, method = "itau")
    cop0 = normalCopula(param = unname(coef(fit0.cop)), dim = J, dispstr = "un")
    
    #estimated copula densities of group 1 and 0 for all observations in test set
    c1.mvalue = dCopula(U1.test, cop1)
    c0.mvalue = dCopula(U0.test, cop0)
    
  }else{
    
    fit1.cop <- fitCopula(tCopula(dim = J, dispstr = "un"), U1, method = "itau.mpl")
    coe1.v = unname(coef(fit1.cop))
    cop1 = tCopula(param = coe1.v[1:(length(coe1.v) - 1)], dim = J, dispstr = "un", df = coe1.v[length(coe1.v)])
    
    fit0.cop <- fitCopula(tCopula(dim = J, dispstr = "un"), U0, method = "itau.mpl")
    coe0.v = unname(coef(fit0.cop))
    cop0 = tCopula(param = coe0.v[1:(length(coe0.v) - 1)], dim = J, dispstr = "un", df = coe0.v[length(coe0.v)])
    
    c1.mvalue = dCopula(U1.test, cop1)
    c0.mvalue = dCopula(U0.test, cop0)
    
  }
  
  vv.R = rep(4, nrow(c_test))
  
  for(i in 1:nrow(c_test)){
    test.k1 = rep(4, J + 1) #vector storing J estimated marginal densities, and estimated copula density
    test.k0 = rep(4, J + 1)
    
    for(j in 1:J){
      f1x = myKDE(as.matrix(m1)[,j], c_test[i,j], h1[j])
      f0x = myKDE(as.matrix(m0)[,j], c_test[i,j], h0[j])
      test.k1[j] = f1x
      test.k0[j] = f0x
    }
    test.k1[J + 1] = c1.mvalue[i]
    test.k0[J + 1] = c0.mvalue[i]
    
    if(sum(test.k1 == 0) > 0 & sum(test.k0 == 0)>0){
      vv.R[i] = -1
    }else if(sum(test.k1 == 0) > 0 & sum(test.k0 == 0) == 0){
      vv.R[i] = 0
    }else if(sum(test.k1 == 0) == 0 & sum(test.k0 == 0) > 0){
      vv.R[i] = 1
    }else{
      ratio.ind = sum(log(test.k1 / test.k0)) + log(sum(Y.t) / (length(Y.t) - sum(Y.t)))
      vv.R[i] = as.numeric(ratio.ind > 0)
    }
  }
  
  return(vv.R)
}



######################## Functional Centroid Classifier ########################

CenC<-function(X.t, Y.t, x.tt, r){ #r is number of basis chosen 
  X.c = as.matrix(sweep(X.t, 2, colMeans(X.t))) #centered X
  S.R = svd(X.c)
  basis.m = S.R$v[, 1 : r] #discrete r eigenvectors of design matrix X
  eig.v = (S.R$d[1 : r]) ^ 2 / (nrow(X.c) * ptv) #true functional eigenvalue list
  fscore.m = (X.t %*% basis.m) * (1 / sqrt(ptv)) #conversion to actual functional scores from discrete scores
  umean = colMeans(as.matrix(fscore.m[Y.t == 1,])) - colMeans(as.matrix(fscore.m[Y.t == 0,])) #first r mean scores
  phi.coef = (1 / eig.v) * umean  
  
  #Centroid classifier projects test function case x to function phi.tk as built below. phi.tk evaluated at tlist as mentioned before
  if(r == 1){
    phi.tk = basis.m * phi.coef * sqrt(ptv)
  }else{
    phi.tk = basis.m %*% phi.coef * sqrt(ptv)
  }
  
  G1.mean = colMeans(X.t[Y.t == 1,]) #function mean at each group
  G0.mean = colMeans(X.t[Y.t == 0,])
  
  class.v = rep(0,nrow(x.tt)) #vector to store classification result for each test case
  
  for(i in 1 : nrow(x.tt)){
    TX1 = as.numeric(((x.tt[i,] - G1.mean) %*% phi.tk / ptv) ^ 2)
    TX0 = as.numeric(((x.tt[i,] - G0.mean) %*% phi.tk / ptv) ^ 2)
    class.v[i] = as.numeric(TX1 < TX0) #comparing test case x's distance to each group mean to decide label
  }
  
  return(class.v)
}



#################################### PLSDA #####################################

PLS.CEN<-function(X.t, Y.t, x.tt, r){ 
  Y.f = factor(Y.t)
  X.t = X.t * (1 / sqrt(ptv))
  model.pls = caret:::plsda(X.t, Y.f, ncomp = r)
  
  class.v = predict(model.pls, x.tt, type = "class")
  if(is.factor(class.v)){
    class.v = as.numeric(levels(class.v)[class.v])
  }
  
  return(class.v)
}




######################## Functional Logistic Regression ########################

logistic.FC<-function(C_train, Y_train, c_test, J){
  C_train = C_train[, 1:J]
  c_test = c_test[, 1:J]
  Y.f = factor(Y_train)
  dd.fr = data.frame("group" = Y.f,"fscore" = C_train)
  
  model1.glm = glm(group~., family = binomial(link = 'logit'), data = dd.fr)
  fitted.results <- predict(model1.glm, newdata = data.frame("fscore" = c_test), type = 'response')
  pre.class=as.numeric(fitted.results > 0.5)
  return(pre.class)
}



################################################################################
#                                                                              
#            Cross Validation of Optimal Method and Hyperparameters            
# Arguments:
#   - pc_h0, pc_h1, pls_h0, pls_h1: bandwidths for scores at group 0 and 1, based on PC and PLS
#   - J.t: ceiling J value to consider
################################################################################

BCCV<-function(X.t, Y.t, pc_h0, pc_h1, pls_h0, pls_h1, J.t){ #BCCV is for one data scenario
  folds <- cut(seq(1, nrow(X.t)), breaks = 10, labels = FALSE)
  named_fun1 = rep(c("PC", "PLS"), 2)
  named_fun2 = c(rep("Gaussian", 2), rep("t", 2))
  rate.estimators = matrix(rep(0, 8 * J.t), 8) #update when number of classifiers changed
  
  for(k in 1:10){
    Xtest = X.t[folds == k,]
    Ytest = Y.t[folds == k]
    Xtrain = X.t[folds != k,]
    Ytrain = Y.t[folds != k]
    
    #decompose PC scores
    X_train = as.matrix(sweep(Xtrain, 2, colMeans(Xtrain))) #centered X
    basis.m = svd(X_train)$v[, 1:J.t] #discrete J.t eigenvectors of design matrix X
    Ctrain_PC = (Xtrain %*% basis.m) * (1 / sqrt(ptv)) #conversion to actual functional scores from discrete scores
    Ctest_PC = (Xtest %*% basis.m) * (1 / sqrt(ptv)) #actual test cases functional scores
    
    #decompose PLS scores
    train.set = data.frame(Ytrain, Xtrain)
    pls.model = cppls(Ytrain ~., J.t, data = train.set, validation="none")
    pls.add.mean = as.vector(matrix(pls.model$Xmeans, nrow = 1) %*% pls.model$projection)
    dscore.m = as.matrix(sweep(pls.model$scores, 2, pls.add.mean,"+"))  #PLS first J scores
    Ctrain_PLS = dscore.m * (1 / sqrt(ptv)) #scaled to actual functional scores
    
    dtest.score = Xtest %*% pls.model$projection
    Ctest_PLS = dtest.score * (1 / sqrt(ptv))
    
    #when j=1
    Yest.BC.k = BC(Ctrain_PC, Ytrain, Ctest_PC, 1, pc_h0, pc_h1)
    rate.estimators[1,1] = rate.estimators[1,1] + sum(Yest.BC.k!=Ytest) / nrow(X.t)
    
    rate.estimators[2:5,1]=10 # value 10 as a sign of wrong output. Classifiers 2 ~ 5 are copula based. J >= 2
    
    Yest.Cen.k = CenC(Xtrain, Ytrain, Xtest, 1)
    rate.estimators[6,1] = rate.estimators[6,1] + sum(Yest.Cen.k != Ytest) / nrow(X.t)
    
    Yest.PLSDA.k = PLS.CEN(Xtrain, Ytrain, Xtest, 1)
    rate.estimators[7,1] = rate.estimators[7,1] + sum(Yest.PLSDA.k != Ytest) / nrow(X.t)
    
    Yest.logistic.k = logistic.FC(Ctrain_PC, Ytrain, Ctest_PC, 1)
    rate.estimators[8,1] = rate.estimators[8,1] + sum(Yest.logistic.k != Ytest) / nrow(X.t)
    
    for(j in 2:J.t){
      BC_j_est = BC(Ctrain_PC, Ytrain, Ctest_PC, j, pc_h0, pc_h1)
      rate.estimators[1, j] = rate.estimators[1, j] + sum(BC_j_est != Ytest) / nrow(X.t)
      
      copula_res = mapply(BC.copula, score_type = named_fun1, copula_type = named_fun2, 
                          MoreArgs = list(C_train_PC = Ctrain_PC, C_train_PLS = Ctrain_PLS, Y_train = Ytrain,
                                          c_test_pc = Ctest_PC, c_test_pls = Ctest_PLS, J = j, h0 = pc_h0,
                                          h1 = pc_h1, hh0 = pls_h0, hh1 = pls_h1), SIMPLIFY = FALSE)
      
      rate.estimators[2, j] = rate.estimators[2, j] + sum(copula_res[[1]] != Ytest) / nrow(X.t)
      rate.estimators[3, j] = rate.estimators[3, j] + sum(copula_res[[2]] != Ytest) / nrow(X.t)
      rate.estimators[4, j] = rate.estimators[4, j] + sum(copula_res[[3]] != Ytest) / nrow(X.t)
      rate.estimators[5, j] = rate.estimators[5, j] + sum(copula_res[[4]] != Ytest) / nrow(X.t)
      
      Cen_j_est = CenC(Xtrain, Ytrain, Xtest, j)
      rate.estimators[6, j] = rate.estimators[6, j] + sum(Cen_j_est != Ytest) / nrow(X.t)
      
      PLSDA_j_est = PLS.CEN(Xtrain, Ytrain, Xtest, j)
      rate.estimators[7, j] = rate.estimators[7, j] + sum(PLSDA_j_est != Ytest) / nrow(X.t)
      
      logistic_j_est = logistic.FC(Ctrain_PC, Ytrain, Ctest_PC, j)
      rate.estimators[8, j] = rate.estimators[8, j] + sum(logistic_j_est != Ytest) / nrow(X.t)
      
    }
  }
  
  BC.J=which(as.vector(rate.estimators[1,])==min(as.vector(rate.estimators[1,])))[1]
  BCG.J=which(as.vector(rate.estimators[2,])==min(as.vector(rate.estimators[2,])))[1]
  BCGPLS.J=which(as.vector(rate.estimators[3,])==min(as.vector(rate.estimators[3,])))[1]
  BCt.J=which(as.vector(rate.estimators[4,])==min(as.vector(rate.estimators[4,])))[1]
  BCtPLS.J=which(as.vector(rate.estimators[5,])==min(as.vector(rate.estimators[5,])))[1]
  Cen.J=which(as.vector(rate.estimators[6,])==min(as.vector(rate.estimators[6,])))[1]
  PLSDA.J=which(as.vector(rate.estimators[7,])==min(as.vector(rate.estimators[7,])))[1]
  logistic.J=which(as.vector(rate.estimators[8,])==min(as.vector(rate.estimators[8,])))[1]
  
  
  min.indicator = which(rate.estimators == min(rate.estimators), arr.ind = TRUE)
  min.row = min.indicator[, 1]
  
  return(list(J.BC = BC.J, J.BCG = BCG.J, J.BCGPLS = BCGPLS.J, J.BCt = BCt.J, J.BCtPLS = BCtPLS.J,
              J.CEN = Cen.J, J.PLSDA = PLSDA.J, J.logistic = logistic.J, CV.selected = unname(min.row[1]),
              JCVmatrix = rate.estimators))
  
}



############################ PC bandwidth generator ############################

bd.gen <- function(X.train, Y.train, J.t){
  X.tc = as.matrix(sweep(X.train, 2, colMeans(X.train))) #centered X
  basis.m = svd(X.tc)$v[, 1:J.t] #discrete J eigenvectors of design matrix X
  fscore.m = (X.train %*% basis.m) / sqrt(ptv) #conversion to actual functional scores from discrete scores
  m1 = fscore.m[Y.train == 1,] #group 1 training cases
  m0 = fscore.m[Y.train == 0,] #group 0 training cases
  
  h1.v = rep(1, J.t)
  h0.v = rep(1, J.t)
  
  for(j in 1:J.t){
    h1.v[j] = dpik(as.vector(m1[, j])) #direct plug-in methodology
    h0.v[j] = dpik(as.vector(m0[, j]))
  }
  
  return(list(h1.list = h1.v, h0.list = h0.v))
}



########################### PLS bandwidth generator ############################

bdpls.gen <- function(X.train, Y.train, J.t){ 
  train.set = data.frame(Y.train, X.train)
  pls.model = cppls(Y.train~., J.t, data = train.set, validation="none")
  score.m = pls.model$scores #PLS first J.t scores
  
  m1 = score.m[Y.train == 1,] #group 1 training cases
  m0 = score.m[Y.train == 0,] #group 0 training cases
  
  h1.v = rep(1, J.t)
  h0.v = rep(1, J.t)
  
  for(j in 1:J.t){
    h1.v[j] = dpik(as.vector(m1[, j])) #direct plug-in methodology
    h0.v[j] = dpik(as.vector(m0[, j]))
  }
  
  return(list(h1.list = h1.v, h0.list = h0.v))
}


################################################################################
#                                                                              
#                               Data Preparation                               
#                                                                              
################################################################################
data(DTI)
DTI.visit1 = DTI[DTI$visit == 1,] #using only visit = 1 cases
DTI.visit1.clean = DTI.visit1[-59,] #removing the case with missing data
cca.m = DTI.visit1.clean$cca #X matrix of functional data
cca.case = DTI.visit1.clean$case #total Y vec

cca.smoothed = matrix(rep(0, nrow(cca.m) * ncol(cca.m)), nrow = nrow(cca.m)) #matrix to store smoothed data

#smoothing each observation with local linear regression
for(i in 1:nrow(cca.m)){
  cca.smooth.h = dpill(1:93, as.vector(cca.m[i,]), gridsize = 93)
  cca.smoothed[i,] = locpoly(1:93, as.vector(cca.m[i,]), bandwidth = cca.smooth.h, gridsize = 93)$y
}

Jtotal = 30
#generate overall bandwidth vectors for group 1 and 0, keeping 20 functional basis at maximum
H.list = bd.gen(cca.smoothed, cca.case, Jtotal)
H.list.pls = bdpls.gen(cca.smoothed, cca.case, Jtotal)

tlist = 1:93 #argument t values for functional data
pt = length(tlist)
ptv = 1 #inverse of bin size to scale functional scores


################################################################################
#                                                                              
#                   Wrapper Function for Parallel Computing                    
# Returned Values:
#   - rate.CV.vec: error rate vector of each classifier using cross validated J
#   - BCCV.list: selected optimal classifier for the dataset
#   - J.CV.vec: cross validated J of each method
#   - J.trend: error rates of each classifier by J values ranged from 1 to Jtotal
################################################################################

wrapper.cca=function(iter){
  #shuffle the dataset
  order.cca = sample(nrow(cca.m))
  cca.m.shuffled = cca.smoothed[order.cca,]
  cca.case.shuffled = cca.case[order.cca]
  rate.CV.vec = rep(-1, length(classifier.list))
  
  
  BCCV.data = BCCV(cca.m.shuffled, cca.case.shuffled, H.list$h0.list, H.list$h1.list, H.list.pls$h0.list,
                 H.list.pls$h1.list, Jtotal)
  BCCV.result = unlist(BCCV.data[1:9])
  names(BCCV.result) = NULL
  BCCV.list = BCCV.result[length(classifier.list) + 1] #selected classifier index
  J.CV.vec = BCCV.result[1:length(classifier.list)]
  
  J.trend.matrix = BCCV.data$JCVmatrix
  J.trend = as.vector(t(J.trend.matrix))
  
  for(k in 1:8){
    rate.CV.vec[k] = J.trend.matrix[k, J.CV.vec[k]]
  }
  
  
  return(list(c(rate.CV.vec, BCCV.list, J.CV.vec), J.trend))
  
}



############################## Parallel Computing ##############################

library(parallel)
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores,type = "FORK")
clusterSetRNGStream(cl,iseed = 961210)

tic()
CL.result=parLapply(cl,1:1000,wrapper.cca) #simulate 1000 times
CL.total=do.call(rbind,CL.result)
CL.CVfinal=do.call(rbind,CL.total[,1])
CL.Jtrend=do.call(rbind,CL.total[,2])
stopCluster(cl)
toc()


