rm(list=ls())
getwd()
library(MASS)
data<-Boston
dim(data)
n=nrow(data)
p=ncol(data)
'p is Number of unknown coefficient to be estimated for regression line '
Y=(data$medv)                      #Response
X=as.matrix(cbind(rep(1,n),data[,-which(colnames(data)=='medv')]))    #Design Matrix or Predictor matrix 
colnames(X)[1]='intercept'
#determinant of (t(X)%*%X) is not zero . Hence we conclude the matrix is full column rank matrix
beta_hat=solve(t(X)%*%X)%*%t(X)%*%Y         #Least square solution

#------------Testing Significance of Regressors for predicting mean hosing price-------------------------------------------------------------

#Ho=beta(i)=0 for all i=1(1)13 vs H1:- beta(i) not equalto 0 for atleast one  i 
"FULL MODEL  is 
Y=beta0 + beta1*CRIM + beta2*ZN + beta3*INDUS + beta4*CHAS + beta5*NOX + beta6*rm + beta7*AGE + beta8*DIS + beta9*RAD + beta10*TAX 
+ beta11*PTRATIO + beta12*BLACK + beta13*LSTAT"
SSreg=sum((X%*%beta_hat-mean(Y))^2) #Residual sum of square under null hypthesis
SSres= sum((Y-(X%*%beta_hat))^2)                    #Residual sum of square for full model 
Fcal=(SSreg*(n-p))/((p-1)*SSres)   #Test statistics 
Fcal
alpha=0.05
Ftab=qf(1-alpha,p-1,n-p)

"(Decision Rule reject Ho if Fcal > Ftab   otherwise accept it) 
(Since Fcal=108.0767  is greater than Ftab=1.740074 under 5% level of significance so we reject Ho)
Thus we know there exist atleast one regressor which is significant"

#---------------------------------------------'sub-part end'-----------------------------------------------------------------------------------

  
#----------------------------------'Performing Individual regressor test'----------------------------------------------------------------------

C=solve(t(X)%*%X)
ftab=qf(1-alpha,1,n-p)  #tabulatted value of t dist. with (n-p-1)dof at 5% level of significance 
fcal=rep(0,p-1)
for (i in 1:(p-1))
  {
  fcal[i]=((beta_hat[i+1])^2)/((SSres/(n-p))*C[i+1,i+1])
  if(fcal[i]>ftab) print(" regressor is significant")
  else print(paste(colnames(X)[i+1],"regressor is not significant"))  
}

"As we see regressor corresponding to INDUS and AGE is not signifiacnt. Now , we consider new model with 
all regressor except the regressor corrresponding to INDUS and AGE i.e. beta3 and beta7 and including 
intercept term"



#-----------------------------------Model M1 is True or not -----------------------------------------------------
"Our New model is 
Y=beta0 + beta1*CRIM + beta2*ZN + beta4*CHAS + beta5*NOX + beta6*rm + beta8*DIS + beta9*RAD + beta10*TAX 
  + beta11*PTRATIO + beta12*BLACK + beta13*LSTAT"
#Hypothesis of interest is Ho:-"M1 is true" at 5% level of significance  
#Testing Ho:M1 is true  BETA2=0  where BETA=(BETA1,BETA2)'
Z=X[,c(-which(colnames(X)=='indus'),-which(colnames(X)=='age'))]
W=X[,c(which(colnames(X)=='indus'),which(colnames(X)=='age'))]
H=X%*%solve(t(X)%*%X)%*%t(X)
H1=Z%*%solve(t(Z)%*%Z)%*%t(Z)
SSregg=t(Y)%*%(H-H1)%*%Y    # Extra Sum of square due to BETA2
FcalN=(SSregg)*(n-p)/(SSres*2)
FtabN=qf(1-alpha,2,n-p)
if (FcalN<FtabN) print("M1 is TRUE MODEL")  else print(" M1 is FALSE MODEL ")



#----------------------CONFIDENCE INTERVAL OF NOX FROM FULL MODEL AND SLRM---------------------------------------------------------
cnox=which(colnames(X)=='nox')
cl=beta_hat[cnox]-sqrt(qf(1-alpha,1,n-p)*C[cnox,cnox]*(SSres/(n-p)))
cu=beta_hat[cnox]+sqrt(qf(1-alpha,1,n-p)*C[cnox,cnox]*(SSres/(n-p)))
CONFIDENCE_INTERVAL=c(cl,cu)
CONFIDENCE_INTERVAL

"Considering SLRM with predictor as NOX variable and response MEDV"
#model Y=b0 + b1*nox
cintercept=which(colnames(X)=='intercept')
X1=as.matrix(X[,c(cintercept,cnox)])
b_hat_nox=solve(t(X1)%*%X1)%*%t(X1)%*%Y
c1=solve(t(X1)%*%X1)
RSS1=t(Y-X1%*%b_hat_nox)%*%(Y-X1%*%b_hat_nox)
c1l=b_hat_nox[2]-(qt(0.975,n-2)*sqrt(c1[2,2]*(RSS1/(n-2))))
c1u=b_hat_nox[2]+(qt(0.975,n-2)*sqrt(c1[2,2]*(RSS1/(n-2))))
CON._INT=c(c1l,c1u)
CON._INT
#Since estimate of regressor lie in confidence interval so we say NOX is a significant estimator of predictor 

#--------------------------Testing the euality of regressor of CRIM and AGE -------------------------------------------------------------------------------
# "Considering Full Model"
# "Test Hypothesis Ho: beta(CRIRM)=beta(AGE) "
cld=(beta_hat[2]-beta_hat[8])-(qt(0.998,n-p)*(sqrt((SSres/(n-p))*(C[2,2]+C[8,8]-2*C[2,8]))))
cud=(beta_hat[2]-beta_hat[8])+(qt(0.998,n-p)*(sqrt((SSres/(n-p))*(C[2,2]+C[8,8]-2*C[2,8]))))
CON.INTd=c(cld,cud)
CON.INTd
# " This interval does not contain zero. So, there is  significant difference in the
# regressor coefficient  of CRIM AND AGE at 5% level of significance "

# "Considering model with standardized value of predictors"
Xd=matrix(0,n,p-1)
for(i in 1:p-1)
  {
    Xd[,i]=(X[,i+1]-mean(X[,i+1]))/sd(X[,i+1])
    }
Xd=cbind(rep(1,n),Xd)
beta_hatd=solve(t(Xd)%*%Xd)%*%t(Xd)%*%Y 
Cd=solve(t(Xd)%*%Xd)
RSSd= t(Y-Xd%*%beta_hatd)%*%(Y-Xd%*%beta_hatd)  
clsd=(beta_hatd[2]-beta_hatd[8])-(qt(0.975,n-p)*(sqrt((RSSd/(n-p))*(Cd[2,2]+Cd[8,8]-2*Cd[2,8]))))
cusd=(beta_hatd[2]-beta_hatd[8])+(qt(0.975,n-p)*(sqrt((RSSd/(n-p))*(Cd[2,2]+Cd[8,8]-2*Cd[2,8]))))
CON.INTsd=c(clsd,cusd)
CON.INTsd
# "Confidence interval does not contain value 0 .So, there is no change in decision also
# after standardizing the model there is  significant difference in the
# regressor coefficient  of CRIM AND AGE at 5% level of significance"


