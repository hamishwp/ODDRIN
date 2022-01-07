# Number of time the LM is run with different samples of training vs test data
Nb<-20
# Percentage of observations used for training vs test dataset
Ns = floor(0.80*nrow(output))
# Variables in the output data.frame
vars<-c("CC.INF.PHY","CC.INS.DRR","CC.INS.GOV.GE","HA.NAT.EQ","Ik","VU.SEV.AD","VU.SEV.PD",
            "p10p100","p20p100","p30p100","p40p100","p50p100","p60p100","p70p100","p80p100","p90p100")

# Relative error
output$Y<-log(1e-10+output$gmax/output$predictor)
# Log Likelihood
output$iLL<-abs(output$predictor-output$gmax)/output$gmax
# The components of the LM that never change
eqn_base<-c("Y ~ ","+0")
   
prediction<-data.frame(eqn=NULL,LL_Ex=NULL,LL_Lo=NULL,LL_Up=NULL,AIC=NULL,BIC=NULL)
# For each number of variables to be used in the LM (combinations of 1 to 6 different variables)
for(n in 1:6){
  # How many different combinations can you make with n different variables?
  eqn<-combn(vars,n)
  # Write up the equation by combining into a purely additive model
  if(n==1) {eqn%<>%as.character()
  }else eqn%<>%apply(2,function(x) strcat(x,collapse = "+"))
  # For each equation applied
  for(eee in eqn){
    # Write out the equation
    equation<-paste0(eqn_base[1],eee,eqn_base[2])
    predictors<-data.frame(LL_Ex=NULL,LL_Lo=NULL,LL_Up=NULL,AIC=NULL,BIC=NULL)
    # Split test and training and run LM many times to ensure accurate parameterisation
    for (i in 1:Nb){
      # Split the test from training data by randomly sampling
      ind = sample(seq_len(nrow(output)),size = Ns)
      train<-output[ind,]
      test<-output[-ind,]
      # Run the LM MLE
      predy<-lm(formula = as.formula(equation),data = train)
      # Extract confidence intervals
      predme<-predict(predy,test, interval = 'confidence')
      # Convert output prediction from log-likelihood to observed error form
      me<-colMeans(apply(predme,2,function(x) abs(exp(x)-exp(test$Y))/(exp(test$Y)+exp(x))),na.rm = T)
      # me<-colMeans(apply(predme,2,function(x) log(1+(exp(x)-exp(test$Y))^2/(exp(test$Y)+exp(x)))),na.rm = T)
      
      predictors%<>%rbind(data.frame(LL_Ex=me[1],LL_Lo=me[2],LL_Up=me[3],AIC=AIC(predy),BIC=BIC(predy)))
    }
    
    prediction<-rbind(cbind(data.frame(eqn=rep(equation)),t(colMeans(predictors))),prediction)
    
  }
  print(n)
  print(min(prediction$BIC))
  print(as.character(prediction$eqn[which.min(prediction$BIC)]))
}

# names(prediction)<-c("LL_Ex","LL_Lo","LL_Up","AIC","BIC")

# predy<-lm(formula = as.formula("Y ~ CC.INS.DRR+CC.INS.GOV.GE+0"),data = output) # NO CHN OR IDN
predy<-lm(formula = as.formula("Y ~ CC.INF.PHY+CC.INS.DRR+VU.SEV.AD+0"),data = output) # USING log-OLS
# predy<-lm(formula = as.formula("Y ~ CC.INF.PHY+CC.INS.DRR+CC.INS.GOV.GE+VU.SEV.AD+0"),data = output) # USING ROOT-OLS
# predy<-lm(formula = as.formula("Y ~ CC.INF.PHY+CC.INS.DRR+p80p100+0"),data = output)
# predy<-lm(formula = as.formula("Y ~ CC.INS.DRR+CC.INS.GOV.GE+p60p100+0"),data = output)
print(paste0("BIC = ",3*log(length(output$gmax))-
               2*log(mean((exp(predy$fitted.values)*output$predictor-output$gmax)/
                             (exp(predy$fitted.values)*output$predictor+output$gmax)))))
# predy<-lm(formula = as.formula("Y ~ p70p100+0"),data = output)
# print(paste0("BIC = ",log(length(output$gmax))-
#                2*log(mean((exp(predy$fitted.values)*output$predictor-output$gmax)/
#                              (exp(predy$fitted.values)*output$predictor+output$gmax)))))

# min(prediction$BIC)->409.7246
# prediction$eqn[which.min(prediction$BIC)]  ->   "Y ~ CC.INF.PHY+CC.INS.DRR+p80p100+0"
# predy$coefficients
# CC.INF.PHY CC.INS.DRR  VU.SEV.AD 
# 25.55795  -17.84694  -31.56208
# CC.INF.PHY    CC.INS.DRR CC.INS.GOV.GE     VU.SEV.AD 
# 16.83771     -38.08718      21.87292     -33.28555
### CC.INF.PHY CC.INS.DRR    p80p100 
### 24.986209 -34.787344   1.667839 

output$predmod<-output$predictor*exp(25.55795*output$CC.INF.PHY + -17.84694*output$CC.INS.DRR +-31.56208*output$VU.SEV.AD) #LOG-OLS
# output$predmod<-output$predictor*exp(-37.32581*output$CC.INS.DRR + 26.72798*output$CC.INS.GOV.GE) # NO CHN OR IDN 
output$iLL_mod<-abs(output$predmod-output$gmax)/output$gmax

ggplot(output,aes(predmod,gmax))+geom_point() + geom_abline(slope=1) +
  scale_y_log10() + scale_x_log10() + xlab("Prediction Displaced") + ylab("Observed Displaced")

ggplot(filter(output,iso3!="CHN"),aes(predmod,gmax))+geom_point() + geom_abline(slope=1) +
  scale_y_log10() + scale_x_log10() + xlab("Prediction Displaced") + ylab("Observed Displaced")

ggplot(output,aes(predictor,gmax))+geom_point() + geom_abline(slope=1) +
  scale_y_log10() + scale_x_log10() + xlab("Prediction Displaced") + ylab("Observed Displaced")
