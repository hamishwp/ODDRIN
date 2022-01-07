library(mlbench)
library(caret)
library(tidyverse)
library(dplyr)
library(magrittr)

source("./RCode/GetSocioEconomic.R")
# source("./RCode/GetINFORM.R")

DispData<-readRDS("./IIDIPUS_Input/DispData_EQ_V2.Rdata")

# val<-GetWB_Vals(DispData)
val<-readRDS("./IIDIPUS_Input/val_WB.Rdata")
indies<-read_csv("./IIDIPUS_Input/REDUCED_WB-WorldDevelopment_Indicators.csv")
val%<>%filter(indicator%in%indies$indicator_id)

############################################################################################################
modifiers<-readRDS("./IIDIPUS_Results/ODDobjects/modifiers_extra.Rdata")
modifiers$eventid<-vapply(1:nrow(modifiers),function(i) as.numeric(strsplit(as.character(modifiers$event[i]),split = "_")[[1]][2]),numeric(1))
modifiers$modifier%<>%as.numeric()
modifiers$LL<-NULL
ids<-unique(val$eventid)[!unique(val$eventid)%in%unique(modifiers$eventid)]
ufiles<-list.files(path=paste0(dir,"IIDIPUS_Results/ODDobjects/"),pattern=Model$haz,recursive = T,ignore.case = T)
ufiles2<-vapply(1:length(ufiles),function(i) as.numeric(strsplit(ufiles[i],split = "_")[[1]][2]),numeric(1))
ufiles<-ufiles[ufiles2%in%ids]

for(i in 1:length(ufiles)){
  ODDy<-readRDS(paste0(dir,"IIDIPUS_Results/ODDobjects/",ufiles[i]))
  modifiers%<>%rbind(data.frame(iso3=ODDy@predictDisp$iso3,modifier=unlist(ODDy@modifier),
                                event=ufiles[i],gmax=ODDy@gmax$gmax,eventid=as.numeric(strsplit(ufiles[i],split = "_")[[1]][2]),
                                predictor=ODDy@predictDisp$predictor))
}

val%<>%filter(eventid%in%unique(modifiers$eventid))

modifiers<-val%>%group_by(eventid)%>%summarise(date=unique(date))%>%merge(modifiers,by="eventid")
modifiers$indicator<-"modifier"
modifiers%<>%transmute(eventid=eventid,iso3=iso3,nearval=modifier,indicator=indicator)
modifiers$iso3%<>%as.character()  
modifiers%<>%filter(!(eventid==310  & iso3=="UGA") &
                    !(eventid==660  & iso3=="PAK") &
                    !(eventid==1545 & iso3=="IRQ"))

val%<>%dplyr::select(eventid,iso3,nearval,indicator)%>%rbind(modifiers)
############################################################################################################
# x1<-filter(val,indicator=="NY.GDP.PCAP.CD") #; sum(is.na(x1$intrpval)) ; sum(is.na(x1$nearval)) ; print("...")
# x2<-filter(val,indicator=="SP.URB.TOTL.IN.ZS") #; sum(is.na(x2$intrpval)) ; sum(is.na(x2$nearval)) ; print("...")
# x3<-filter(val,indicator=="SL.UEM.TOTL.NE.ZS" ) #; sum(is.na(x3$intrpval)) ; sum(is.na(x3$nearval)) ; print("...")
# x4<-filter(val,indicator=="SM.POP.NETM" ) #; sum(is.na(x4$intrpval)) ; sum(is.na(x4$nearval)) ; print("...")
# modifiers<- x1
# modifiers$intrpval<-modifiers$intrpval*x2$intrpval/x3$intrpval*x4$intrpval
# modifiers$nearval<-modifiers$nearval*x2$nearval/x3$nearval*x4$nearval  
# modifiers$indicator<-modifiers$grouping<-"modifier"

# val%<>%rbind(modifiers)
############################################################################################################

nearval<-modifiers%>%dplyr::select(eventid,iso3,nearval); colnames(nearval)[3]<-"modifier"
# intrpval<-modifiers%>%dplyr::select(eventid,iso3,intrpval); colnames(intrpval)[3]<-"modifier"

for(inds in unique(val$indicator[val$indicator!="modifier"])){
  
  tmp<-filter(val,indicator==inds)
  tn<-data.frame(value=tmp$nearval); colnames(tn)<-inds
  nearval%<>%cbind(tn)
  # ti<-data.frame(value=tmp$intrpval); colnames(ti)<-inds
  # intrpval%<>%cbind(ti)
  
}
rm(tn,ti,tmp)

nearval[3:ncol(nearval)]<-scale(nearval[3:ncol(nearval)])
# intrpval[3:ncol(intrpval)]<-scale(intrpval[3:ncol(intrpval)])

tnearval<-nearval[,which(unname(apply(nearval,2,function(x) sum(!is.na(x))))>=160)]
tnearval<-tnearval[!apply(tnearval,1,function(x) any(is.na(x))),]
# tintrpval<-intrpval[,which(unname(apply(intrpval,2,function(x) sum(!is.na(x))))>163)]
# tintrpval<-tintrpval[!apply(tintrpval,1,function(x) any(is.na(x))),]
ncor <- cor(tnearval[,4:ncol(tnearval)])
# icor <- cor(tintrpval[,4:ncol(tintrpval)])
tnearval<-tnearval[,c(1,3,which(!1:ncol(tnearval)%in%caret::findCorrelation(ncor, cutoff=0.75)))]; tnearval<-tnearval[,c(1,3,2,4:ncol(tnearval))]
# tintrpval<-tintrpval[,!1:ncol(tintrpval)%in%findCorrelation(icor, cutoff=0.75)]
indepfeatures<-cor(tnearval[,4:ncol(tnearval)],tnearval$modifier)
ix<-sort(abs(indepfeatures),decreasing = T,index.return=T)$ix
tmp<-data.frame(cor=indepfeatures[ix],indicator_id=row.names(indepfeatures)[ix]); indepfeatures<-tmp; rm(tmp)
WBindies<-read_csv("./IIDIPUS_Input/WB-WorldDevelopment_Indicators.csv")[,c(1:5)]
indepfeatures%<>%merge(WBindies,by="indicator_id")
View(indepfeatures)
# rm(indepfeatures)
spearman <- apply(tnearval[,4:ncol(tnearval)],2,
                  function(x) (cor.test(x,
                                        tnearval$modifier,
                                        method = 'spearman'))$estimate)
spearman <- merge(data.frame(cor=spearman,
                             indicator_id=names(spearman)),
                  WBindies,by="indicator_id")
View(spearman)







# indepfeatures<-cor(tnearval[,5:ncol(tnearval)],tnearval$modifier)
# View(indepfeatures)

# final<-tnearval%>%dplyr::select(c(2,which(colnames(tnearval)%in%as.character(bestinds$indicator_id[1:30]))))
final<-tnearval%>%dplyr::select(-c(eventid,iso3))
colnames(final)[1]<-"Y"
weights<-1./Model$HighLevelPriors(Omega,Model,modifier = final$Y)+1;weights[weights>0]<-0; weights%<>%exp()

set.seed(998)
control <- trainControl(method="repeatedcv", number=5, repeats=15)

nbrnn <- train(modifier~., data=tnearval[,3:ncol(tnearval)], method='brnn',trControl=control,weights=weights)
nimportance<-varImp(nbrnn, scale=FALSE)
nimportance<-data.frame(cor=nimportance$importance,
                       indicator_id=row.names(nimportance$importance))%>%
                         merge(dplyr::select(WBindies,indicator,indicator_id),by="indicator_id")
View(nimportance)

nicr <- train(modifier~., data=tnearval[,3:ncol(tnearval)], method="icr",trControl=control,weights=weights)
nimportance<-varImp(nicr, scale=FALSE)
View(nimportance$importance)

nlm <- train(modifier~., data=tnearval[,3:ncol(tnearval)], method="lmStepAIC",trControl=control,weights=weights)
nimportance<-varImp(nlm, scale=FALSE)
View(nimportance$importance)

bestinds<-data.frame(indicator_id=row.names(nimportance$importance),rank=nimportance$importance$Overall,stringsAsFactors = F)%>%arrange(desc(rank))
bestinds%<>%merge(fresh_indicators,by="indicator_id")%>%arrange(desc(rank))
View(dplyr::select(bestinds,indicator_id,rank))

LMFeatureSelection<-function(output,Nb=12,intercept=F,fn="+",nlim=3,weights=NULL){
  # Nb - Number of times LM is run with different samples of training vs test data
  # intercept - do we include an intercept in the equation?
  
  # Use 80% of the observations as training and 20% for testing
  Ns <- floor(0.80*nrow(output))
  # List of variables in the output data.frame
  vars<-colnames(dplyr::select(output,-Y))
  # The parts of the LM that never change
  if(!intercept) eqn_base<-c("Y ~ ","+0") else eqn_base<-c("Y ~ ","")
  # Store the LM outputs
  prediction<-data.frame(eqn=NULL,AIC=NULL,BIC=NULL,adjR2=NULL)
  # How many variables to use at a time?
  for(n in 1:nlim){
    print(n)
    eqn<-combn(vars,n)
    if(n==1) {eqn%<>%as.character()
    }else eqn<-apply(eqn,2,function(x) pracma::strcat(x,collapse = fn))
    
    for(eee in eqn){
      equation<-paste0(eqn_base[1],eee,eqn_base[2])
      predictors<-data.frame(AIC=NULL,BIC=NULL)
      for (i in 1:Nb){
        ind = sample(seq_len(nrow(output)),size = Ns)
        train<-output[ind,]
        test<-output[-ind,]
        if(!is.null(weights)) predy<-lm(formula = as.formula(equation),
                                        data = train, weights = weights[ind])
        else predy<-lm(formula = as.formula(equation),data = train)
        # print(paste0(ceiling(AIC(predy)),", ",ceiling(BIC(predy))))
        # predme<-predict(predy,test, interval = 'confidence')
        predictors%<>%rbind(data.frame(AIC=AIC(predy),BIC=BIC(predy),adjR2=summary(predy)$adj.r.squared))
      }
      prediction<-rbind(cbind(data.frame(eqn=rep(equation)),t(colMeans(predictors))),prediction)
      
    }
  }
  
  return(prediction)
  
}

padd<-LMFeatureSelection(final,nlim = 3,weights = weights)
padd$eqn%<>%as.character()
padd$indicator_id<-vapply(1:nrow(padd),
                           function(i) trimws(strsplit(strsplit(padd$eqn[i],split = "+",fixed = T)[[1]][1],split="~",fixed = T)[[1]][2]),
                           character(1))
padd%<>%merge(dplyr::select(WBindies,indicator,indicator_id),by="indicator_id")
View(padd)

pprod<-LMFeatureSelection(final,nlim = 3,weights = weights,fn = ":",intercept = T)
pprod$eqn%<>%as.character()
# pprod$indicator_id<-vapply(1:nrow(pprod),
#                            function(i) trimws(strsplit(strsplit(pprod$eqn[i],split = "*",fixed = T)[[1]][1],split="~",fixed = T)[[1]][2]),
#                            character(1))
# pprod%<>%merge(dplyr::select(WBindies,indicator,indicator_id),by="indicator_id")
View(pprod)

pprodVIF<-LMFeatureSelection(final,nlim = 4,weights = weights,fn = ":",intercept = T)


# ilm <- train(modifier~., data=tintrpval[,2:ncol(tintrpval)], method="lmStepAIC",trControl=control)


# iimportance<-varImp(ilm, scale=FALSE)
View(nimportance$importance)
# View(iimportance$importance)
# plot(nimportance$importance$Overall[which(row.names(nimportance$importance)%in%row.names(iimportance$importance))],
#      iimportance$importance$Overall[which(row.names(iimportance$importance)%in%row.names(nimportance$importance))])
# Calculate the AUC_i*AUC_n and use these variables in the presentation?

summary(lm(obs~pred,data.frame(pred=predict(nbrnn,tnearval[,4:ncol(tnearval)]), obs=tnearval$modifier)))$adj.r.squared

summary(lm(data.frame(pred=predict(nlm, tnearval[,3:ncol(tnearval)]), obs=tnearval$modifier)))$adj.r.squared
# summary(lm(data.frame(pred=predict(ilm, tintrpval[,3:ncol(tintrpval)]), obs=tintrpval$modifier)))$adj.r.squared

nmodifier<-predict(nlm, tnearval[,2:ncol(tnearval)])
imodifier<-predict(ilm, tintrpval[,2:ncol(tintrpval)])

ggplot(data.frame(pred=predict(nbrnnREDVIF,tnearval[,5:ncol(tnearval)]),obs=tnearval$modifier),aes(pred,obs))+
  geom_point()+ylim(-2,2)+geom_smooth(method=lm , color="red", fill="blue",se=F)+xlab("Predicted")+ylab("Observed")


output<-readRDS("./IIDIPUS_Results/output_V2_nomodifier.Rdata")
output%<>%merge(tnearval,by=c("eventid","iso3"))
weights<-1./Model$HighLevelPriors(Omega,Model,modifier = output$modifier)+1;weights[weights>0]<-0; weights%<>%exp()

posthocDF<-dplyr::select(output,-c("eventid","iso3","qualifier","namer","ISO3C"))
posthocDF%<>%mutate(gmax=log(gmax+1),predictor=log(predictor+1))
Posthoc <- train(gmax~., data=posthocDF,
                 method='brnn',trControl=control,weights=weights)

nimportance<-varImp(Posthoc, scale=FALSE)
nimportance<-data.frame(cor=nimportance$importance,
                        indicator_id=row.names(nimportance$importance))%>%
  merge(dplyr::select(WBindies,indicator,indicator_id),by="indicator_id")
View(nimportance)

ggplot(data.frame(pred=exp(predict(Posthoc,posthocDF)), obs=exp(posthocDF$gmax)),aes(pred,obs))+
  geom_point()+geom_abline(slope=1)+xlab("Predicted")+ylab("Observed")+
  scale_y_log10()+scale_x_log10()

ggplot(output,aes(predictor,gmax))+
  geom_point()+geom_abline(slope=1)+xlab("Predicted")+ylab("Observed")+
  scale_y_log10()+scale_x_log10()

PosthocNoW <- train(gmax~., data=posthocDF,
                 method='brnn',trControl=control)

ggplot(data.frame(pred=exp(predict(PosthocNoW,posthocDF)), obs=exp(posthocDF$gmax)),aes(pred,obs))+
  geom_point()+geom_abline(slope=1)+xlab("Predicted")+ylab("Observed")+
  scale_y_log10()+scale_x_log10()

compositeI<-as.data.frame(cbind(predict(nbrnn,tnearval[,4:ncol(tnearval)]),tnearval$iso3)); names(compositeI)<-c("Index","iso3"); compositeI$Index%<>%as.character()%>%as.numeric()
compositeI%<>%group_by(iso3)%>%summarise(index=mean(Index,na.rm = T),Error=sd(Index,na.rm = T))
View(compositeI)
write_csv(compositeI,"./IIDIPUS_Results/BRNN_CompositeIndicator.csv")
