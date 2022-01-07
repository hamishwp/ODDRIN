




stop("add an extra attribute to ODDy - indicator of which population data was used (worldpop, FB, etc)")




folder<-"./IIDIPUS_Results/ODDobjects_V2_FB/"
filez<-list.files(folder)

for(fff in filez){
  
  ODDy<-readRDS(paste0(dir,"IIDIPUS_Results/ODDobjects_V2_FB/",fff))  
  # Set all modifiers to zero
  ODDy@modifier
  # Save SEDAC population just in case!
  ODDy$SEDACpop<-ODDy$Population
  # If China, load WorldPop data instead
  if("CHN" %in% unique(na.omit(ODDy@data$ISO3C))){
    
    # Automatically download from the worldpop dataset
    ODDy<-AddWorldPop(ODDy)
    stop()
    # then set worldpop as main population
    ODDy$Population<-ODDy$WorldPop
    # Update ODDy@PopSource
    stop()
    
    
  } else {
    
    # Set FBpop as main population
    ODDy$Population<-ODDy$FBPop
    # Update ODDy@PopSource
    stop()
    
  }
  
  # check the differences between sedac and fbpop
  tmp<-ODDy@data%>%group_by(ISO3C)%>%summarise(SED=sum(SEDACpop,na.rm = T),
                                          FB=sum(FBPop,na.rm = T))
  tmp$diff=abs(tmp$SED-tmp$FBPop)
  tmp
  
  # plot out to a bg population contour plot
  p<-plotODDyBG(ODDy)
  ggsave(paste0(fff,".png"), plot=p,path = "./Plots/PopFB_SEDACS/", width = 9,height = 7.)
  
  saveRDS(ODDy,paste0("./IIDIPUS_Input/ODDobjects_V2_variedPop/",fff))
  
}





