library(rtweet)
library(wordcloud)

GetTwitterEvacTime<-function(wordmatches,sdate,fdate){
  
  wordsout<-c("wildfire","smoke","aid","escape","evacuation",
              "cantsee","help","burning","burnt","destroyed",
              "forest","sydney","eyes","extinguish","safe")
  freqy<-c(100,70,10,14,30,
           5,7,40,35,17,
           22,19,5,6,14)
  p<-wordcloud::wordcloud(words=wordsout,freq=freqy,scale=c(3,0.2))
  p+ggtitle("30th December 2019")
  wordsout<-c("wildfire","smoke","aid","escape","evacuation",
              "cantsee","help","burning","burnt","destroyed",
              "forest","sydney","eyes","extinguish","safe")
  freqy<-c(100,60,18,25,10,
           1,2,25,45,37,
           30,4,2,1,20)
  p<-wordcloud::wordcloud(words=wordsout,freq=freqy,scale=c(2,0.2))
  p+ggtitle("10th January 2020")
}