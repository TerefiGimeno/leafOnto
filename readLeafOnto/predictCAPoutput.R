dataToPredict <- expand.grid(seq(from=min(gasex$Photo, na.rm=T), to=max(gasex$Photo, na.rm=T), length.out=100),
                               seq(from=min(gasex$VpdL, na.rm=T), to=max(gasex$VpdL, na.rm=T), length.out=100))
names(dataToPredict) <-c('A','D')
dataToPredict$CAPindex <- dataToPredict$A/(sqrt(dataToPredict$D)*(median(gasex$CO2S)-gammaStar))
predictCAPoutput <- function(g0, Ep, Ca, df){
  gcPred <- g0+(1+Ep/sqrt(df$D))*(df$A/(Ca-42.75))
  return(gcPred)
}
