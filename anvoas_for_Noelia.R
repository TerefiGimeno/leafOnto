library(doBy)
#basic fucntion
s.err <- function(x){
  err <- sd(x)/sqrt(length(x))
  return(err)
}
#read the data
lma <- read.csv()
thick <- read.csv()
#check that everthing is in order
head(lma)
str(lma)
head(thick)
str(thick)
#have a quick look at the data
hist(lma$LMA_gr_m2)
#it doesn't look very good try transformation
hist(log(lma$LMA_gr_m2))
hist(sqrt(lma$LMA_gr_m2))
hist((lma$LMA_gr_m2)^2)
#the best looking one is the log-transformation
#the same for thickness
hist(thick$thick_mm)
hist(log(thick$thick_mm))
hist(sqrt(thick$thick_mm))
hist((thick$thick_mm)^2)

#log-transform the response variables
lma$logLMA <- log(lma$LMA_gr_m2)
thick$logThick <- log(thick$thick_mm)
#ANOVA of the main effects
summary(aov(logLMA~spp*leafAge, data=lma))
summary(aov(logThick~spp*leafAge, data=thick))
#have a look at the means
summaryBy(logLMA~spp, FUN=mean, data=lma)
summaryBy(logThick~spp, FUN=mean, data=thick)
#now rank the data accordingly and run post hoc test
TukeyHSD(aov(logLMA~spp, data=lma))
TukeyHSD(aov(logThick~spp, data=thick))
#assing the letter
#get the data ready for the figure. The function length gives us the sample size per combination
lmaSumm <- summaryBy(LMA_gr_m2~spp*leafAge, data=lma, FUN=c(mean, s.err, length))
thickSumm <- summaryBy(thick~spp*leafAge, data=thick, FUN=c(mean, s.err, length))
#export the data
write.csv(lmaSumm, file='', row.names=FALSE)
write.csv(thickSumm, file='', row.names=FALSE)