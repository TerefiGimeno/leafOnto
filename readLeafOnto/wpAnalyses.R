source('readLeafOnto/readLibsLeafOnto.R')
source('readLeafOnto/basicFun.R')
wp <- read.csv('dataLeafOnto/WPglashouseDrought.csv')
wp$wpNew <- wp$WP_predawn
wp$wpNew[which(is.na(wp$WP_predawn))] <- -2#this is the WP_cri on table 1 in Dewar et al. 2017
wp$wpNew[which(wp$Notes=='soil is wet!!!')] <- NA
wpSummary <- as.data.frame(summarise(group_by(wp, Spp, Treatment), wpMean=mean.na(wpNew),wpSE=s.err(wpNew)))
wpSummary2 <- wpSummary
wpSummary$all <- paste0(round(wpSummary$wpMean, digits=2), " ± ", round(wpSummary$wpSE, digits=2))
wpSummary <- wpSummary[,c('Spp','Treatment','all')]
wpSummaryWide <- reshape(wpSummary, v.names='all', idvar='Spp', timevar='Treatment', direction='wide')
names(wpSummaryWide) <- c('Spp', 'wpControl','wpDrought')

wpTest <- read.csv('dataLeafOnto/wpTest.csv')
wpTestNames <- unique(wpTest[,c('Spp','id')])
wpTestSumm <- as.data.frame(summarise(group_by(wpTest, id), min=min(wp_MPa, na.rm=T), max=max(wp_MPa, na.rm=T)))
wpTestSumm <- merge(wpTestSumm, wpTestNames, by='id', all=T)
wpTestS <- as.data.frame(summarise(group_by(wpTestSumm, Spp),minM=mean(min), minSE=s.err(min),
                                   maxM=mean(max), maxSE=s.err(max)))
wpSummary2[(nrow(wpSummary2)+1),] <- c('EUC','Control', wpTestS[which(wpTestS$Spp=='EUC'),'maxM'],
                                       wpTestS[which(wpTestS$Spp=='EUC'),'maxSE'])
wpSummary2$Treatment <- ifelse(wpSummary2$Treatment=='Control','well-watered','drought')
names(wpSummary2)[1] <- 'spp'
wpTestS$maxAll <- paste0(round(wpTestS$maxM, digits=2), " ± ", round(wpTestS$maxSE, digits=2))
wpTestS$minAll <- paste0(round(wpTestS$minM, digits=2), " ± ", round(wpTestS$minSE, digits=2))
wpTestS <- wpTestS[,c('Spp','maxAll','minAll')]

Table1 <- data.frame(row.names=1:7)
Table1$Spp <- c('PTE','BET','QUE','EUC','PIN','MOL','ZEA')
Table1$Species <- c('P. aquilinum', "B. pendula", "Q. robur","E. gunnii", "P. pinaster", "M. caerulea", "Z. mays")
Table1 <- merge(merge(Table1, wpTestS, by='Spp', all=T),wpSummaryWide, by='Spp', all=T)
write.csv(Table1[,2:ncol(Table1)], row.names=FALSE, file='outputLeafOnto/Table1.csv')
summary(aov(wpNew~Spp*Treatment, data=wp))

