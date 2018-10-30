k <- subset(aug, leafAge=='mature' & treatment=='drought')
model <- gnls(gc ~ g0+(1+xi/sqrt(Dmmol))*(Photo/(CO2S-42.75)), start=list(g0=c(rep.int(0,14)), xi=c(rep.int(3,14))),
              param=list(xi~spp*treatment, g0~spp*treatment), data=k)
summary(model)

