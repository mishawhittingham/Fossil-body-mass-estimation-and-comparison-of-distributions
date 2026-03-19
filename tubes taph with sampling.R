library(moments)
library(multimode)
library(Matching)

BM_mamms = read.csv("C:\\Wirc\\Projekt3\\Bigness\\subsets\\TetBMs_tubesTaphready.csv", header = TRUE)
BM_mamms20 = read.csv("C:\\Wirc\\Projekt3\\Bigness\\subsets\\TetBMs_tubesTaphready_size20.csv", header = TRUE)

sites = 20
sites20 = 10


Foss = read.csv("C:\\Wirc\\Projekt3\\Bigness\\FossilBMs.csv", header = TRUE)

allsmallmodes = c()
allsmallkurts = c()
allsmallskews = c()
allsmallmeds = c()
allbothmodes = c()
allbothkurts = c()
allbothskews = c()
allbothmeds = c()
allbasemodes = c()
allbasekurts = c()
allbaseskews = c()
allbasemeds = c()

allsmallmodes20 = c()
allsmallkurts20 = c()
allsmallskews20 = c()
allsmallmeds20 = c() 
allbothmodes20 = c()
allbothkurts20 = c()
allbothskews20 = c()
allbothmeds20  = c()
allbasemodes20 = c()
allbasekurts20 = c()
allbaseskews20 = c()
allbasemeds20 = c()

allsmallJOGmaxps = c()
allsmallJOGmaxds = c()

allbothPAmaxps = c()
allbothPAmaxds = c()
allbasePAmaxps = c()
allbasePAmaxds = c()

allsmallJOGmaxps20 = c()
allsmallJOGmaxds20 = c()

allbothPAmaxps20 = c()
allbothPAmaxds20 = c()
allbasePAmaxps20 = c()
allbasePAmaxds20 = c()


#now for CSAmax estimates

for(i in 1:sites){
  baseset = na.omit(BM_mamms[,i])
  smallmaxfit5 = ifelse(BM_mamms[,i] < min(Foss[,11]), 0.05, 0.95)
  bothmaxfit5 = ifelse(BM_mamms[,i] <  max(na.omit(Foss[,6])) & BM_mamms[,i] > min(Foss[,11]),0.95, 0.05)
  smallps = c()
  smallds = c()
  smallmodes = c()
  smallkurts = c()
  smallskews = c()
  smallmeds = c()
  bothmodes = c()
  bothkurts = c()
  bothskews = c()
  bothmeds = c()
  basemodes = c()
  basekurts = c()
  baseskews = c()
  basemeds = c()
  for(j in 1:100){
    smallsampleset = sample(baseset, 12, replace = FALSE, prob = na.omit(smallmaxfit5))
    smallresult = ks.boot(log10(na.omit(Foss[,4])), log10(smallsampleset), nboots = 1000)
    smallmode = as.numeric(modetest(log10(smallsampleset), mod0 = 1)$statistic)
    smallkurt = as.numeric(anscombe.test(log10(smallsampleset))$statistic[1])
    smallskew = as.numeric(agostino.test(log10(smallsampleset), alternative = "two.sided")$statistic[1])
    smallmodes = c(smallmodes, smallmode)
    smallkurts = c(smallkurts, smallkurt)
    smallskews = c(smallskews, smallskew)
    smallmeds = c(smallmeds, median(smallsampleset))
    smallpval = smallresult$ks$p.value
    smalldval = smallresult$ks$statistic
    smallps = c(smallps, smallpval)
    smallds = c(smallds, smalldval)
    bothsampleset = sample(baseset, 12, replace = FALSE, prob = na.omit(bothmaxfit5))
    bothmode = as.numeric(modetest(log10(bothsampleset), mod0 = 1)$statistic)
    bothkurt = as.numeric(anscombe.test(log10(bothsampleset))$statistic[1])
    bothskew = as.numeric(agostino.test(log10(bothsampleset), alternative = "two.sided")$statistic[1])
    bothmodes = c(bothmodes, bothmode)
    bothkurts = c(bothkurts, bothkurt)
    bothskews = c(bothskews, bothskew)
    bothmeds = c(bothmeds, median(bothsampleset))
    basesampleset = sample(baseset, 12, replace = FALSE)
    basemode = as.numeric(modetest(log10(basesampleset), mod0 = 1)$statistic)
    basekurt = as.numeric(anscombe.test(log10(basesampleset))$statistic[1])
    baseskew = as.numeric(agostino.test(log10(basesampleset), alternative = "two.sided")$statistic[1])
    basemodes = c(basemodes, basemode)
    basekurts = c(basekurts, basekurt)
    baseskews = c(baseskews, baseskew)
    basemeds = c(basemeds, median(basesampleset))
  }
  allsmallJOGmaxps = c(allsmallJOGmaxps, smallps)
  allsmallJOGmaxds = c(allsmallJOGmaxds, smallds)
  allsmallmodes = c(allsmallmodes, smallmodes)
  allsmallkurts = c(allsmallkurts, smallkurts)
  allsmallskews = c(allsmallskews, smallskews)
  allbothmodes = c(allbothmodes, bothmodes)
  allbothkurts = c(allbothkurts, bothkurts)
  allbothskews = c(allbothskews, bothskews)
  allbasemodes = c(allbasemodes, basemodes)
  allbasekurts = c(allbasekurts, basekurts)
  allbaseskews = c(allbaseskews, baseskews)
  allsmallmeds = c(allsmallmeds, smallmeds)
  allbothmeds = c(allbothmeds, bothmeds)
  allbasemeds = c(allbasemeds, basemeds)
}

#now do PA
for(i in 1:sites){
  baseset = na.omit(BM_mamms[,i])
  bothmaxfit5 = ifelse(BM_mamms[,i] <  max(na.omit(Foss[,6])) & BM_mamms[,i] > min(Foss[,11]),0.95, 0.05)
  bothps = c()
  bothds = c()
  baseps = c()
  baseds = c()
  for(j in 1:100){
    bothsampleset = sample(baseset, 9, replace = FALSE, prob = na.omit(bothmaxfit5))
    bothresult = ks.boot(log10(na.omit(Foss[,10])), log10(bothsampleset), nboots = 1000)
    bothpval = bothresult$ks$p.value
    bothdval = bothresult$ks$statistic
    bothps = c(bothps, bothpval)
    bothds = c(bothds, bothdval)
    basesampleset = sample(baseset, 9, replace = FALSE)
    baseresult = ks.boot(log10(na.omit(Foss[,10])), log10(basesampleset), nboots = 1000)
    basepval = baseresult$ks$p.value
    basedval = baseresult$ks$statistic
    baseps = c(baseps, basepval)
    baseds = c(baseds, basedval)
  }
  allbothPAmaxps = c(allbothPAmaxps, bothps)
  allbothPAmaxds = c(allbothPAmaxds, bothds)
  allbasePAmaxps = c(allbasePAmaxps, baseps)
  allbasePAmaxds = c(allbasePAmaxds, baseds)
}

#now at sample size 20

for(i in 1:sites20){
  baseset = na.omit(BM_mamms20[,i])
  smallmaxfit5 = ifelse(BM_mamms20[,i] < min(Foss[,11]), 0.05, 0.95)
  bothmaxfit5 = ifelse(BM_mamms20[,i] <  max(na.omit(Foss[,6])) & BM_mamms20[,i] > min(Foss[,11]),0.95, 0.05)
  smallps = c()
  smallds = c()
  smallmodes = c()
  smallkurts = c()
  smallskews = c()
  smallmeds = c()
  bothmodes = c()
  bothkurts = c()
  bothskews = c()
  bothmeds = c()
  basemodes = c()
  basekurts = c()
  baseskews = c()
  basemeds = c()
  for(j in 1:100){
    smallsampleset = sample(baseset, 20, replace = FALSE, prob = na.omit(smallmaxfit5))
    smallresult = ks.boot(log10(na.omit(Foss[,4])), log10(smallsampleset), nboots = 1000)
    smallmode = as.numeric(modetest(log10(smallsampleset), mod0 = 1)$statistic)
    smallkurt = as.numeric(anscombe.test(log10(smallsampleset))$statistic[1])
    smallskew = as.numeric(agostino.test(log10(smallsampleset), alternative = "two.sided")$statistic[1])
    smallmodes = c(smallmodes, smallmode)
    smallkurts = c(smallkurts, smallkurt)
    smallskews = c(smallskews, smallskew)
    smallmeds = c(smallmeds, median(smallsampleset))
    smallpval = smallresult$ks$p.value
    smalldval = smallresult$ks$statistic
    smallps = c(smallps, smallpval)
    smallds = c(smallds, smalldval)
    bothsampleset = sample(baseset, 20, replace = FALSE, prob = na.omit(bothmaxfit5))
    bothmode = as.numeric(modetest(log10(bothsampleset), mod0 = 1)$statistic)
    bothkurt = as.numeric(anscombe.test(log10(bothsampleset))$statistic[1])
    bothskew = as.numeric(agostino.test(log10(bothsampleset), alternative = "two.sided")$statistic[1])
    bothmodes = c(bothmodes, bothmode)
    bothkurts = c(bothkurts, bothkurt)
    bothskews = c(bothskews, bothskew)
    bothmeds = c(bothmeds, median(bothsampleset))
    basesampleset = sample(baseset, 20, replace = FALSE)
    basemode = as.numeric(modetest(log10(basesampleset), mod0 = 1)$statistic)
    basekurt = as.numeric(anscombe.test(log10(basesampleset))$statistic[1])
    baseskew = as.numeric(agostino.test(log10(basesampleset), alternative = "two.sided")$statistic[1])
    basemodes = c(basemodes, basemode)
    basekurts = c(basekurts, basekurt)
    baseskews = c(baseskews, baseskew)
    basemeds = c(basemeds, median(basesampleset))
  }
  allsmallJOGmaxps20 = c(allsmallJOGmaxps20, smallps)
  allsmallJOGmaxds20 = c(allsmallJOGmaxds20, smallds)
  allsmallmodes20 = c(allsmallmodes20, smallmodes)
  allsmallkurts20 = c(allsmallkurts20, smallkurts)
  allsmallskews20 = c(allsmallskews20, smallskews)
  allbothmodes20 = c(allbothmodes20, bothmodes)
  allbothkurts20 = c(allbothkurts20, bothkurts)
  allbothskews20 = c(allbothskews20, bothskews)
  allbasemodes20 = c(allbasemodes20, basemodes)
  allbasekurts20 = c(allbasekurts20, basekurts)
  allbaseskews20 = c(allbaseskews20, baseskews)
  allsmallmeds20 = c(allsmallmeds20, smallmeds)
  allbothmeds20 = c(allbothmeds20, bothmeds)
  allbasemeds20 = c(allbasemeds20, basemeds)
}

#now do PA
for(i in 1:sites20){
  baseset = na.omit(BM_mamms20[,i])
  bothmaxfit5 = ifelse(BM_mamms20[,i] <  max(na.omit(Foss[,6])) & BM_mamms20[,i] > min(Foss[,11]),0.95, 0.05)
  bothps = c()
  bothds = c()
  baseps = c()
  baseds = c()
  for(j in 1:100){
    bothsampleset = sample(baseset, 20, replace = FALSE, prob = na.omit(bothmaxfit5))
    bothresult = ks.boot(log10(na.omit(Foss[,10])), log10(bothsampleset), nboots = 1000)
    bothpval = bothresult$ks$p.value
    bothdval = bothresult$ks$statistic
    bothps = c(bothps, bothpval)
    bothds = c(bothds, bothdval)
    basesampleset = sample(baseset, 20, replace = FALSE)
    baseresult = ks.boot(log10(na.omit(Foss[,10])), log10(basesampleset), nboots = 1000)
    basepval = baseresult$ks$p.value
    basedval = baseresult$ks$statistic
    baseps = c(baseps, basepval)
    baseds = c(baseds, basedval)
  }
  
  allbothPAmaxps20 = c(allbothPAmaxps20, bothps)
  allbothPAmaxds20 = c(allbothPAmaxds20, bothds)
  allbasePAmaxps20 = c(allbasePAmaxps20, baseps)
  allbasePAmaxds20 = c(allbasePAmaxds20, baseds)
}

toprow = c("small-bias JOGmax", "small-bias JOGmax size20", "edge-bias PAmax", "edge-bias PAmax size 20", 
             "unbiased PAmax", "unbiased PAmax size 20")

poutput = data.frame(matrix(nrow = (sites*100), ncol = 6))
poutput[,1] = allsmallJOGmaxps
poutput[1:(sites20*100),2] = allsmallJOGmaxps20
poutput[,3] = allbothPAmaxps
poutput[1:(sites20*100),4] = allbothPAmaxps20
poutput[,5] = allbasePAmaxps
poutput[1:(sites20*100),6] = allbasePAmaxps20

colnames(poutput) = toprow


doutput = data.frame(matrix(nrow = (sites*100), ncol = 6))
doutput[,1] = allsmallJOGmaxds
doutput[1:(sites20*100),2] = allsmallJOGmaxds20
doutput[,3] = allbothPAmaxds
doutput[1:(sites20*100),4] = allbothPAmaxds20
doutput[,5] = allbasePAmaxds
doutput[1:(sites20*100),6] = allbasePAmaxds20

colnames(doutput) = toprow

toprow2 = c("small-bias", "small-bias size 20", "edge-bias", "edge-bias size 20",  
           "unbiased", "unbiased size 20")

moutput = data.frame(matrix(nrow = (sites*100), ncol = 6))
moutput[,1] = allsmallmodes
moutput[1:(sites20*100),2] = allsmallmodes20
moutput[3] = allbothmodes
moutput[1:(sites20*100),4] = allbothmodes20
moutput[,5] = allbasemodes
moutput[1:(sites20*100),6] = allbasemodes20

colnames(moutput) = toprow2

koutput = data.frame(matrix(nrow = (sites*100), ncol = 6))
koutput[,1] = allsmallkurts
koutput[1:(sites20*100),2] = allsmallkurts20
koutput[3] = allbothkurts
koutput[1:(sites20*100),4] = allbothkurts20
koutput[,5] = allbasekurts
koutput[1:(sites20*100),6] = allbasekurts20

colnames(koutput) = toprow2


soutput = data.frame(matrix(nrow = (sites*100), ncol = 6))
soutput[,1] = allsmallskews
soutput[1:(sites20*100),2] = allsmallskews20
soutput[3] = allbothskews
soutput[1:(sites20*100),4] = allbothskews20
soutput[,5] = allbaseskews
soutput[1:(sites20*100),6] = allbaseskews20

colnames(soutput) = toprow2


aoutput = data.frame(matrix(nrow = (sites*100), ncol = 6))
aoutput[,1] = allsmallmeds
aoutput[1:(sites20*100),2] = allsmallmeds20
aoutput[3] = allbothmeds
aoutput[1:(sites20*100),4] = allbothmeds20
aoutput[,5] = allbasemeds
aoutput[1:(sites20*100),6] = allbasemeds20

colnames(aoutput) = toprow2


write.csv(poutput, "C:\\Wirc\\Projekt3\\Bigness\\subsets\\pvals_tubesTaph_pgls.csv")
write.csv(doutput, "C:\\Wirc\\Projekt3\\Bigness\\subsets\\dvals_tubesTaph_pgls.csv")
write.csv(moutput, "C:\\Wirc\\Projekt3\\Bigness\\subsets\\modes_tubesTaph_pgls.csv")
write.csv(koutput, "C:\\Wirc\\Projekt3\\Bigness\\subsets\\kurts_tubesTaph_pgls.csv")
write.csv(soutput, "C:\\Wirc\\Projekt3\\Bigness\\subsets\\skews_tubesTaph_pgls.csv")
write.csv(aoutput, "C:\\Wirc\\Projekt3\\Bigness\\subsets\\medians_tubesTaph_pgls.csv")
