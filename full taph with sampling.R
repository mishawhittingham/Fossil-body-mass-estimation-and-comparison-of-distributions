library(moments)
library(multimode)
library(Matching)

BM_mamms = read.csv("C:\\Wirc\\Projekt3\\Bigness\\subsets\\TetBMs_withmamms.csv", header = TRUE)
BM_mamms20 = read.csv("C:\\Wirc\\Projekt3\\Bigness\\subsets\\TetBMs_withmamms_size20.csv", header = TRUE)
BM_mamms40 = read.csv("C:\\Wirc\\Projekt3\\Bigness\\subsets\\TetBMs_withmamms_size40.csv", header = TRUE)
BM_mamms60 = read.csv("C:\\Wirc\\Projekt3\\Bigness\\subsets\\TetBMs_withmamms_size60.csv", header = TRUE)
Foss = read.csv("C:\\Wirc\\Projekt3\\Bigness\\FossilBMs.csv", header = TRUE)

allsmallmodes = c()
allsmallkurts = c()
allsmallskews = c()
allsmallmeds = c()
allbigmodes = c()
allbigkurts = c()
allbigskews = c()
allbigmeds = c()
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
allbigmodes20 = c()
allbigkurts20 = c()
allbigskews20 = c()
allbigmeds20 = c()
allbothmodes20 = c()
allbothkurts20 = c()
allbothskews20 = c()
allbothmeds20  = c()
allbasemodes20 = c()
allbasekurts20 = c()
allbaseskews20 = c()
allbasemeds20 = c()

allsmallmodes40 = c()
allsmallkurts40 = c()
allsmallskews40 = c()
allsmallmeds40 = c()
allbigmodes40 = c()
allbigkurts40 = c()
allbigskews40 = c()
allbigmeds40 = c()
allbothmodes40 = c()
allbothkurts40 = c()
allbothskews40 = c()
allbothmeds40 = c()
allbasemodes40 = c()
allbasekurts40 = c()
allbaseskews40 = c()
allbasemeds40 = c()

allsmallmodes60 = c()
allsmallkurts60 = c()
allsmallskews60 = c()
allsmallmeds60 = c()
allbigmodes60 = c()
allbigkurts60 = c()
allbigskews60 = c()
allbigmeds60 = c()
allbothmodes60 = c()
allbothkurts60 = c()
allbothskews60 = c()
allbothmeds60 = c()
allbasemodes60 = c()
allbasekurts60 = c()
allbaseskews60 = c()
allbasemeds60 = c()

allsmallJOGmaxps = c()
allsmallJOGmaxds = c()
allbigJOGmaxps = c()
allbigJOGmaxds = c()
allbothJOGmaxps = c()
allbothJOGmaxds = c()
allbaseJOGmaxps = c()
allbaseJOGmaxds = c()

allsmallPAmaxps = c()
allsmallPAmaxds = c()
allbigPAmaxps = c()
allbigPAmaxds = c()
allbothPAmaxps = c()
allbothPAmaxds = c()
allbasePAmaxps = c()
allbasePAmaxds = c()

allsmallJOGmaxps20 = c()
allsmallJOGmaxds20 = c()
allbigJOGmaxps20 = c()
allbigJOGmaxds20 = c()
allbothJOGmaxps20 = c()
allbothJOGmaxds20 = c()
allbaseJOGmaxps20 = c()
allbaseJOGmaxds20 = c()

allsmallPAmaxps20 = c()
allsmallPAmaxds20 = c()
allbigPAmaxps20 = c()
allbigPAmaxds20 = c()
allbothPAmaxps20 = c()
allbothPAmaxds20 = c()
allbasePAmaxps20 = c()
allbasePAmaxds20 = c()

allsmallJOGmaxps40 = c()
allsmallJOGmaxds40 = c()
allbigJOGmaxps40 = c()
allbigJOGmaxds40 = c()
allbothJOGmaxps40 = c()
allbothJOGmaxds40 = c()
allbaseJOGmaxps40 = c()
allbaseJOGmaxds40 = c()

allsmallPAmaxps40 = c()
allsmallPAmaxds40 = c()
allbigPAmaxps40 = c()
allbigPAmaxds40 = c()
allbothPAmaxps40 = c()
allbothPAmaxds40 = c()
allbasePAmaxps40 = c()
allbasePAmaxds40 = c()

allsmallJOGmaxps60 = c()
allsmallJOGmaxds60 = c()
allbigJOGmaxps60 = c()
allbigJOGmaxds60 = c()
allbothJOGmaxps60 = c()
allbothJOGmaxds60 = c()
allbaseJOGmaxps60 = c()
allbaseJOGmaxds60 = c()

allsmallPAmaxps60 = c()
allsmallPAmaxds60 = c()
allbigPAmaxps60 = c()
allbigPAmaxds60 = c()
allbothPAmaxps60 = c()
allbothPAmaxds60 = c()
allbasePAmaxps60 = c()
allbasePAmaxds60 = c()


#now for CSAmax estimates

for(i in 1:29){
  baseset = na.omit(BM_mamms[,i])
  smallmaxfit5 = ifelse(BM_mamms[,i] < min(Foss[,11]), 0.05, 0.95)
  bigmaxfit5 = ifelse(BM_mamms[,i] > max(na.omit(Foss[,6])), 0.05, 0.95)
  bothmaxfit5 = ifelse(BM_mamms[,i] <  max(na.omit(Foss[,6])) & BM_mamms[,i] > min(Foss[,11]),0.95, 0.05)
  smallps = c()
  smallds = c()
  smallmodes = c()
  smallkurts = c()
  smallskews = c()
  smallmeds = c()
  bigps = c()
  bigds = c()
  bigmodes = c()
  bigkurts = c()
  bigskews = c()
  bigmeds = c()
  bothps = c()
  bothds = c()
  bothmodes = c()
  bothkurts = c()
  bothskews = c()
  bothmeds = c()
  baseps = c()
  baseds = c()
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
    bigsampleset = sample(baseset, 12, replace = FALSE, prob = na.omit(bigmaxfit5))
    bigresult = ks.boot(log10(na.omit(Foss[,4])), log10(bigsampleset), nboots = 1000)
    bigmode = as.numeric(modetest(log10(bigsampleset), mod0 = 1)$statistic)
    bigkurt = as.numeric(anscombe.test(log10(bigsampleset))$statistic[1])
    bigskew = as.numeric(agostino.test(log10(bigsampleset), alternative = "two.sided")$statistic[1])
    bigmodes = c(bigmodes, bigmode)
    bigkurts = c(bigkurts, bigkurt)
    bigskews = c(bigskews, bigskew)
    bigmeds = c(bigmeds, median(bigsampleset))
    bigpval = bigresult$ks$p.value
    bigdval = bigresult$ks$statistic
    bigps = c(bigps, bigpval)
    bigds = c(bigds, bigdval)
    bothsampleset = sample(baseset, 12, replace = FALSE, prob = na.omit(bothmaxfit5))
    bothresult = ks.boot(log10(na.omit(Foss[,4])), log10(bothsampleset), nboots = 1000)
    bothmode = as.numeric(modetest(log10(bothsampleset), mod0 = 1)$statistic)
    bothkurt = as.numeric(anscombe.test(log10(bothsampleset))$statistic[1])
    bothskew = as.numeric(agostino.test(log10(bothsampleset), alternative = "two.sided")$statistic[1])
    bothmodes = c(bothmodes, bothmode)
    bothkurts = c(bothkurts, bothkurt)
    bothskews = c(bothskews, bothskew)
    bothmeds = c(bothmeds, median(bothsampleset))
    bothpval = bothresult$ks$p.value
    bothdval = bothresult$ks$statistic
    bothps = c(bothps, bothpval)
    bothds = c(bothds, bothdval)
    basesampleset = sample(baseset, 12, replace = FALSE)
    baseresult = ks.boot(log10(na.omit(Foss[,4])), log10(basesampleset), nboots = 1000)
    basemode = as.numeric(modetest(log10(basesampleset), mod0 = 1)$statistic)
    basekurt = as.numeric(anscombe.test(log10(basesampleset))$statistic[1])
    baseskew = as.numeric(agostino.test(log10(basesampleset), alternative = "two.sided")$statistic[1])
    basemodes = c(basemodes, basemode)
    basekurts = c(basekurts, basekurt)
    baseskews = c(baseskews, baseskew)
    basemeds = c(basemeds, median(basesampleset))
    basepval = baseresult$ks$p.value
    basedval = baseresult$ks$statistic
    baseps = c(baseps, basepval)
    baseds = c(baseds, basedval)
  }
  allsmallJOGmaxps = c(allsmallJOGmaxps, smallps)
  allsmallJOGmaxds = c(allsmallJOGmaxds, smallds)
  allsmallmodes = c(allsmallmodes, smallmodes)
  allsmallkurts = c(allsmallkurts, smallkurts)
  allsmallskews = c(allsmallskews, smallskews)
  allbigJOGmaxps = c(allbigJOGmaxps, bigps)
  allbigJOGmaxds = c(allbigJOGmaxds, bigds)
  allbigmodes = c(allbigmodes, bigmodes)
  allbigkurts = c(allbigkurts, bigkurts)
  allbigskews = c(allbigskews, bigskews)
  allbothJOGmaxps = c(allbothJOGmaxps, bothps)
  allbothJOGmaxds = c(allbothJOGmaxds, bothds)
  allbothmodes = c(allbothmodes, bothmodes)
  allbothkurts = c(allbothkurts, bothkurts)
  allbothskews = c(allbothskews, bothskews)
  allbaseJOGmaxps = c(allbaseJOGmaxps, baseps)
  allbaseJOGmaxds = c(allbaseJOGmaxds, baseds)
  allbasemodes = c(allbasemodes, basemodes)
  allbasekurts = c(allbasekurts, basekurts)
  allbaseskews = c(allbaseskews, baseskews)
  allsmallmeds = c(allsmallmeds, smallmeds)
  allbigmeds = c(allbigmeds, bigmeds)
  allbothmeds = c(allbothmeds, bothmeds)
  allbasemeds = c(allbasemeds, basemeds)
}

#now do PA
for(i in 1:29){
  baseset = na.omit(BM_mamms[,i])
  smallmaxfit5 = ifelse(BM_mamms[,i] < min(Foss[,11]), 0.05, 0.95)
  bigmaxfit5 = ifelse(BM_mamms[,i] > max(na.omit(Foss[,6])), 0.05, 0.95)
  bothmaxfit5 = ifelse(BM_mamms[,i] <  max(na.omit(Foss[,6])) & BM_mamms[,i] > min(Foss[,11]),0.95, 0.05)
  smallps = c()
  smallds = c()
  bigps = c()
  bigds = c()
  bothps = c()
  bothds = c()
  baseps = c()
  baseds = c()
  for(j in 1:100){
    smallsampleset = sample(baseset, 9, replace = FALSE, prob = na.omit(smallmaxfit5))
    smallresult = ks.boot(log10(na.omit(Foss[,10])), log10(smallsampleset), nboots = 1000)
    smallpval = smallresult$ks$p.value
    smalldval = smallresult$ks$statistic
    smallps = c(smallps, smallpval)
    smallds = c(smallds, smalldval)
    bigsampleset = sample(baseset, 9, replace = FALSE, prob = na.omit(bigmaxfit5))
    bigresult = ks.boot(log10(na.omit(Foss[,10])), log10(bigsampleset), nboots = 1000)
    bigpval = bigresult$ks$p.value
    bigdval = bigresult$ks$statistic
    bigps = c(bigps, bigpval)
    bigds = c(bigds, bigdval)
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
  allsmallPAmaxps = c(allsmallPAmaxps, smallps)
  allsmallPAmaxds = c(allsmallPAmaxds, smallds)
  allbigPAmaxps = c(allbigPAmaxps, bigps)
  allbigPAmaxds = c(allbigPAmaxds, bigds)
  allbothPAmaxps = c(allbothPAmaxps, bothps)
  allbothPAmaxds = c(allbothPAmaxds, bothds)
  allbasePAmaxps = c(allbasePAmaxps, baseps)
  allbasePAmaxds = c(allbasePAmaxds, baseds)
}

#now at sample size 20

for(i in 1:26){
  baseset = na.omit(BM_mamms20[,i])
  smallmaxfit5 = ifelse(BM_mamms20[,i] < min(Foss[,11]), 0.05, 0.95)
  bigmaxfit5 = ifelse(BM_mamms20[,i] > max(na.omit(Foss[,6])), 0.05, 0.95)
  bothmaxfit5 = ifelse(BM_mamms20[,i] <  max(na.omit(Foss[,6])) & BM_mamms20[,i] > min(Foss[,11]),0.95, 0.05)
  smallps = c()
  smallds = c()
  smallmodes = c()
  smallkurts = c()
  smallskews = c()
  smallmeds = c()
  bigps = c()
  bigds = c()
  bigmodes = c()
  bigkurts = c()
  bigskews = c()
  bigmeds = c()
  bothps = c()
  bothds = c()
  bothmodes = c()
  bothkurts = c()
  bothskews = c()
  bothmeds = c()
  baseps = c()
  baseds = c()
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
    bigsampleset = sample(baseset, 20, replace = FALSE, prob = na.omit(bigmaxfit5))
    bigresult = ks.boot(log10(na.omit(Foss[,4])), log10(bigsampleset), nboots = 1000)
    bigmode = as.numeric(modetest(log10(bigsampleset), mod0 = 1)$statistic)
    bigkurt = as.numeric(anscombe.test(log10(bigsampleset))$statistic[1])
    bigskew = as.numeric(agostino.test(log10(bigsampleset), alternative = "two.sided")$statistic[1])
    bigmodes = c(bigmodes, bigmode)
    bigkurts = c(bigkurts, bigkurt)
    bigskews = c(bigskews, bigskew)
    bigmeds = c(bigmeds, median(bigsampleset))
    bigpval = bigresult$ks$p.value
    bigdval = bigresult$ks$statistic
    bigps = c(bigps, bigpval)
    bigds = c(bigds, bigdval)
    bothsampleset = sample(baseset, 20, replace = FALSE, prob = na.omit(bothmaxfit5))
    bothresult = ks.boot(log10(na.omit(Foss[,4])), log10(bothsampleset), nboots = 1000)
    bothmode = as.numeric(modetest(log10(bothsampleset), mod0 = 1)$statistic)
    bothkurt = as.numeric(anscombe.test(log10(bothsampleset))$statistic[1])
    bothskew = as.numeric(agostino.test(log10(bothsampleset), alternative = "two.sided")$statistic[1])
    bothmodes = c(bothmodes, bothmode)
    bothkurts = c(bothkurts, bothkurt)
    bothskews = c(bothskews, bothskew)
    bothmeds = c(bothmeds, median(bothsampleset))
    bothpval = bothresult$ks$p.value
    bothdval = bothresult$ks$statistic
    bothps = c(bothps, bothpval)
    bothds = c(bothds, bothdval)
    basesampleset = sample(baseset, 20, replace = FALSE)
    baseresult = ks.boot(log10(na.omit(Foss[,4])), log10(basesampleset), nboots = 1000)
    basemode = as.numeric(modetest(log10(basesampleset), mod0 = 1)$statistic)
    basekurt = as.numeric(anscombe.test(log10(basesampleset))$statistic[1])
    baseskew = as.numeric(agostino.test(log10(basesampleset), alternative = "two.sided")$statistic[1])
    basemodes = c(basemodes, basemode)
    basekurts = c(basekurts, basekurt)
    baseskews = c(baseskews, baseskew)
    basemeds = c(basemeds, median(basesampleset))
    basepval = baseresult$ks$p.value
    basedval = baseresult$ks$statistic
    baseps = c(baseps, basepval)
    baseds = c(baseds, basedval)
  }
  allsmallJOGmaxps20 = c(allsmallJOGmaxps20, smallps)
  allsmallJOGmaxds20 = c(allsmallJOGmaxds20, smallds)
  allsmallmodes20 = c(allsmallmodes20, smallmodes)
  allsmallkurts20 = c(allsmallkurts20, smallkurts)
  allsmallskews20 = c(allsmallskews20, smallskews)
  allbigJOGmaxps20 = c(allbigJOGmaxps20, bigps)
  allbigJOGmaxds20 = c(allbigJOGmaxds20, bigds)
  allbigmodes20 = c(allbigmodes20, bigmodes)
  allbigkurts20 = c(allbigkurts20, bigkurts)
  allbigskews20 = c(allbigskews20, bigskews)
  allbothJOGmaxps20 = c(allbothJOGmaxps20, bothps)
  allbothJOGmaxds20 = c(allbothJOGmaxds20, bothds)
  allbothmodes20 = c(allbothmodes20, bothmodes)
  allbothkurts20 = c(allbothkurts20, bothkurts)
  allbothskews20 = c(allbothskews20, bothskews)
  allbaseJOGmaxps20 = c(allbaseJOGmaxps20, baseps)
  allbaseJOGmaxds20 = c(allbaseJOGmaxds20, baseds)
  allbasemodes20 = c(allbasemodes20, basemodes)
  allbasekurts20 = c(allbasekurts20, basekurts)
  allbaseskews20 = c(allbaseskews20, baseskews)
  allsmallmeds20 = c(allsmallmeds20, smallmeds)
  allbigmeds20 = c(allbigmeds20, bigmeds)
  allbothmeds20 = c(allbothmeds20, bothmeds)
  allbasemeds20 = c(allbasemeds20, basemeds)
}

#now do PA
for(i in 1:26){
  baseset = na.omit(BM_mamms20[,i])
  smallmaxfit5 = ifelse(BM_mamms20[,i] < min(Foss[,11]), 0.05, 0.95)
  bigmaxfit5 = ifelse(BM_mamms20[,i] > max(na.omit(Foss[,6])), 0.05, 0.95)
  bothmaxfit5 = ifelse(BM_mamms20[,i] <  max(na.omit(Foss[,6])) & BM_mamms20[,i] > min(Foss[,11]),0.95, 0.05)
  smallps = c()
  smallds = c()
  bigps = c()
  bigds = c()
  bothps = c()
  bothds = c()
  baseps = c()
  baseds = c()
  for(j in 1:100){
    smallsampleset = sample(baseset, 20, replace = FALSE, prob = na.omit(smallmaxfit5))
    smallresult = ks.boot(log10(na.omit(Foss[,10])), log10(smallsampleset), nboots = 1000)
    smallpval = smallresult$ks$p.value
    smalldval = smallresult$ks$statistic
    smallps = c(smallps, smallpval)
    smallds = c(smallds, smalldval)
    bigsampleset = sample(baseset, 20, replace = FALSE, prob = na.omit(bigmaxfit5))
    bigresult = ks.boot(log10(na.omit(Foss[,10])), log10(bigsampleset), nboots = 1000)
    bigpval = bigresult$ks$p.value
    bigdval = bigresult$ks$statistic
    bigps = c(bigps, bigpval)
    bigds = c(bigds, bigdval)
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
  allsmallPAmaxps20 = c(allsmallPAmaxps20, smallps)
  allsmallPAmaxds20 = c(allsmallPAmaxds20, smallds)
  allbigPAmaxps20 = c(allbigPAmaxps20, bigps)
  allbigPAmaxds20 = c(allbigPAmaxds20, bigds)
  allbothPAmaxps20 = c(allbothPAmaxps20, bothps)
  allbothPAmaxds20 = c(allbothPAmaxds20, bothds)
  allbasePAmaxps20 = c(allbasePAmaxps20, baseps)
  allbasePAmaxds20 = c(allbasePAmaxds20, baseds)
}

#now at sample size 40

for(i in 1:16){
  baseset = na.omit(BM_mamms40[,i])
  smallmaxfit5 = ifelse(BM_mamms40[,i] < min(Foss[,11]), 0.05, 0.95)
  bigmaxfit5 = ifelse(BM_mamms40[,i] > max(na.omit(Foss[,6])), 0.05, 0.95)
  bothmaxfit5 = ifelse(BM_mamms40[,i] <  max(na.omit(Foss[,6])) & BM_mamms40[,i] > min(Foss[,11]),0.95, 0.05)
  smallps = c()
  smallds = c()
  smallmodes = c()
  smallkurts = c()
  smallskews = c()
  smallmeds = c()
  bigps = c()
  bigds = c()
  bigmodes = c()
  bigkurts = c()
  bigskews = c()
  bigmeds = c()
  bothps = c()
  bothds = c()
  bothmodes = c()
  bothkurts = c()
  bothskews = c()
  bothmeds = c()
  baseps = c()
  baseds = c()
  basemodes = c()
  basekurts = c()
  baseskews = c()
  basemeds = c()
  for(j in 1:100){
    smallsampleset = sample(baseset, 40, replace = FALSE, prob = na.omit(smallmaxfit5))
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
    bigsampleset = sample(baseset, 40, replace = FALSE, prob = na.omit(bigmaxfit5))
    bigresult = ks.boot(log10(na.omit(Foss[,4])), log10(bigsampleset), nboots = 1000)
    bigmode = as.numeric(modetest(log10(bigsampleset), mod0 = 1)$statistic)
    bigkurt = as.numeric(anscombe.test(log10(bigsampleset))$statistic[1])
    bigskew = as.numeric(agostino.test(log10(bigsampleset), alternative = "two.sided")$statistic[1])
    bigmodes = c(bigmodes, bigmode)
    bigkurts = c(bigkurts, bigkurt)
    bigskews = c(bigskews, bigskew)
    bigmeds = c(bigmeds, median(bigsampleset))
    bigpval = bigresult$ks$p.value
    bigdval = bigresult$ks$statistic
    bigps = c(bigps, bigpval)
    bigds = c(bigds, bigdval)
    bothsampleset = sample(baseset, 40, replace = FALSE, prob = na.omit(bothmaxfit5))
    bothresult = ks.boot(log10(na.omit(Foss[,4])), log10(bothsampleset), nboots = 1000)
    bothmode = as.numeric(modetest(log10(bothsampleset), mod0 = 1)$statistic)
    bothkurt = as.numeric(anscombe.test(log10(bothsampleset))$statistic[1])
    bothskew = as.numeric(agostino.test(log10(bothsampleset), alternative = "two.sided")$statistic[1])
    bothmodes = c(bothmodes, bothmode)
    bothkurts = c(bothkurts, bothkurt)
    bothskews = c(bothskews, bothskew)
    bothmeds = c(bothmeds, median(bothsampleset))
    bothpval = bothresult$ks$p.value
    bothdval = bothresult$ks$statistic
    bothps = c(bothps, bothpval)
    bothds = c(bothds, bothdval)
    basesampleset = sample(baseset, 40, replace = FALSE)
    baseresult = ks.boot(log10(na.omit(Foss[,4])), log10(basesampleset), nboots = 1000)
    basemode = as.numeric(modetest(log10(basesampleset), mod0 = 1)$statistic)
    basekurt = as.numeric(anscombe.test(log10(basesampleset))$statistic[1])
    baseskew = as.numeric(agostino.test(log10(basesampleset), alternative = "two.sided")$statistic[1])
    basemodes = c(basemodes, basemode)
    basekurts = c(basekurts, basekurt)
    baseskews = c(baseskews, baseskew)
    basemeds = c(basemeds, median(basesampleset))
    basepval = baseresult$ks$p.value
    basedval = baseresult$ks$statistic
    baseps = c(baseps, basepval)
    baseds = c(baseds, basedval)
  }
  allsmallJOGmaxps40 = c(allsmallJOGmaxps40, smallps)
  allsmallJOGmaxds40 = c(allsmallJOGmaxds40, smallds)
  allsmallmodes40 = c(allsmallmodes40, smallmodes)
  allsmallkurts40 = c(allsmallkurts40, smallkurts)
  allsmallskews40 = c(allsmallskews40, smallskews)
  allbigJOGmaxps40 = c(allbigJOGmaxps40, bigps)
  allbigJOGmaxds40 = c(allbigJOGmaxds40, bigds)
  allbigmodes40 = c(allbigmodes40, bigmodes)
  allbigkurts40 = c(allbigkurts40, bigkurts)
  allbigskews40 = c(allbigskews40, bigskews)
  allbothJOGmaxps40 = c(allbothJOGmaxps40, bothps)
  allbothJOGmaxds40 = c(allbothJOGmaxds40, bothds)
  allbothmodes40 = c(allbothmodes40, bothmodes)
  allbothkurts40 = c(allbothkurts40, bothkurts)
  allbothskews40 = c(allbothskews40, bothskews)
  allbaseJOGmaxps40 = c(allbaseJOGmaxps40, baseps)
  allbaseJOGmaxds40 = c(allbaseJOGmaxds40, baseds)
  allbasemodes40 = c(allbasemodes40, basemodes)
  allbasekurts40 = c(allbasekurts40, basekurts)
  allbaseskews40 = c(allbaseskews40, baseskews)
  allsmallmeds40 = c(allsmallmeds40, smallmeds)
  allbigmeds40 = c(allbigmeds40, bigmeds)
  allbothmeds40 = c(allbothmeds40, bothmeds)
  allbasemeds40 = c(allbasemeds40, basemeds)
}

#now do PA
for(i in 1:16){
  baseset = na.omit(BM_mamms40[,i])
  smallmaxfit5 = ifelse(BM_mamms40[,i] < min(Foss[,11]), 0.05, 0.95)
  bigmaxfit5 = ifelse(BM_mamms40[,i] > max(na.omit(Foss[,6])), 0.05, 0.95)
  bothmaxfit5 = ifelse(BM_mamms40[,i] <  max(na.omit(Foss[,6])) & BM_mamms40[,i] > min(Foss[,11]),0.95, 0.05)
  smallps = c()
  smallds = c()
  bigps = c()
  bigds = c()
  bothps = c()
  bothds = c()
  baseps = c()
  baseds = c()
  for(j in 1:100){
    smallsampleset = sample(baseset, 40, replace = FALSE, prob = na.omit(smallmaxfit5))
    smallresult = ks.boot(log10(na.omit(Foss[,10])), log10(smallsampleset), nboots = 1000)
    smallpval = smallresult$ks$p.value
    smalldval = smallresult$ks$statistic
    smallps = c(smallps, smallpval)
    smallds = c(smallds, smalldval)
    bigsampleset = sample(baseset, 40, replace = FALSE, prob = na.omit(bigmaxfit5))
    bigresult = ks.boot(log10(na.omit(Foss[,10])), log10(bigsampleset), nboots = 1000)
    bigpval = bigresult$ks$p.value
    bigdval = bigresult$ks$statistic
    bigps = c(bigps, bigpval)
    bigds = c(bigds, bigdval)
    bothsampleset = sample(baseset, 40, replace = FALSE, prob = na.omit(bothmaxfit5))
    bothresult = ks.boot(log10(na.omit(Foss[,10])), log10(bothsampleset), nboots = 1000)
    bothpval = bothresult$ks$p.value
    bothdval = bothresult$ks$statistic
    bothps = c(bothps, bothpval)
    bothds = c(bothds, bothdval)
    basesampleset = sample(baseset, 40, replace = FALSE)
    baseresult = ks.boot(log10(na.omit(Foss[,10])), log10(basesampleset), nboots = 1000)
    basepval = baseresult$ks$p.value
    basedval = baseresult$ks$statistic
    baseps = c(baseps, basepval)
    baseds = c(baseds, basedval)
  }
  allsmallPAmaxps40 = c(allsmallPAmaxps40, smallps)
  allsmallPAmaxds40 = c(allsmallPAmaxds40, smallds)
  allbigPAmaxps40 = c(allbigPAmaxps40, bigps)
  allbigPAmaxds40 = c(allbigPAmaxds40, bigds)
  allbothPAmaxps40 = c(allbothPAmaxps40, bothps)
  allbothPAmaxds40 = c(allbothPAmaxds40, bothds)
  allbasePAmaxps40 = c(allbasePAmaxps40, baseps)
  allbasePAmaxds40 = c(allbasePAmaxds40, baseds)
}

#now at sample size 60

for(i in 1:10){
  baseset = na.omit(BM_mamms60[,i])
  smallmaxfit5 = ifelse(BM_mamms60[,i] < min(Foss[,11]), 0.05, 0.95)
  bigmaxfit5 = ifelse(BM_mamms60[,i] > max(na.omit(Foss[,6])), 0.05, 0.95)
  bothmaxfit5 = ifelse(BM_mamms60[,i] <  max(na.omit(Foss[,6])) & BM_mamms60[,i] > min(Foss[,11]),0.95, 0.05)
  smallps = c()
  smallds = c()
  smallmodes = c()
  smallkurts = c()
  smallskews = c()
  smallmeds = c()
  bigps = c()
  bigds = c()
  bigmodes = c()
  bigkurts = c()
  bigskews = c()
  bigmeds = c()
  bothps = c()
  bothds = c()
  bothmodes = c()
  bothkurts = c()
  bothskews = c()
  bothmeds = c()
  baseps = c()
  baseds = c()
  basemodes = c()
  basekurts = c()
  baseskews = c()
  basemeds = c()
  for(j in 1:100){
    smallsampleset = sample(baseset, 60, replace = FALSE, prob = na.omit(smallmaxfit5))
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
    bigsampleset = sample(baseset, 60, replace = FALSE, prob = na.omit(bigmaxfit5))
    bigresult = ks.boot(log10(na.omit(Foss[,4])), log10(bigsampleset), nboots = 1000)
    bigmode = as.numeric(modetest(log10(bigsampleset), mod0 = 1)$statistic)
    bigkurt = as.numeric(anscombe.test(log10(bigsampleset))$statistic[1])
    bigskew = as.numeric(agostino.test(log10(bigsampleset), alternative = "two.sided")$statistic[1])
    bigmodes = c(bigmodes, bigmode)
    bigkurts = c(bigkurts, bigkurt)
    bigskews = c(bigskews, bigskew)
    bigmeds = c(bigmeds, median(bigsampleset))
    bigpval = bigresult$ks$p.value
    bigdval = bigresult$ks$statistic
    bigps = c(bigps, bigpval)
    bigds = c(bigds, bigdval)
    bothsampleset = sample(baseset, 60, replace = FALSE, prob = na.omit(bothmaxfit5))
    bothresult = ks.boot(log10(na.omit(Foss[,4])), log10(bothsampleset), nboots = 1000)
    bothmode = as.numeric(modetest(log10(bothsampleset), mod0 = 1)$statistic)
    bothkurt = as.numeric(anscombe.test(log10(bothsampleset))$statistic[1])
    bothskew = as.numeric(agostino.test(log10(bothsampleset), alternative = "two.sided")$statistic[1])
    bothmodes = c(bothmodes, bothmode)
    bothkurts = c(bothkurts, bothkurt)
    bothskews = c(bothskews, bothskew)
    bothmeds = c(bothmeds, median(bothsampleset))
    bothpval = bothresult$ks$p.value
    bothdval = bothresult$ks$statistic
    bothps = c(bothps, bothpval)
    bothds = c(bothds, bothdval)
    basesampleset = sample(baseset, 60, replace = FALSE)
    baseresult = ks.boot(log10(na.omit(Foss[,4])), log10(basesampleset), nboots = 1000)
    basemode = as.numeric(modetest(log10(basesampleset), mod0 = 1)$statistic)
    basekurt = as.numeric(anscombe.test(log10(basesampleset))$statistic[1])
    baseskew = as.numeric(agostino.test(log10(basesampleset), alternative = "two.sided")$statistic[1])
    basemodes = c(basemodes, basemode)
    basekurts = c(basekurts, basekurt)
    baseskews = c(baseskews, baseskew)
    basemeds = c(basemeds, median(basesampleset))
    basepval = baseresult$ks$p.value
    basedval = baseresult$ks$statistic
    baseps = c(baseps, basepval)
    baseds = c(baseds, basedval)
  }
  allsmallJOGmaxps60 = c(allsmallJOGmaxps60, smallps)
  allsmallJOGmaxds60 = c(allsmallJOGmaxds60, smallds)
  allsmallmodes60 = c(allsmallmodes60, smallmodes)
  allsmallkurts60 = c(allsmallkurts60, smallkurts)
  allsmallskews60 = c(allsmallskews60, smallskews)
  allbigJOGmaxps60 = c(allbigJOGmaxps60, bigps)
  allbigJOGmaxds60 = c(allbigJOGmaxds60, bigds)
  allbigmodes60 = c(allbigmodes60, bigmodes)
  allbigkurts60 = c(allbigkurts60, bigkurts)
  allbigskews60 = c(allbigskews60, bigskews)
  allbothJOGmaxps60 = c(allbothJOGmaxps60, bothps)
  allbothJOGmaxds60 = c(allbothJOGmaxds60, bothds)
  allbothmodes60 = c(allbothmodes60, bothmodes)
  allbothkurts60 = c(allbothkurts60, bothkurts)
  allbothskews60 = c(allbothskews60, bothskews)
  allbaseJOGmaxps60 = c(allbaseJOGmaxps60, baseps)
  allbaseJOGmaxds60 = c(allbaseJOGmaxds60, baseds)
  allbasemodes60 = c(allbasemodes60, basemodes)
  allbasekurts60 = c(allbasekurts60, basekurts)
  allbaseskews60 = c(allbaseskews60, baseskews)
  allsmallmeds60 = c(allsmallmeds60, smallmeds)
  allbigmeds60 = c(allbigmeds60, bigmeds)
  allbothmeds60 = c(allbothmeds60, bothmeds)
  allbasemeds60 = c(allbasemeds60, basemeds)
}

#now do PA
for(i in 1:10){
  baseset = na.omit(BM_mamms60[,i])
  smallmaxfit5 = ifelse(BM_mamms60[,i] < min(Foss[,11]), 0.05, 0.95)
  bigmaxfit5 = ifelse(BM_mamms60[,i] > max(na.omit(Foss[,6])), 0.05, 0.95)
  bothmaxfit5 = ifelse(BM_mamms60[,i] <  max(na.omit(Foss[,6])) & BM_mamms60[,i] > min(Foss[,11]),0.95, 0.05)
  smallps = c()
  smallds = c()
  bigps = c()
  bigds = c()
  bothps = c()
  bothds = c()
  baseps = c()
  baseds = c()
  for(j in 1:100){
    smallsampleset = sample(baseset, 60, replace = FALSE, prob = na.omit(smallmaxfit5))
    smallresult = ks.boot(log10(na.omit(Foss[,10])), log10(smallsampleset), nboots = 1000)
    smallpval = smallresult$ks$p.value
    smalldval = smallresult$ks$statistic
    smallps = c(smallps, smallpval)
    smallds = c(smallds, smalldval)
    bigsampleset = sample(baseset, 60, replace = FALSE, prob = na.omit(bigmaxfit5))
    bigresult = ks.boot(log10(na.omit(Foss[,10])), log10(bigsampleset), nboots = 1000)
    bigpval = bigresult$ks$p.value
    bigdval = bigresult$ks$statistic
    bigps = c(bigps, bigpval)
    bigds = c(bigds, bigdval)
    bothsampleset = sample(baseset, 60, replace = FALSE, prob = na.omit(bothmaxfit5))
    bothresult = ks.boot(log10(na.omit(Foss[,10])), log10(bothsampleset), nboots = 1000)
    bothpval = bothresult$ks$p.value
    bothdval = bothresult$ks$statistic
    bothps = c(bothps, bothpval)
    bothds = c(bothds, bothdval)
    basesampleset = sample(baseset, 60, replace = FALSE)
    baseresult = ks.boot(log10(na.omit(Foss[,10])), log10(basesampleset), nboots = 1000)
    basepval = baseresult$ks$p.value
    basedval = baseresult$ks$statistic
    baseps = c(baseps, basepval)
    baseds = c(baseds, basedval)
  }
  allsmallPAmaxps60 = c(allsmallPAmaxps60, smallps)
  allsmallPAmaxds60 = c(allsmallPAmaxds60, smallds)
  allbigPAmaxps60 = c(allbigPAmaxps60, bigps)
  allbigPAmaxds60 = c(allbigPAmaxds60, bigds)
  allbothPAmaxps60 = c(allbothPAmaxps60, bothps)
  allbothPAmaxds60 = c(allbothPAmaxds60, bothds)
  allbasePAmaxps60 = c(allbasePAmaxps60, baseps)
  allbasePAmaxds60 = c(allbasePAmaxds60, baseds)
}

toprow = c("small-bias JOGmax", "small-bias JOGmax size20", "small-bias JOGmax size40", "small-bias JOGmax size60", "small-bias PAmax", "small-bias PAmax size 20", "small-bias PAmax size 40", "small-bias PAmax size 60",
             "large-bias JOGmax", "large-bias JOGmax size 20", "large-bias JOGmax size 40", "large-bias JOGmax size 60", "large-bias PAmax", "large-bias PAmax size 20", "large-bias PAmax size 40", "large-bias PAmax size 60",
             "edge-bias JOGmax", "edge-bias JOGmax size 20", "edge-bias JOGmax size 40", "edge-bias JOGmax size 60", "edge-bias PAmax", "edge-bias PAmax size 20", "edge-bias PAmax size 40", "edge-bias PAmax size 60", 
             "unbiased JOGmax", "unbiased JOGmax size 20", "unbiased JOGmax size 40", "unbiased JOGmax size 60", "unbiased PAmax", "unbiased PAmax size 20", "unbiased PAmax size 40", "unbiased PAmax size 60")

poutput = data.frame(matrix(nrow = 2900, ncol = 32))
poutput[,1] = allsmallJOGmaxps
poutput[1:2600,2] = allsmallJOGmaxps20
poutput[1:1600,3] = allsmallJOGmaxps40
poutput[1:1000,4] = allsmallJOGmaxps60
poutput[,5] = allsmallPAmaxps
poutput[1:2600,6] = allsmallPAmaxps20
poutput[1:1600,7] = allsmallPAmaxps40
poutput[1:1000,8] = allsmallPAmaxps60
poutput[,9] = allbigJOGmaxps
poutput[1:2600,10] = allbigJOGmaxps20
poutput[1:1600,11] = allbigJOGmaxps40
poutput[1:1000,12] = allbigJOGmaxps60
poutput[,13] = allbigPAmaxps
poutput[1:2600,14] = allbigPAmaxps20
poutput[1:1600,15] = allbigPAmaxps40
poutput[1:1000,16] = allbigPAmaxps60
poutput[,17] = allbothJOGmaxps
poutput[1:2600,18] = allbothJOGmaxps20
poutput[1:1600,19] = allbothJOGmaxps40
poutput[1:1000,20] = allbothJOGmaxps60
poutput[,21] = allbothPAmaxps
poutput[1:2600,22] = allbothPAmaxps20
poutput[1:1600,23] = allbothPAmaxps40
poutput[1:1000,24] = allbothPAmaxps60
poutput[,25] = allbaseJOGmaxps
poutput[1:2600,26] = allbaseJOGmaxps20
poutput[1:1600,27] = allbaseJOGmaxps40
poutput[1:1000,28] = allbaseJOGmaxps60
poutput[,29] = allbasePAmaxps
poutput[1:2600,30] = allbasePAmaxps20
poutput[1:1600,31] = allbasePAmaxps40
poutput[1:1000,32] = allbasePAmaxps60

colnames(poutput) = toprow


doutput = data.frame(matrix(nrow = 2900, ncol = 32))
doutput[,1] = allsmallJOGmaxds
doutput[1:2600,2] = allsmallJOGmaxds20
doutput[1:1600,3] = allsmallJOGmaxds40
doutput[1:1000,4] = allsmallJOGmaxds60
doutput[,5] = allsmallPAmaxds
doutput[1:2600,6] = allsmallPAmaxds20
doutput[1:1600,7] = allsmallPAmaxds40
doutput[1:1000,8] = allsmallPAmaxds60
doutput[,9] = allbigJOGmaxds
doutput[1:2600,10] = allbigJOGmaxds20
doutput[1:1600,11] = allbigJOGmaxds40
doutput[1:1000,12] = allbigJOGmaxds60
doutput[,13] = allbigPAmaxds
doutput[1:2600,14] = allbigPAmaxds20
doutput[1:1600,15] = allbigPAmaxds40
doutput[1:1000,16] = allbigPAmaxds60
doutput[,17] = allbothJOGmaxds
doutput[1:2600,18] = allbothJOGmaxds20
doutput[1:1600,19] = allbothJOGmaxds40
doutput[1:1000,20] = allbothJOGmaxds60
doutput[,21] = allbothPAmaxds
doutput[1:2600,22] = allbothPAmaxds20
doutput[1:1600,23] = allbothPAmaxds40
doutput[1:1000,24] = allbothPAmaxds60
doutput[,25] = allbaseJOGmaxds
doutput[1:2600,26] = allbaseJOGmaxds20
doutput[1:1600,27] = allbaseJOGmaxds40
doutput[1:1000,28] = allbaseJOGmaxds60
doutput[,29] = allbasePAmaxds
doutput[1:2600,30] = allbasePAmaxds20
doutput[1:1600,31] = allbasePAmaxds40
doutput[1:1000,32] = allbasePAmaxds60

colnames(doutput) = toprow

toprow2 = c("small-bias JOGmax", "small-bias JOGmax size20", "small-bias JOGmax size40", "small-bias JOGmax size60",
           "large-bias JOGmax", "large-bias JOGmax size 20", "large-bias JOGmax size 40", "large-bias JOGmax size 60", 
           "edge-bias JOGmax", "edge-bias JOGmax size 20", "edge-bias JOGmax size 40", "edge-bias JOGmax size 60",  
           "unbiased JOGmax", "unbiased JOGmax size 20", "unbiased JOGmax size 40", "unbiased JOGmax size 60")

moutput = data.frame(matrix(nrow = 2900, ncol = 16))
moutput[,1] = allsmallmodes
moutput[1:2600,2] = allsmallmodes20
moutput[1:1600,3] = allsmallmodes40
moutput[1:1000,4] = allsmallmodes60
moutput[,5] = allbigmodes
moutput[1:2600,6] = allbigmodes20
moutput[1:1600,7] = allbigmodes40
moutput[1:1000,8] = allbigmodes60
moutput[,9] = allbothmodes
moutput[1:2600,10] = allbothmodes20
moutput[1:1600,11] = allbothmodes40
moutput[1:1000,12] = allbothmodes60
moutput[,13] = allbasemodes
moutput[1:2600,14] = allbasemodes20
moutput[1:1600,15] = allbasemodes40
moutput[1:1000,16] = allbasemodes60[1:1000]

colnames(moutput) = toprow2

koutput = data.frame(matrix(nrow = 2900, ncol = 16))
koutput[,1] = allsmallkurts
koutput[1:2600,2] = allsmallkurts20
koutput[1:1600,3] = allsmallkurts40
koutput[1:1000,4] = allsmallkurts60
koutput[,5] = allbigkurts
koutput[1:2600,6] = allbigkurts20
koutput[1:1600,7] = allbigkurts40
koutput[1:1000,8] = allbigkurts60
koutput[,9] = allbothkurts
koutput[1:2600,10] = allbothkurts20
koutput[1:1600,11] = allbothkurts40
koutput[1:1000,12] = allbothkurts60
koutput[,13] = allbasekurts
koutput[1:2600,14] = allbasekurts20
koutput[1:1600,15] = allbasekurts40
koutput[1:1000,16] = allbasekurts60[1:1000]

colnames(koutput) = toprow2

soutput = data.frame(matrix(nrow = 2900, ncol = 16))
soutput[,1] = allsmallskews
soutput[1:2600,2] = allsmallskews20
soutput[1:1600,3] = allsmallskews40
soutput[1:1000,4] = allsmallskews60
soutput[,5] = allbigskews
soutput[1:2600,6] = allbigskews20
soutput[1:1600,7] = allbigskews40
soutput[1:1000,8] = allbigskews60
soutput[,9] = allbothskews
soutput[1:2600,10] = allbothskews20
soutput[1:1600,11] = allbothskews40
soutput[1:1000,12] = allbothskews60
soutput[,13] = allbaseskews
soutput[1:2600,14] = allbaseskews20
soutput[1:1600,15] = allbaseskews40
soutput[1:1000,16] = allbaseskews60[1:1000]

colnames(soutput) = toprow2

aoutput = data.frame(matrix(nrow = 2900, ncol = 16))
aoutput[,1] = allsmallmeds
aoutput[1:2600,2] = allsmallmeds20
aoutput[1:1600,3] = allsmallmeds40
aoutput[1:1000,4] = allsmallmeds60
aoutput[,5] = allbigmeds
aoutput[1:2600,6] = allbigmeds20
aoutput[1:1600,7] = allbigmeds40
aoutput[1:1000,8] = allbigmeds60
aoutput[,9] = allbothmeds
aoutput[1:2600,10] = allbothmeds20
aoutput[1:1600,11] = allbothmeds40
aoutput[1:1000,12] = allbothmeds60
aoutput[,13] = allbasemeds
aoutput[1:2600,14] = allbasemeds20
aoutput[1:1600,15] = allbasemeds40
aoutput[1:1000,16] = allbasemeds60[1:1000]

colnames(aoutput) = toprow2

write.csv(poutput, "C:\\Wirc\\Projekt3\\Bigness\\subsets\\pvals_fullTaph_pgls.csv")
write.csv(doutput, "C:\\Wirc\\Projekt3\\Bigness\\subsets\\dvals_fullTaph_pgls.csv")
write.csv(moutput, "C:\\Wirc\\Projekt3\\Bigness\\subsets\\modes_fullTaph_pgls.csv")
write.csv(koutput, "C:\\Wirc\\Projekt3\\Bigness\\subsets\\kurts_fullTaph_pgls.csv")
write.csv(soutput, "C:\\Wirc\\Projekt3\\Bigness\\subsets\\skews_fullTaph_pgls.csv")
write.csv(aoutput, "C:\\Wirc\\Projekt3\\Bigness\\subsets\\medians_fullTaph_pgls.csv")
