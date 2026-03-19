library(moments)
library(multimode)
library(Matching)

BM_mamms = read.csv("C:\\Wirc\\Projekt3\\Bigness\\subsets\\TetBMs_justmammsTaphready.csv", header = TRUE)
BM_mamms20 = read.csv("C:\\Wirc\\Projekt3\\Bigness\\subsets\\TetBMs_justmammsTaphready_size20.csv", header = TRUE)
BM_mamms40 = read.csv("C:\\Wirc\\Projekt3\\Bigness\\subsets\\TetBMs_justmammsTaphready_size40.csv", header = TRUE)

sites = 18
sites20 = 10
sites40 = 4


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

allsmallmodes40 = c()
allsmallkurts40 = c()
allsmallskews40 = c()
allsmallmeds40 = c()
allbothmodes40 = c()
allbothkurts40 = c()
allbothskews40 = c()
allbothmeds40 = c()
allbasemodes40 = c()
allbasekurts40 = c()
allbaseskews40 = c()
allbasemeds40 = c()

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

allsmallJOGmaxps40 = c()
allsmallJOGmaxds40 = c()

allbothPAmaxps40 = c()
allbothPAmaxds40 = c()
allbasePAmaxps40 = c()
allbasePAmaxds40 = c()


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

#now at sample size 40

for(i in 1:sites40){
  baseset = na.omit(BM_mamms40[,i])
  smallmaxfit5 = ifelse(BM_mamms40[,i] < min(Foss[,11]), 0.05, 0.95)
  bothmaxfit5 = ifelse(BM_mamms40[,i] <  max(na.omit(Foss[,6])) & BM_mamms40[,i] > min(Foss[,11]),0.95, 0.05)
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
    bothsampleset = sample(baseset, 40, replace = FALSE, prob = na.omit(bothmaxfit5))
    bothmode = as.numeric(modetest(log10(bothsampleset), mod0 = 1)$statistic)
    bothkurt = as.numeric(anscombe.test(log10(bothsampleset))$statistic[1])
    bothskew = as.numeric(agostino.test(log10(bothsampleset), alternative = "two.sided")$statistic[1])
    bothmodes = c(bothmodes, bothmode)
    bothkurts = c(bothkurts, bothkurt)
    bothskews = c(bothskews, bothskew)
    bothmeds = c(bothmeds, median(bothsampleset))
    basesampleset = sample(baseset, 40, replace = FALSE)
    basemode = as.numeric(modetest(log10(basesampleset), mod0 = 1)$statistic)
    basekurt = as.numeric(anscombe.test(log10(basesampleset))$statistic[1])
    baseskew = as.numeric(agostino.test(log10(basesampleset), alternative = "two.sided")$statistic[1])
    basemodes = c(basemodes, basemode)
    basekurts = c(basekurts, basekurt)
    baseskews = c(baseskews, baseskew)
    basemeds = c(basemeds, median(basesampleset))
  }
  allsmallJOGmaxps40 = c(allsmallJOGmaxps40, smallps)
  allsmallJOGmaxds40 = c(allsmallJOGmaxds40, smallds)
  allsmallmodes40 = c(allsmallmodes40, smallmodes)
  allsmallkurts40 = c(allsmallkurts40, smallkurts)
  allsmallskews40 = c(allsmallskews40, smallskews)
  allbothmodes40 = c(allbothmodes40, bothmodes)
  allbothkurts40 = c(allbothkurts40, bothkurts)
  allbothskews40 = c(allbothskews40, bothskews)
  allbasemodes40 = c(allbasemodes40, basemodes)
  allbasekurts40 = c(allbasekurts40, basekurts)
  allbaseskews40 = c(allbaseskews40, baseskews)
  allsmallmeds40 = c(allsmallmeds40, smallmeds)
  allbothmeds40 = c(allbothmeds40, bothmeds)
  allbasemeds40 = c(allbasemeds40, basemeds)
}

#now do PA
for(i in 1:sites40){
  baseset = na.omit(BM_mamms40[,i])
  bothmaxfit5 = ifelse(BM_mamms40[,i] <  max(na.omit(Foss[,6])) & BM_mamms40[,i] > min(Foss[,11]),0.95, 0.05)
  bothps = c()
  bothds = c()
  baseps = c()
  baseds = c()
  for(j in 1:100){
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
  allbothPAmaxps40 = c(allbothPAmaxps40, bothps)
  allbothPAmaxds40 = c(allbothPAmaxds40, bothds)
  allbasePAmaxps40 = c(allbasePAmaxps40, baseps)
  allbasePAmaxds40 = c(allbasePAmaxds40, baseds)
}


toprow = c("small-bias JOGmax", "small-bias JOGmax size20", "small-bias JOGmax size40", 
             "edge-bias PAmax", "edge-bias PAmax size 20", "edge-bias PAmax size 40",  
             "unbiased PAmax", "unbiased PAmax size 20", "unbiased PAmax size 40")

poutput = data.frame(matrix(nrow = (sites*100), ncol = 9))
poutput[,1] = allsmallJOGmaxps
poutput[1:(sites20*100),2] = allsmallJOGmaxps20
poutput[1:(sites40*100),3] = allsmallJOGmaxps40
poutput[,4] = allbothPAmaxps
poutput[1:(sites20*100),5] = allbothPAmaxps20
poutput[1:(sites40*100),6] = allbothPAmaxps40
poutput[,7] = allbasePAmaxps
poutput[1:(sites20*100),8] = allbasePAmaxps20
poutput[1:(sites40*100),9] = allbasePAmaxps40

colnames(poutput) = toprow


doutput = data.frame(matrix(nrow = (sites*100), ncol = 9))
doutput[,1] = allsmallJOGmaxds
doutput[1:(sites20*100),2] = allsmallJOGmaxds20
doutput[1:(sites40*100),3] = allsmallJOGmaxds40
doutput[,4] = allbothPAmaxds
doutput[1:(sites20*100),5] = allbothPAmaxds20
doutput[1:(sites40*100),6] = allbothPAmaxds40
doutput[,7] = allbasePAmaxds
doutput[1:(sites20*100),8] = allbasePAmaxds20
doutput[1:(sites40*100),9] = allbasePAmaxds40

colnames(doutput) = toprow

toprow2 = c("small-bias", "small-bias size 20", "small-bias size 40",
           "edge-bias", "edge-bias size 20", "edge-bias size 40",  
           "unbiased", "unbiased size 20", "unbiased size 40")

moutput = data.frame(matrix(nrow = (sites*100), ncol = 9))
moutput[,1] = allsmallmodes
moutput[1:(sites20*100),2] = allsmallmodes20
moutput[1:(sites40*100),3] = allsmallmodes40
moutput[,4] = allbothmodes
moutput[1:(sites20*100),5] = allbothmodes20
moutput[1:(sites40*100),6] = allbothmodes40
moutput[,7] = allbasemodes
moutput[1:(sites20*100),8] = allbasemodes20
moutput[1:(sites40*100),9] = allbasemodes40


colnames(moutput) = toprow2

koutput = data.frame(matrix(nrow = (sites*100), ncol = 9))
koutput[,1] = allsmallkurts
koutput[1:(sites20*100),2] = allsmallkurts20
koutput[1:(sites40*100),3] = allsmallkurts40
koutput[,4] = allbothkurts
koutput[1:(sites20*100),5] = allbothkurts20
koutput[1:(sites40*100),6] = allbothkurts40
koutput[,7] = allbasekurts
koutput[1:(sites20*100),8] = allbasekurts20
koutput[1:(sites40*100),9] = allbasekurts40

colnames(koutput) = toprow2

soutput = data.frame(matrix(nrow = (sites*100), ncol = 9))
soutput[,1] = allsmallskews
soutput[1:(sites20*100),2] = allsmallskews20
soutput[1:(sites40*100),3] = allsmallskews40
soutput[,4] = allbothskews
soutput[1:(sites20*100),5] = allbothskews20
soutput[1:(sites40*100),6] = allbothskews40
soutput[,7] = allbaseskews
soutput[1:(sites20*100),8] = allbaseskews20
soutput[1:(sites40*100),9] = allbaseskews40

colnames(soutput) = toprow2

aoutput = data.frame(matrix(nrow = (sites*100), ncol = 9))
aoutput[,1] = allsmallmeds
aoutput[1:(sites20*100),2] = allsmallmeds20
aoutput[1:(sites40*100),3] = allsmallmeds40
aoutput[,4] = allbothmeds
aoutput[1:(sites20*100),5] = allbothmeds20
aoutput[1:(sites40*100),6] = allbothmeds40
aoutput[,7] = allbasemeds
aoutput[1:(sites20*100),8] = allbasemeds20
aoutput[1:(sites40*100),9] = allbasemeds40

colnames(aoutput) = toprow2

write.csv(poutput, "C:\\Wirc\\Projekt3\\Bigness\\subsets\\pvals_justmammsTaph_pgls.csv")
write.csv(doutput, "C:\\Wirc\\Projekt3\\Bigness\\subsets\\dvals_justmammsTaph_pgls.csv")
write.csv(moutput, "C:\\Wirc\\Projekt3\\Bigness\\subsets\\modes_justmammsTaph_pgls.csv")
write.csv(koutput, "C:\\Wirc\\Projekt3\\Bigness\\subsets\\kurts_justmammsTaph_pgls.csv")
write.csv(soutput, "C:\\Wirc\\Projekt3\\Bigness\\subsets\\skews_justmammsTaph_pgls.csv")
write.csv(aoutput, "C:\\Wirc\\Projekt3\\Bigness\\subsets\\medians_justmammsTaph_pgls.csv")
