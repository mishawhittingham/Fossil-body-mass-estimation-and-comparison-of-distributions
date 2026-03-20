fullset = read.csv("C:\\...\\Full set.csv")
fullset$Clade <- as.factor(fullset$Clade)
library(MLmetrics)
library(AICcmodavg)

comparison = data.frame(matrix(nrow = 3, ncol = 10))
rownames(comparison) = c("meanOSA", "minOSA", "maxOSA")
colnames(comparison) = c("intercept", "slope", "r-squared", "SEE", "AICc", "Akaike Weight", "Correction Factor", "PPE", "PPE CI Lower", "PPE CI Upper")

BMs = fullset$BM.log
meanOSA = fullset$Mean.OSA
minOSA = fullset$Min.OSA
maxOSA = fullset$Max.OSA
rawBMs = fullset$BM.g


meanOSA_lm = lm(BMs~meanOSA)
sum_meanOSA = summary.lm(meanOSA_lm)
formula_meanOSA = meanOSA_lm$coefficients
comparison[1,1] = formula_meanOSA[1]
comparison[1,2] = formula_meanOSA[2]
comparison[1,3] = sum_meanOSA$adj.r.squared
comparison[1,4] = sum_meanOSA$sigma
resids_meanOSA = abs(residuals(meanOSA_lm))
detrans_resids_meanOSA = 10^(residuals(meanOSA_lm))
comparison[1,7] = mean(detrans_resids_meanOSA)
comparison[1,8] = mean(resids_meanOSA/(BMs))
comparison[1,9] = t.test(resids_meanOSA/(BMs))$conf.int[1]
comparison[1,10] = t.test(resids_meanOSA/BMs)$conf.int[2]

trymean=10^(predict(meanOSA_lm, data.frame(meanOSA = fullset$Mean.OSA)))
unlogppemean =  median(abs(rawBMs-trymean)/rawBMs)

minOSA_lm = lm(BMs~minOSA)
sum_minOSA = summary.lm(minOSA_lm)
formula_minOSA = minOSA_lm$coefficients
comparison[2,1] = formula_minOSA[1]
comparison[2,2] = formula_minOSA[2]
comparison[2,3] = sum_minOSA$adj.r.squared
comparison[2,4] = sum_minOSA$sigma
resids_minOSA = abs(residuals(minOSA_lm))
detrans_resids_minOSA = 10^(residuals(minOSA_lm))
comparison[2,7] = mean(detrans_resids_minOSA)
comparison[2,8] = mean(resids_minOSA/(BMs))
comparison[2,9] = t.test(resids_minOSA/(BMs))$conf.int[1]
comparison[2,10] = t.test(resids_minOSA/(BMs))$conf.int[2]

trymin=10^(predict(minOSA_lm, data.frame(minOSA = fullset$Min.OSA)))
unlogppemin =  median(abs(rawBMs-trymin)/rawBMs)


maxOSA_lm = lm(BMs~maxOSA)
sum_maxOSA = summary.lm(maxOSA_lm)
formula_maxOSA = maxOSA_lm$coefficients
comparison[3,1] = formula_maxOSA[1]
comparison[3,2] = formula_maxOSA[2]
comparison[3,3] = sum_maxOSA$adj.r.squared
comparison[3,4] = sum_maxOSA$sigma
resids_maxOSA = abs(residuals(maxOSA_lm))
detrans_resids_maxOSA = 10^(residuals(maxOSA_lm))
comparison[3,7] = mean(detrans_resids_maxOSA)
comparison[3,8] = mean(resids_maxOSA/(BMs))
comparison[3,9] = t.test(resids_maxOSA/(BMs))$conf.int[1]
comparison[3,10] = t.test(resids_maxOSA/(BMs))$conf.int[2]

trymax=10^(predict(maxOSA_lm, data.frame(maxOSA = fullset$Max.OSA)))
unlogppemax =  median(abs(rawBMs-trymax)/rawBMs)

write.csv(comparison, "C:\\...\\OLS results summary.csv")

header = c("Clade", "resids")
residsmean = data.frame(matrix(nrow=82,ncol=2))
residsmean[,1] = fullset$Clade
residsmean[,2] = residuals(meanOSA_lm)
residsmax = data.frame(matrix(nrow=82,ncol=2))
residsmax[,1] = fullset$Clade
residsmax[,2] = residuals(maxOSA_lm)
residsmin = data.frame(matrix(nrow=82,ncol=2))
residsmin[,1] = fullset$Clade
residsmin[,2] = residuals(minOSA_lm)
colnames(residsmean) = header
colnames(residsmax) = header
colnames(residsmin) = header

romerdata = read.csv("C:\\...\\romerverts.csv", header = TRUE)
mydata = read.csv("C:\\...\\mydata.csv", header = TRUE)

romermaxcorrect = predict(maxOSA_lm, data.frame(maxOSA = romerdata$MaxOSA), interval = "confidence")
detransromer = (10^(romermaxcorrect[,1]))*comparison[3,7]
romerrangeU = 10^(romermaxcorrect[,1]*(1+comparison[3,8]))*comparison[3,7]
romerrangeL = 10^(romermaxcorrect[,1]*(1-comparison[3,8]))*comparison[3,7]
row.names(romermaxcorrect) = romerdata[,1]
romermaxcorrect = cbind(detransromer, romerrangeL, romerrangeU, romerdata[,2])

mymaxcorrect = predict(maxOSA_lm, data.frame(maxOSA = mydata$MaxOSA), interval = "confidence")
mydetrans = (10^(mymaxcorrect[,1]))*comparison[3,7]
myrangeU = 10^(mymaxcorrect[,1]*(1+comparison[3,8]))*comparison[3,7]
myrangeL = 10^(mymaxcorrect[,1]*(1-comparison[3,8]))*comparison[3,7]
row.names(mymaxcorrect) = mydata[,1]
mymaxcorrect = cbind(mydetrans, myrangeL, myrangeU, mydata$Locality, mydata$Name)

write.csv(romermaxcorrect, "C:\\...\\romermaxcorrected update2.csv")
write.csv(mymaxcorrect, "C:\\...\\Do DVCs even work\\new BM ests (maxOSA).csv")


