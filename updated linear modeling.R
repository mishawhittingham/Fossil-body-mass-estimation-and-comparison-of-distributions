fullset = read.csv("C:\\Wirc\\Projekt3\\Bigness\\Do DVCs even work\\testing\\Full set.csv")
fullset$Clade <- as.factor(fullset$Clade)
library(MLmetrics)
library(AICcmodavg)

comparison = data.frame(matrix(nrow = 3, ncol = 10))
rownames(comparison) = c("meanCSA", "minCSA", "maxCSA")
colnames(comparison) = c("intercept", "slope", "r-squared", "SEE", "AICc", "Akaike Weight", "Correction Factor", "PPE", "PPE CI Lower", "PPE CI Upper")

BMs = fullset$BM.log
meanCSA = fullset$Mean.CSA
minCSA = fullset$Min.CSA
maxCSA = fullset$Max.CSA
rawBMs = fullset$BM.g


meanCSA_lm = lm(BMs~meanCSA)
sum_meanCSA = summary.lm(meanCSA_lm)
formula_meanCSA = meanCSA_lm$coefficients
comparison[1,1] = formula_meanCSA[1]
comparison[1,2] = formula_meanCSA[2]
comparison[1,3] = sum_meanCSA$adj.r.squared
comparison[1,4] = sum_meanCSA$sigma
resids_meanCSA = abs(residuals(meanCSA_lm))
detrans_resids_meanCSA = 10^(residuals(meanCSA_lm))
comparison[1,7] = mean(detrans_resids_meanCSA)
comparison[1,8] = mean(resids_meanCSA/(BMs))
comparison[1,9] = t.test(resids_meanCSA/(BMs))$conf.int[1]
comparison[1,10] = t.test(resids_meanCSA/BMs)$conf.int[2]

trymean=10^(predict(meanCSA_lm, data.frame(meanCSA = fullset$Mean.CSA)))
unlogppemean =  median(abs(rawBMs-trymean)/rawBMs)

minCSA_lm = lm(BMs~minCSA)
sum_minCSA = summary.lm(minCSA_lm)
formula_minCSA = minCSA_lm$coefficients
comparison[2,1] = formula_minCSA[1]
comparison[2,2] = formula_minCSA[2]
comparison[2,3] = sum_minCSA$adj.r.squared
comparison[2,4] = sum_minCSA$sigma
resids_minCSA = abs(residuals(minCSA_lm))
detrans_resids_minCSA = 10^(residuals(minCSA_lm))
comparison[2,7] = mean(detrans_resids_minCSA)
comparison[2,8] = mean(resids_minCSA/(BMs))
comparison[2,9] = t.test(resids_minCSA/(BMs))$conf.int[1]
comparison[2,10] = t.test(resids_minCSA/(BMs))$conf.int[2]

trymin=10^(predict(minCSA_lm, data.frame(minCSA = fullset$Min.CSA)))
unlogppemin =  median(abs(rawBMs-trymin)/rawBMs)


maxCSA_lm = lm(BMs~maxCSA)
sum_maxCSA = summary.lm(maxCSA_lm)
formula_maxCSA = maxCSA_lm$coefficients
comparison[3,1] = formula_maxCSA[1]
comparison[3,2] = formula_maxCSA[2]
comparison[3,3] = sum_maxCSA$adj.r.squared
comparison[3,4] = sum_maxCSA$sigma
resids_maxCSA = abs(residuals(maxCSA_lm))
detrans_resids_maxCSA = 10^(residuals(maxCSA_lm))
comparison[3,7] = mean(detrans_resids_maxCSA)
comparison[3,8] = mean(resids_maxCSA/(BMs))
comparison[3,9] = t.test(resids_maxCSA/(BMs))$conf.int[1]
comparison[3,10] = t.test(resids_maxCSA/(BMs))$conf.int[2]

trymax=10^(predict(maxCSA_lm, data.frame(maxCSA = fullset$Max.CSA)))
unlogppemax =  median(abs(rawBMs-trymax)/rawBMs)

write.csv(comparison, "C:\\Wirc\\Projekt3\\Bigness\\Do DVCs even work\\testing\\OLS results summary.csv")

header = c("Clade", "resids")
residsmean = data.frame(matrix(nrow=82,ncol=2))
residsmean[,1] = fullset$Clade
residsmean[,2] = residuals(meanCSA_lm)
residsmax = data.frame(matrix(nrow=82,ncol=2))
residsmax[,1] = fullset$Clade
residsmax[,2] = residuals(maxCSA_lm)
residsmin = data.frame(matrix(nrow=82,ncol=2))
residsmin[,1] = fullset$Clade
residsmin[,2] = residuals(minCSA_lm)
colnames(residsmean) = header
colnames(residsmax) = header
colnames(residsmin) = header

pdf(file="C:\\Wirc\\Projekt3\\Bigness\\Do DVCs even work\\testing\\lm residuals.pdf", width = 10, height = 10)
par(mfrow = c(3,3))

ggplot(residsmean, aes(group=Clade, Clade, resids))+
  geom_boxplot()+
  stat_summary(fun.y=mean, geom= 'point', size=4, color = 'red')+
  labs(title = "Difference in residuals among clades for mean OSA vs Body Mass")

ggplot(residsmax, aes(group=Clade, Clade, resids))+
  geom_boxplot()+
  stat_summary(fun.y=mean, geom= 'point', size=4, color = 'red')+
  labs(title = "Difference in residuals among clades for max OSA vs Body Mass")

ggplot(residsmin, aes(group=Clade, Clade, resids))+
  geom_boxplot()+
  stat_summary(fun.y=mean, geom= 'point', size=4, color = 'red')+
  labs(title = "Difference in residuals among clades for min OSA vs Body Mass")

n = dev.off()

romerdata = read.csv("C:\\Wirc\\Projekt3\\Bigness\\Do DVCs even work\\testing\\romerverts.csv", header = TRUE)
mydata = read.csv("C:\\Wirc\\Projekt3\\Bigness\\Do DVCs even work\\mydata.csv", header = TRUE)

romermaxcorrect = predict(maxCSA_lm, data.frame(maxCSA = romerdata$MaxCSA), interval = "confidence")
detransromer = (10^(romermaxcorrect[,1]))*comparison[3,7]
romerrangeU = 10^(romermaxcorrect[,1]*(1+comparison[3,8]))*comparison[3,7]
romerrangeL = 10^(romermaxcorrect[,1]*(1-comparison[3,8]))*comparison[3,7]
row.names(romermaxcorrect) = romerdata[,1]
romermaxcorrect = cbind(detransromer, romerrangeL, romerrangeU, romerdata[,2])

mymaxcorrect = predict(maxCSA_lm, data.frame(maxCSA = mydata$MaxCSA), interval = "confidence")
mydetrans = (10^(mymaxcorrect[,1]))*comparison[3,7]
myrangeU = 10^(mymaxcorrect[,1]*(1+comparison[3,8]))*comparison[3,7]
myrangeL = 10^(mymaxcorrect[,1]*(1-comparison[3,8]))*comparison[3,7]
row.names(mymaxcorrect) = mydata[,1]
mymaxcorrect = cbind(mydetrans, myrangeL, myrangeU, mydata$Locality, mydata$Name)

write.csv(romermaxcorrect, "C:\\Wirc\\Projekt3\\Bigness\\Do DVCs even work\\testing\\romermaxcorrected update2.csv")
write.csv(mymaxcorrect, "C:\\Wirc\\Projekt3\\Bigness\\Do DVCs even work\\new BM ests (maxCSA) update2.csv")


