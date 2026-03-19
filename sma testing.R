library(smatr)
library(ggplot2)
fullset = read.csv("C:\\Wirc\\Projekt3\\Bigness\\Do DVCs even work\\testing\\Full set.csv")
fullset$Clade = as.factor(fullset$Clade)
fullset$Terrestriality= as.factor(fullset$Terrestriality)
fullset$BM.Category = as.factor(fullset$BM.Category)

meanCSA_sma = sma(fullset$BM.g~fullset$Mean.CSA, data = fullset, method = "SMA", alpha = 0.05)
minCSA_sma = sma(fullset$BM.g~fullset$Min.CSA, data = fullset, method = "SMA", alpha = 0.05)
maxCSA_sma = sma(fullset$BM.g~fullset$Max.CSA, data = fullset, method = "SMA", alpha = 0.05)


SMAcomparison = data.frame(matrix(nrow = 33, ncol = 17))
SMAcomparison[,1] = rep(c("meanCSA", "minCSA", "maxCSA"), each = 11)
colnames(SMAcomparison) = c("test variable","subset","n","r-squared", "p-value", "slope","slope CI lower", "slope CI upper", "intercept", "intercept CI lower", "intercept CI upper", "clade slopes p", "terrestriality slopes p", "BM category slopes p", "clade intercepts p", "terrestriality intercepts p", "BM category intercepts p")


pdf(file="C:\\Wirc\\Projekt3\\Bigness\\Do DVCs even work\\testing\\SMAs.pdf", width = 10, height = 10)
par(mfrow = c(3,3))

for(i in 1:3){
  result = sma(fullset$BM.g~fullset[,(9+i)], data = fullset, method = "SMA", alpha = 0.05)
  clades = sma(fullset$BM.g~fullset[,(9+i)]*fullset$Clade, data = fullset, method = "SMA", alpha = 0.05)
  terrestrials = sma(fullset$BM.g~fullset[,(9+i)]*fullset$Terrestriality, data = fullset, method = "SMA", alpha = 0.05)
  BMcategories = sma(fullset$BM.g~fullset[,(9+i)]*fullset$BM.Category, data = fullset, method = "SMA", alpha = 0.05)
  clades_e = sma(fullset$BM.g~fullset[,(9+i)]+fullset$Clade, data = fullset, method = "SMA", type = "elevation", alpha = 0.05)
  terrestrials_e = sma(fullset$BM.g~fullset[,(9+i)]+fullset$Terrestriality, data = fullset, method = "SMA", type = "elevation", alpha = 0.05)
  BMcategories_e = sma(fullset$BM.g~fullset[,(9+i)]+fullset$BM.Category, data = fullset, method = "SMA", type = "elevation", alpha = 0.05)
  
  SMAcomparison[(11*(i-1))+1,2:11] = result$groupsummary[1:10]
  SMAcomparison[(11*(i-1))+1,12] = clades$commoncoef$p
  SMAcomparison[(11*(i-1))+1,13] = terrestrials$commoncoef$p
  SMAcomparison[(11*(i-1))+1,14] = BMcategories$commoncoef$p
  SMAcomparison[(11*(i-1))+1,15] = clades_e$gtr$p
  SMAcomparison[(11*(i-1))+1,16] = terrestrials_e$gtr$p
  SMAcomparison[(11*(i-1))+1,17] = BMcategories_e$gtr$p
  SMAcomparison[(11*(i-1))+2,2:11] = clades$groupsummary[1,1:10]
  SMAcomparison[(11*(i-1))+3,2:11] = clades$groupsummary[2,1:10]
  SMAcomparison[(11*(i-1))+4,2:11] = terrestrials$groupsummary[1,1:10]
  SMAcomparison[(11*(i-1))+5,2:11] = terrestrials$groupsummary[2,1:10]
  SMAcomparison[(11*(i-1))+6,2:11] = BMcategories$groupsummary[1,1:10]
  SMAcomparison[(11*(i-1))+7,2:11] = BMcategories$groupsummary[5,1:10]
  SMAcomparison[(11*(i-1))+8,2:11] = BMcategories$groupsummary[3,1:10]
  SMAcomparison[(11*(i-1))+9,2:11] = BMcategories$groupsummary[4,1:10]
  SMAcomparison[(11*(i-1))+10,2:11] = BMcategories$groupsummary[6,1:10]
  SMAcomparison[(11*(i-1))+11,2:11] = BMcategories$groupsummary[2,1:10]
 
  plot(sma(fullset$BM.g~fullset[,(9+i)]*fullset$Clade, data = fullset, method = "SMA", alpha = 0.05), 
       main=substitute(paste('log BMs vs ', b, " grouped by clade"), list(b=colnames(fullset)[9+i])), xlab = "log Surface Area", ylab = "log BM")+
    abline(result$groupsummary[,8], meanCSA_sma$groupsummary[,5])
  plot(sma(fullset$BM.g~fullset[,(9+i)]*fullset$Terrestriality, data = fullset, method = "SMA", alpha = 0.05),
       main=substitute(paste('log BMs vs ', b, " grouped by terrestriality"), list(b=colnames(fullset)[9+i])), xlab = "log Surface Area", ylab = "log BM")+
    abline(result$groupsummary[,8], meanCSA_sma$groupsummary[,5])
  plot(sma(fullset$BM.g~fullset[,(9+i)]*fullset$BM.Category, data = fullset, method = "SMA", alpha = 0.05),
       main=substitute(paste('log BMs vs ', b, " grouped by body mass"), list(b=colnames(fullset)[9+i])), xlab = "log Surface Area", ylab = "log BM")+
    abline(result$groupsummary[,8], meanCSA_sma$groupsummary[,5])
   
}

n = dev.off()

write.csv(SMAcomparison, "C:\\Wirc\\Projekt3\\Bigness\\Do DVCs even work\\testing\\SMA results summary.csv")


baseSMAcomparison[1,] = meanCSA_sma$groupsummary[2:10]
baseSMAcomparison[2,] = minCSA_sma$groupsummary[2:10]
baseSMAcomparison[3,] = maxCSA_sma$groupsummary[2:10]


meanCSA_sma_clades = sma(fullset$BM.g~fullset$Mean.CSA*fullset$Clade, data = fullset, method = "SMA", alpha = 0.05)
minCSA_sma_clades = sma(fullset$BM.g~fullset$Min.CSA*fullset$Clade, data = fullset, method = "SMA", alpha = 0.05)
maxCSA_sma_clades = sma(fullset$BM.g~fullset$Max.CSA*fullset$Clade, data = fullset, method = "SMA", alpha = 0.05)
