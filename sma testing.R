library(smatr)
fullset = read.csv("C:\\...\\Full set.csv")
fullset$Clade = as.factor(fullset$Clade)
fullset$Terrestriality= as.factor(fullset$Terrestriality)
fullset$BM.Category = as.factor(fullset$BM.Category)

meanOSA_sma = sma(fullset$BM.g~fullset$Mean.OSA, data = fullset, method = "SMA", alpha = 0.05)
minOSA_sma = sma(fullset$BM.g~fullset$Min.OSA, data = fullset, method = "SMA", alpha = 0.05)
maxOSA_sma = sma(fullset$BM.g~fullset$Max.OSA, data = fullset, method = "SMA", alpha = 0.05)


SMAcomparison = data.frame(matrix(nrow = 33, ncol = 17))
SMAcomparison[,1] = rep(c("meanOSA", "minOSA", "maxOSA"), each = 11)
colnames(SMAcomparison) = c("test variable","subset","n","r-squared", "p-value", "slope","slope CI lower", "slope CI upper", "intercept", "intercept CI lower", "intercept CI upper", "clade slopes p", "terrestriality slopes p", "BM category slopes p", "clade intercepts p", "terrestriality intercepts p", "BM category intercepts p")

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
 
}


meanOSA_sma_clades = sma(fullset$BM.g~fullset$Mean.OSA*fullset$Clade, data = fullset, method = "SMA", alpha = 0.05)
minOSA_sma_clades = sma(fullset$BM.g~fullset$Min.OSA*fullset$Clade, data = fullset, method = "SMA", alpha = 0.05)
maxOSA_sma_clades = sma(fullset$BM.g~fullset$Max.OSA*fullset$Clade, data = fullset, method = "SMA", alpha = 0.05)

meanOSA_sma_BMcats = sma(fullset$BM.g~fullset$Mean.OSA*fullset$BM.Category, data = fullset, method = "SMA", alpha = 0.05)
minOSA_sma_BMcats = sma(fullset$BM.g~fullset$Min.OSA*fullset$BM.Category, data = fullset, method = "SMA", alpha = 0.05)
maxOSA_sma_BMcats = sma(fullset$BM.g~fullset$Max.OSA*fullset$BM.Category, data = fullset, method = "SMA", alpha = 0.05)

meanOSA_sma_terr = sma(fullset$BM.g~fullset$Mean.OSA*fullset$Terrestriality, data = fullset, method = "SMA", alpha = 0.05)
minOSA_sma_terr = sma(fullset$BM.g~fullset$Min.OSA*fullset$Terrestriality, data = fullset, method = "SMA", alpha = 0.05)
maxOSA_sma_terr = sma(fullset$BM.g~fullset$Max.OSA*fullset$Terrestriality, data = fullset, method = "SMA", alpha = 0.05)
