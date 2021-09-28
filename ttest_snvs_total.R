###################################################################################
##paired t-test for total SNVs
library(dplyr)
setwd("D:/maulana/GSE_second_paper/september_version/raw_images")
data = read.table("SNVs_total.txt", col.names = c("breed", "group", "ars", "uoa"))
#taurine
taurine = filter(data, group == "Taurine")
# Shapiro-Wilk normality test for the differences
#shapiro.test(taurine$ars - taurine$uoa)
t.test(taurine$ars, taurine$uoa, paired = TRUE)

#indicine
indicine = filter(data, group == "Indicine")
t.test(indicine$ars, indicine$uoa, paired = TRUE)

#when Bohai is assigned as taurine
data[7,2] = as.factor("Taurine")
#taurine
taurine = filter(data, group == "Taurine")
t.test(taurine$ars, taurine$uoa, paired = TRUE)

#indicine
indicine = filter(data, group == "Indicine")
t.test(indicine$ars, indicine$uoa, paired = TRUE)

##################################################################################
##paired t-test for total mapped and retained reads - group dataset (from table 1 in manuscript)
setwd("D:/maulana/GSE_second_paper/september_version/raw_images")
list.files()
data = read.table("total_retain_groups.txt", header = TRUE)
t.test(data$Map.ars, data$Map.uoa, paired = TRUE, alternative = "greater") #mapped reads
t.test(data$Ret.ars, data$Ret.uoa, paired = TRUE, alternative = "greater") #retained reads

##################################################################################
##ggplot for admixture CV plot
library(dplyr)
library(ggplot2)
#read data
setwd("D:/maulana/GSE_second_paper/september_version/raw_images")
list.files()
data = read.table("CV_error_admixture.txt", header = TRUE)
# plot
ggplot(data = data, aes(x=K, y=CV.error)) + geom_line(aes(colour=SNVs))
