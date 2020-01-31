library(data.table)
library(purrr)
library(magrittr)
library(stringi)

meta = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30870&targ=gsm&form=text&view=brief"
meta = meta %>% map(readLines) %>% unlist
meta = split(meta,cumsum(meta %like% "^\\^SAMPLE = GSM"))
names(meta) = map(meta,1) %>% stri_match_first(regex="GSM\\d+")

imap(meta,function(s,acc){
s = strsplit(s,split=" = ",fixed=TRUE)
data.table(gsm=acc,variable=map_chr(s,1),value=map_chr(s,2))
}) -> meta

meta = rbindlist(meta)

### tidy the data
meta = meta[variable %chin% c("!Sample_characteristics_ch1")]
i = meta[variable == "!Sample_characteristics_ch1",which=TRUE]
ch = meta$value[i] %>% stri_split(fixed=": ")
meta$variable[i] = map_chr(ch,1)
meta$value   [i] = map_chr(ch,2)
rm(ch)

#  1: GSM765860      age 103 years
#  2: GSM765861      age   Newborn
#  3: GSM765862      age  97 years
#  4: GSM765863      age  95 years
#  5: GSM765864      age  97 years
#  6: GSM765865      age  97 years
#  7: GSM765866      age  98 years
#  8: GSM765867      age  96 years
#  9: GSM765868      age 100 years
# 10: GSM765869      age  90 years
# 11: GSM765870      age  91 years
# 12: GSM765871      age  92 years
# 13: GSM765872      age  89 years
# 14: GSM765873      age  90 years
# 15: GSM765874      age  89 years
# 16: GSM765875      age  90 years
# 17: GSM765876      age  91 years
# 18: GSM765877      age  89 years
# 19: GSM765878      age  89 years
# 20: GSM765879      age  90 years
# 21: GSM765880      age  90 years
# 22: GSM765881      age   Newborn
# 23: GSM765882      age   Newborn
# 24: GSM765883      age   Newborn
# 25: GSM765884      age   Newborn
# 26: GSM765885      age   Newborn
# 27: GSM765886      age   Newborn
# 28: GSM765887      age   Newborn
# 29: GSM765888      age   Newborn
# 30: GSM765889      age   Newborn
# 31: GSM765890      age   Newborn
# 32: GSM765891      age   Newborn
# 33: GSM765892      age   Newborn
# 34: GSM765893      age   Newborn
# 35: GSM765894      age   Newborn
# 36: GSM765895      age   Newborn
# 37: GSM765896      age   Newborn
# 38: GSM765897      age   Newborn
# 39: GSM765898      age   Newborn
# 40: GSM765899      age   Newborn

# GSM765878 89
# GSM765893 Newborn
# GSM765880 90

beta = fread("source.txt")
beta = beta[,.(ID_REF,GSM765878,GSM765893,GSM765880)]
names(beta) = c("ID_REF","A","B","C")
write.csv(beta,file="beta_matrix.csv",row.names=FALSE)

### Possible solution
library(data.table)
beta = read.csv("beta_matrix.csv",row.names=1)
markers = fread("13059_2013_3156_MOESM3_ESM.csv",skip=2)
markers = markers[,.(ID_REF=CpGmarker,coef=CoefficientTraining)]
intercept = markers[ID_REF=="(Intercept)"]$coef
markers = markers[ID_REF!="(Intercept)"]
i = match(markers$ID_REF,rownames(beta))
beta = beta[i,]
age = intercept+colSums(beta * markers$coef,na.rm=TRUE)
inverse_transformation <- function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }
inverse_transformation(age)

#          A          B          C 
# 87.7932361  0.2184761 88.3745105 