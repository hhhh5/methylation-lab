## This one did not lead to the top hits being ...
meta = meta[c(
 # these are the first 15 smokers and non-smokers
 "GSM2260480","GSM2260481","GSM2260482","GSM2260483","GSM2260484"
,"GSM2260485","GSM2260486","GSM2260487","GSM2260488","GSM2260489"
,"GSM2260490","GSM2260491","GSM2260492","GSM2260493","GSM2260494"
,"GSM2260495","GSM2260496","GSM2260497","GSM2260498","GSM2260499"
,"GSM2260500","GSM2260501","GSM2260511","GSM2260514","GSM2260516"
,"GSM2260519","GSM2260525","GSM2260528","GSM2260530","GSM2260532"
,"GSM2260543" # same person as GSM2260485
,"GSM2260653" # this is the potentially contaminated sample
,"GSM1185585" # unrelated sample from another GSE, granulocytes instead of whole blood
,"GSM2219539" # unrelated sample of lung tissue
)]


# Load the necessary libraries

library(stringi)
library(magrittr)
library(data.table)
library(svd)
devtools::install_github("hhhh5/ewastools@v1.4") # not on CRAN
library(ewastools)

# Use the `?` operator to access the documentation of a function `fun`: `?fun`

# Importing the data -----------------------------
# 1. Read in the file `meta.csv` using `fread` from the data.table package, save it as object named `meta`.
# 2. Import the methylation data using the function `read_idats`.
#    You only have to provide the first part of the file name without `_Red.idat.gz` or `_Grn.idat.gz`.
#    Save as object named `meta`

meta = fread("meta.csv")
meth = read_idats(meta$id)

# Control metrics --------------------------------
# 1. Use the functions `control_metrics` and `sample_failure`
# 2. Are there any failed samples?

meth %>% control_metrics %>% sample_failure -> meta$failed

# Sex check --------------------------------------
# 1. Apply the function `check_sex` on the `meth` object to compute the normalized average total fluorescence intensity of probes targeting the X and Y chromosomes.
# 2. Use the function `predict_sex` to infer the sex of the sample donors from the methylation data.
# 3. Add column `predicted_sex` to `meta`
# 4. Are there samples where `sex != predicted_sex`? What is the sample ID?
# 5. Flag the problematic samples with a logical variable `exclude` in `meta`

meta[,c("X","Y"):=check_sex(meth)]

meta$predicted_sex = predict_sex(meta$X,meta$Y,which(meta$sex=="m"),which(meta$sex=="f"))

meta[sex!=predicted_sex]
plot(Y~X,data=meta,type="n")
text(Y~X,labels=id,col=ifelse(sex=="m",2,1),data=meta)

meta[,exclude:=FALSE]
meta[sex!=predicted_sex,exclude:=TRUE]

# Preprocessing ----------------------------------
# Includes detection p-values, dye bias correction and conversion to beta-values

meth %>% detectionP %>% mask(0.01) %>% correct_dye_bias %>% dont_normalize -> beta

# SNP outliers -----------------------------------

snps = which(rownames(beta) %like% "^rs")
gt = call_genotypes(beta[snps,])
meta$outlier = snp_outliers(gt)

stripchart(meta$outlier,method="jitter",pch=4)
abline(v=-4,lty="dotted",col=2)

meta[outlier>-4]
meta[outlier>-4,exclude:=TRUE]

# Principal component analysis -------------------
# 1. Get of subset of beta without probes on the X and Y chromosome

chrXY = meth$manifest$chr %in% c("X","Y") | meth$manifest$probe_type == "rs"
pcs = beta[-chrXY,]
pcs = pcs - rowMeans(pcs)
pcs = na.omit(pcs)
pcs = t(pcs)
pcs = trlan.svd(pcs,neig=2)
pcs = pcs$u

meta$pc1 = pcs[,1]
meta$pc2 = pcs[,2]

plot(pc1~pc2,col=ifelse(sex=="m",2,1),data=meta)
text(pc1~pc2,labels="  "%s+%meta$id,col=ifelse(sex=="m",2,1),data=meta,pos=4,offset=0)

meta[id=="SUBJ19",exclude:=TRUE]

# Leukocyte composition --------------------------
# 1. Estimate the leukocyte composition using the function `estimateLC`. What does the function return?
# 2. Add the cell proportion estimates as new columns to the `meta` data.table using the function `cbind`
# 3. Plot the proportions of granulocytes versus (the column named `GR`) using
# 4. Set the `exclude` flag TRUE for the conspicuous sample

meta = cbind(meta,estimateLC(beta))

plot(meta$GR,ylim=c(0,1))
text(meta$GR,labels=meta$id %s+% "  ",srt=90,pos=2,offset=0)

meta[id=="SUBJ4",exclude:=TRUE]








# EWAS -------------------------------------------

library(purrr)
library(furrr)

meta = meta[exclude==FALSE]

meta = meta[,.(
	 id
	,sex = factor(sex)
	,smoker = factor(smoker,levels=c("non-smoker","smoker"))
	,GR,MO,B,CD4,CD8,NK,nRBC
	)]

meta$id %>% read_idats %>% detectionP %>% mask(0.01) %>% correct_dye_bias %>% dont_normalize -> beta

plan("multiprocess",workers=8)

f = function(i){
	meta$cpg = beta[i,]
	m = lm(cpg~1+sex+smoker+GR+MO+B+CD4+CD8+NK+nRBC,data=meta)
	coef(summary(m))["smokersmoker",c(1,4)]
}

f = possibly(f,otherwise=c(NA,NA))

ewas_results = future_map(1:nrow(beta),f,.progress=TRUE)

ewas_results = do.call("rbind",ewas_results) %>% data.table
setnames(ewas_results,1:2,c("coef","pval"))
ewas_results$probe_id = rownames(beta)
ewas_results %<>% na.omit

ewas_results[,fdr:=p.adjust(pval,m="fdr")]
ewas_results2[fdr<0.05]

         coef         pval   probe_id        fdr
1: -0.2150846 4.157428e-08 cg21566642 0.01977623
2: -0.2847536 8.146008e-08 cg05575921 0.01977623



> ewas_results[fdr<0.05]
          coef         pval   probe_id          fdr
1: -0.18275357 1.978601e-07 cg03636183 1.371459e-02
2: -0.19395647 2.908239e-08 cg05951221 2.822167e-03
3: -0.29539155 6.653004e-11 cg21566642 3.228051e-05
4: -0.16953551 2.487618e-10 cg01940273 6.034986e-05
5: -0.07466515 1.278921e-07 cg16587010 1.034225e-02
6: -0.34860918 2.823757e-09 cg05575921 4.566976e-04
7: -0.12393670 1.501231e-08 cg21161138 1.821000e-03
