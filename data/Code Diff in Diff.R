# This script contains all necessary commands to reproduce the paper
# Greenstone and Hannah(2014), and extend it. The reproduction takes as input the 
# data cleaned by their do-file. (Dataset extracted just before their two step approach).
# We reproduce the two step approach with a built-in method and from scratch,
# as well as the second stage regression with their first stage stata output,
# and check that all these method yield the same output (lines 42-103)
# We propose a first extension with a different clustering and control choice
# (lines 109 -134)
# We propose a second extension with Sun and Abraham(2021) method (lines 148-218)
# To run this script, you need to have the "R-ready_do.dta", "R-ready_od.dta"
# and "R-ready_lnfcoli.dta" in your folder, and replace folder1 (line 24)
# by your actual folder name. Comparison with stata outputs requires files "sigmas_bod.dta"
# "sigmas_do.dta" and "sigmas_lnfcoli.dta" as well.
# Authors : Maël Laoufi and Yanis Zeghal - ENSAE/MIE

library(haven)        # To load .dta files
library(fixest)       # Built-in fixed effects regression with Sun & Abraham method
library(fastDummies)  # To split categorical variables into dummies when re-coding 
library(lmtest)       # To compute heteroskedastic robust Std errors. (Not used in final version)
library(sandwich)     # To compute any kind of asymptotic variance matrix (Not used un final version)
library(coefplot)     # To easily plot Sun & Abraham method outputs
library(texreg)       # To convert fitted models to LaTex arrays.

folder1 = "YourPath/" # Your folder name ending with a /
dir.create(file.path(folder1,"Outputs"), showWarnings = TRUE) # Creates the outputs folder
folder = paste(folder1,"Outputs/", sep = "")  # Name of the output folder
Output = list()                               # Will contain fitted models
pollutants = c("bod", "lnfcoli", "do")        # The three pollution measures used in the paper

for (pol in pollutants){
  loc = paste(folder1,"R-ready_",pol,".dta", sep = "")

print(paste("Pollutant : ", pol)) #Tells the advancement of the loop
A = read_dta(loc)
# Split the year and cityriver into categorical variables. We work with C
A$year = as.integer(A$year)
A$cityriver = as.integer(A$cityriver)
B = dummy_cols(A,"cityriver")
C = dummy_cols(B,"year")


# Automatically select the columns for year fe, cityriver fe, controls and Ds
indexyear = which(sapply(colnames(C), FUN = function(x) gregexpr("year_", x)[1]>0))
year = as.matrix(C[indexyear])
colnames(year)= gsub("year_","",colnames(year))
indexcityriver = which(sapply(colnames(C), FUN = function(x) gregexpr("cityriver_", x)[1]>0))
cityriver = as.matrix(C[indexcityriver])
colnames(cityriver)= gsub("cityriver_","",colnames(cityriver))
indexcontrols = which(sapply(colnames(C), FUN = function(x) x %in% c("lit_urban","pce")))
controls = as.matrix(C[indexcontrols])
indexD = which(sapply(colnames(C), FUN = function(x) gregexpr("tau",x)[1]>0))[-c(1,2)]
D = as.matrix(C[indexD])
colnames(D)= sub("tau","",colnames(D))
colnames(D)= sub("m","-",colnames(D))

# Reproduction of coefficients: 2 way fixed effects 
reprod_scratch = lm(C[[pol]]~D+controls+year+cityriver, weights = A$pop_urban)
# See which variables were dropped
names(which(is.na(reprod_scratch$coefficients)))
# Note that std errors are not clustered (yet), we only read coefficients 

# Reproduction of coefficients: with built-in method

# Adapting the formula to the name of pollutant
# /!\ Do not \n mid next line !
strform=paste(pol," ~ taum7 + taum6  + taum5 +taum4 +taum3 + taum2 + taum1 + tau0+ tau1 + tau2 + tau3 + tau4 + tau5 + tau6 + tau7 + tau8 + tau9 + tau10 + tauL + tauR+ lit_urban + pce | year + cityriver")
form = as.formula(strform) # This formula is an argument to feols below.

# Reproducing the 2-way fe OLS
reprod_builtin = feols(form, data = A, weights = A$pop_urban)

# The difference between built-in and from scratch coefficients should be a translation
# If so, d below should be 0. Here, computation errors sum up to e-5, which is huge
model_drift = as.matrix(reprod_scratch$coefficients[2:21] - reprod_builtin$coefficients[1:20])
d=norm(model_drift - mean(model_drift), type = 'I')
print(paste("Difference built-in and scratch-made : ",d ))

#Check that coefficients are indeed a translation of those obtained by Stata.
#The translation is due to the choice of which variables to drop (collinearity)
stata_pol = read_dta(paste(folder1, "sigmas_",pol,".dta", sep =""))
stata_coef = stata_pol$taub[1:18]

# Norm infinity distance between outputs to get a hint of computation discrepancy:
d2 = norm(as.matrix((reprod_scratch$coefficients[2:19]- stata_coef)-mean((reprod_scratch$coefficients[2:19]- stata_coef))),type = "I")
d3 = norm(as.matrix((reprod_builtin$coefficients[1:18]- stata_coef)-mean((reprod_builtin$coefficients[1:18]- stata_coef))),type = "I")
print(paste("Difference Stata and scratch-made : ",d2))
print(paste("Difference Stata and built-in     : ",d3))

# Paper 2nd stage reproduction output:
# Fed with R output
a2paper_reproduced=lm(reprod_scratch$coefficients[2:19]~stata_pol$nrcp, weights = 1/summary(reprod_scratch)$coefficients[2:19,2])
b2paper_reproduced=lm(reprod_scratch$coefficients[2:19]~stata_pol$nrcp + stata_pol$tau, weights = 1/summary(reprod_scratch)$coefficients[2:19,2])
c2paper_reproduced=lm(reprod_scratch$coefficients[2:19]~stata_pol$nrcp + stata_pol$tau + stata_pol$nrcp_trend, weights = 1/summary(reprod_scratch)$coefficients[2:19,2])
# Fed with Stata output
a2paper_real=lm(stata_pol$taub~stata_pol$nrcp, weights = 1/stata_pol$tause)
b2paper_real=lm(stata_pol$taub~stata_pol$nrcp + stata_pol$tau, weights = 1/stata_pol$tause)
c2paper_real=lm(stata_pol$taub~stata_pol$nrcp + stata_pol$tau + stata_pol$nrcp_trend, weights = 1/stata_pol$tause)

# Make it LaTex
texreg(list(a2paper_reproduced, b2paper_reproduced, c2paper_reproduced,a2paper_real,b2paper_real,c2paper_real),custom.model.names = paste("fit", 1:6),custom.coef.map = list("stata_pol$tau" = "time trend","stata_pol$nrcp"="1(policy)", "stata_pol$nrcp_trend"="1(policy)*trend"),
       omit.coef = "(cityriver)|(year)",stars = c(0.01, 0.05, 0.1),reorder.gof = c(4,3,1,2),custom.gof.rows = list("Equation fitted"=rep(c("2a","2b","2c"),2),"Coefficient and Se origin" = c(rep("R",3),rep("Stata",3))),
       include.adjrs = FALSE, caption = paste("Paper reproduction:  Estimates of the Effect of the NRCP on ", pol,".", sep = ""), caption.above =TRUE,
       custom.note = "%stars. Methodology",file = paste(folder,"LaTex_Paper_",pol,".txt",sep = ""))

# Extensions
###############################################################################

# Reproducing the 2-way fe OLS with different clustering
# Cluster on River
reprod_builtinRiver = feols(form, data = A, weights = A$pop_urban, cluster="river")
# Multi-Way-Clustering : River, City 
reprod_builtin2w = feols(form, data = A, weights = A$pop_urban, cluster = c("city","river"))

# Multi-Way-Clustering : River, City , no controls
strformNoCtrl=paste(pol," ~ taum7 + taum6  + taum5 +taum4 +taum3 + taum2 + taum1 + tau0+ tau1 + tau2 + tau3 + tau4 + tau5 + tau6 + tau7 + tau8 + tau9 + tau10 + tauL + tauR| year + cityriver")
formNoCtrl = as.formula(strformNoCtrl) # make the line above a formula for feols
reprod_builtin2wNoCtrl = feols(formNoCtrl, data = A, weights = A$pop_urban, cluster = c("city","river"))

# Second step regressions (note that standard errors are meaningless here, there is independance nowhere)
a2=lm(reprod_builtin2w$coefficients[1:18]~stata_pol$nrcp, weights = 1/reprod_builtin2w$se[1:18])
b2=lm(reprod_builtin2w$coefficients[1:18]~stata_pol$nrcp + stata_pol$tau, weights = 1/reprod_builtin2w$se[1:18])
c2=lm(reprod_builtin2w$coefficients[1:18]~stata_pol$nrcp + stata_pol$tau + stata_pol$nrcp_trend, weights = 1/reprod_builtin2w$se[1:18])

a2NC=lm(reprod_builtin2wNoCtrl$coefficients[1:18]~stata_pol$nrcp, weights = 1/reprod_builtin2wNoCtrl$se[1:18])
b2NC=lm(reprod_builtin2wNoCtrl$coefficients[1:18]~stata_pol$nrcp + stata_pol$tau, weights = 1/reprod_builtin2wNoCtrl$se[1:18])
c2NC=lm(reprod_builtin2wNoCtrl$coefficients[1:18]~stata_pol$nrcp + stata_pol$tau + stata_pol$nrcp_trend, weights = 1/reprod_builtin2wNoCtrl$se[1:18])


texreg(list(a2, b2, c2, a2NC, b2NC, c2NC),custom.coef.map = list("stata_pol$tau" = "time trend","stata_pol$nrcp"="1(policy)", "stata_pol$nrcp_trend"="1(policy)*trend"),
       omit.coef = "(cityriver)|(year)",stars = c(0.01, 0.05, 0.1),reorder.gof = c(5,4,3,2,1),custom.gof.rows = list("Equation fitted"=rep(c("2a","2b","2c"),2),"Multi-way clustering in first stage" = rep("yes",6),
       "Controls in first stage" = c(rep("No",3),rep("Yes",3))),
       include.adjrs = FALSE, caption = paste("Trend Break Estimates of the Effect of the NRCP on ", pol,".", sep = ""), caption.above =TRUE,
       custom.note = "%stars. COPY-PASTE note from paper please !",file = paste(folder,"LaTex_",pol,".txt",sep = ""))



# Observation : stata and R's non heteroskedastic robust asymptotic variances
# are not the same, but std errors remain (almost) proportionnal, so the
# weighting method in 2nd step is still workable : (uncomment line below)
#print(summary(reprod_scratch)$coefficients[2:19,2]/as.vector(stata_pol$tause))


# Clearly, the built-in method is a little bit "black-box".
# The from-scratch method is much more similar to Stata.

# We checked that none of the D s between -7 and 10 was omitted 

# Extension: Sun & Abraham method for heterogeneous treatment effects
# Consistency check: the policy is never adopted or for some unique D, each date.
all((apply(D, FUN = sum, MARGIN  = 1))+ A$neveradopt)

# To use the feols library efficiently, we need to find the cohorts from taus and years,
# and check consistency of the dataset
C$cohort = 0
for (group in sort(unique(C$cityriver))){
  Subset = C[C$cityriver == group,]
  #Check that either tau is missing everywhere or never missing for each group
  if ((!all(is.na(Subset$tau)))&(any(is.na(Subset$tau)))){
    print(c(i, "has partially NA"))
  }
  #Check that tau is missing iff the group is never treated
  if (!all(is.na(Subset$tau)==Subset$neveradopt)){
    print(c(i, " : missing taus and neveradopt do not coincide "))
  }
  Subset$cohort = Subset$year-Subset$tau
  if(all(!is.na(Subset$tau))){
    if (!(all(Subset$cohort==Subset$cohort[1]))){
      print(c(i, ' : cohort is not the same in the whole group'))
    }
  }
  C$cohort[C$cityriver == group] = Subset$cohort
}
# Never treated :
C$cohort[is.na(C$cohort)] = 1e7
# Consistency check  (we want 3 false in the lines below)
any(unique(C$cityriver[C$cohort == 1993])%in% unique(C$cityriver[C$cohort == 1995])) 
any(unique(C$cityriver[C$cohort == 1993])%in% unique(C$cityriver[C$cohort == 1e7]))
any(unique(C$cityriver[C$cohort == 1995])%in% unique(C$cityriver[C$cohort == 1e7]))

# Regression itself :
strformSA=paste(pol,"~lit_urban +pce + sunab(cohort, year)| cityriver + year")
formSA = as.formula(strformSA) # This formula is an argument to feols below.

SunAbr = feols(formSA, C, weights = C$pop_urban, cluster = c("city", "river"))
Output[[pol]] = SunAbr
# Estimation at -1 year is dropped (it is a reference). We rreinsert it 
# to fit the 2nd stage OLS.
coef = summary(SunAbr)$coeftable[5:21,1]
coef = c(coef[1:6],0,coef[7:17])
weights = summary(SunAbr)$coeftable[5:21,2]
weights = c(weights[1:6], mean(weights), weights[7:17])
# Second stage regression
a2SA=lm(coef~stata_pol$nrcp, weights = 1/weights)
b2SA=lm(coef~stata_pol$nrcp + stata_pol$tau, weights = 1/weights)
c2SA=lm(coef~stata_pol$nrcp + stata_pol$tau + stata_pol$nrcp_trend, weights = 1/weights)
# SunAbraham without controls :
strformSANC=paste(pol,"~sunab(cohort, year)| cityriver + year")
formSANC = as.formula(strformSANC) # This formula is an argument to feols below.

SunAbrNoCtrl = feols(formSANC, C, weights = C$pop_urban, cluster = c("city", "river"))
coef = summary(SunAbrNoCtrl)$coeftable[3:19,1]
coef = c(coef[1:6],0,coef[7:17])
weights = summary(SunAbr)$coeftable[3:19,2]
weights = c(weights[1:6], mean(weights), weights[7:17])

a2SANC=lm(coef~stata_pol$nrcp, weights = 1/weights)
b2SANC=lm(coef~stata_pol$nrcp + stata_pol$tau, weights = 1/weights)
c2SANC=lm(coef~stata_pol$nrcp + stata_pol$tau + stata_pol$nrcp_trend, weights = 1/weights)
# Export output to LaTex
texreg(list(a2SA, b2SA, c2SA, a2SANC, b2SANC, c2SANC),custom.coef.map = list("stata_pol$tau" = "time trend","stata_pol$nrcp"="1(policy)", "stata_pol$nrcp_trend"="1(policy)*trend"),
       omit.coef = "(cityriver)|(year)",stars = c(0.01, 0.05, 0.1),reorder.gof = c(5,4,3,2,1),custom.gof.rows = list("Equation fitted"=rep(c("2a","2b","2c"),2),"Multi-way clustering in first stage" = rep("yes",6),
                                                                                                                     "Controls in first stage" = c(rep("No",3),rep("Yes",3))),
       include.adjrs = FALSE, caption = paste("Extension - Trend Break Estimates of the Effect of the NRCP on ", pol,".", sep = ""), caption.above =TRUE,
       custom.note = "%stars. This table contains 2nd step regressions of our extension with Sun and Abraham method.",file = paste(folder,"SunAbraham_",pol,".txt",sep = ""))
#Plot estimated effects
png(filename=paste(folder,"iplot_",pol,".png", sep =""),width =500,height = 400 )
iplot(SunAbr, xlim=c(-7,10),main = "", xlab = "Year (lag)", ylab = "Estimated effect and 95% Conf.Int")
dev.off()
}


