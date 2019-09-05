setwd("~/Box Sync/amd338@cornell.edu - Annise Dobson's Files and Folders/Root architecture")
library(lavaan)
library(qgraph)
library(semPlot)
library(lme4)
library(mice)
library(devtools)
#install_github("jslefche/piecewiseSEM")
library(piecewiseSEM)
# versions::available.versions("piecewiseSEM")
versions::install.versions("piecewiseSEM", "1.2.1")
library(ape)
library(caper)
library(nlme)
library(glmm)
library(pscl)
library(dplyr)
library(AICcmodavg)

x=read.csv("Dobson.Data.csv",header=T) %>%
  mutate(Deer = as.numeric(Deer))%>%
  mutate(Plot = as.factor(Plot))

##Look at missing data
#aggr_plot <- aggr(x, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))
#tempData <- mice(x,m=5,maxit=50,meth='pmm',seed=500)
#summary(tempData)
#imputed <- complete(tempData,1)
#write.csv(imputed, file = "ImputedData.csv",row.names=FALSE, na="",col.names=T, sep=",")

############Survival/Concentration models
########
byspecies <- x[ which(x$Species=='Actaea'), ]

m=glmer(Survival.Binomial ~Dry.Year.gm2 + Deer + (1 | Site), data = byspecies,
        family = binomial(link = "logit"), nAGQ = 25)
ma=glm(Survival.Binomial ~Dry.Year.gm2 + Deer, data = byspecies,
        family = binomial(link = "logit"))
m1 <- glmer(Survival.Binomial ~Dry.Year.gm2 + Deer + SoilN + P.Total.Conc + P.Exch.Conc + (1 | Site), data = byspecies,
             family = binomial(link = "logit"), nAGQ = 25)
m1a <- glm(Survival.Binomial ~Dry.Year.gm2 + Deer + SoilN + P.Total.Conc + P.Exch.Conc, data = byspecies,
            family = binomial(link = "logit"))
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m3=lme(L.4.L~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m4=lme(SoilN~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m5=lme(P.Total.Conc~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m6=lme(P.Exch.Conc~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)

mA = list(
  m1a,
  m2,
  m3,
  m4,
  m5,
  m6
)

# Explore individual model fits
sem.model.fits(mA)

# Run goodness-of-fit tests
sem.fit(mA, byspecies,conditional=T,
        corr.errors=c("P.Exch.Conc~~P.Total.Conc",
                      "P.Exch.Conc~~All.Mycorhizae",
                      "P.Exch.Conc~~L.4.L",
                      "P.Total.Conc~~All.Mycorhizae",
                      "P.Total.Conc~~L.4.L",
                      "All.Mycorhizae~~L.4.L",
                      "SoilN~~P.Exch.Conc",
                      "SoilN~~P.Total.Conc",
                      "SoilN~~All.Mycorhizae",
                      "SoilN~~L.4.L"
        ))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("P.Exch.Conc~~P.Total.Conc",
                        "P.Exch.Conc~~All.Mycorhizae",
                        "P.Exch.Conc~~L.4.L",
                        "P.Total.Conc~~All.Mycorhizae",
                        "P.Total.Conc~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "SoilN~~P.Exch.Conc",
                        "SoilN~~P.Total.Conc",
                        "SoilN~~All.Mycorhizae",
                        "SoilN~~L.4.L"
          ))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("P.Exch.Conc~~P.Total.Conc",
                        "P.Exch.Conc~~All.Mycorhizae",
                        "P.Exch.Conc~~L.4.L",
                        "P.Total.Conc~~All.Mycorhizae",
                        "P.Total.Conc~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "SoilN~~P.Exch.Conc",
                        "SoilN~~P.Total.Conc",
                        "SoilN~~All.Mycorhizae",
                        "SoilN~~L.4.L"
          ))



coef.table2=sem.coefs(mA,data=byspecies,standardize="scale",
                      corr.errors=c("P.Exch.Conc~~P.Total.Conc",
                                    "P.Exch.Conc~~All.Mycorhizae",
                                    "P.Exch.Conc~~L.4.L",
                                    "P.Total.Conc~~All.Mycorhizae",
                                    "P.Total.Conc~~L.4.L",
                                    "All.Mycorhizae~~L.4.L",
                                    "SoilN~~P.Exch.Conc",
                                    "SoilN~~P.Total.Conc",
                                    "SoilN~~All.Mycorhizae",
                                    "SoilN~~L.4.L"
                      ))
summary(ma)
summary(m2)
summary(m3)
summary(m4)
summary(m5)
summary(m6)

sem.plot(mA, byspecies,coef.table2, corr.errors = NULL,
         show.nonsig = TRUE, scaling = 15, alpha = 0.05)

##Aquilegia####################################################
byspecies <- x[ which(x$Species=='Aquilegia'), ]

m=glmer(Survival.Binomial ~Dry.Year.gm2 + Deer + (1 | Site), data = byspecies,
        family = binomial(link = "logit"), nAGQ = 25)
ma=glm(Survival.Binomial ~Dry.Year.gm2 + Deer, data = byspecies, family = binomial(link = "logit"))
m1 <- glmer(Survival.Binomial ~Dry.Year.gm2 + Deer + SoilN + P.Total.Conc + P.Exch.Conc + (1 | Site), data = byspecies,
            family = binomial(link = "logit"), nAGQ = 25)
m1a <- glm(Survival.Binomial ~Dry.Year.gm2 + Deer + SoilN + P.Total.Conc + P.Exch.Conc, data = byspecies,
           family = binomial(link = "logit"))
m1b <- glmer(Survival.Binomial ~Dry.Year.gm2 + Deer + SoilN + P.Total.Conc + P.Exch.Conc + (1 | Site), data = byspecies,
             family = binomial(link = "logit"), nAGQ = 25)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m3=lme(L.4.L~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m4=lme(SoilN~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m5=lme(P.Total.Conc~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m6=lme(P.Exch.Conc~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)

mA = list(
  m1a,
  m2,
  m3,
  m4,
  m5,
  m6
)

# Explore individual model fits
sem.model.fits(mA)

# Run goodness-of-fit tests
sem.fit(mA, byspecies,conditional=T,
        corr.errors=c("P.Exch.Conc~~P.Total.Conc",
                      "P.Exch.Conc~~All.Mycorhizae",
                      "P.Exch.Conc~~L.4.L",
                      "P.Total.Conc~~All.Mycorhizae",
                      "P.Total.Conc~~L.4.L",
                      "All.Mycorhizae~~L.4.L",
                      "SoilN~~P.Exch.Conc",
                      "SoilN~~P.Total.Conc",
                      "SoilN~~All.Mycorhizae",
                      "SoilN~~L.4.L"))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("P.Exch.Conc~~P.Total.Conc",
                        "P.Exch.Conc~~All.Mycorhizae",
                        "P.Exch.Conc~~L.4.L",
                        "P.Total.Conc~~All.Mycorhizae",
                        "P.Total.Conc~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "SoilN~~P.Exch.Conc",
                        "SoilN~~P.Total.Conc",
                        "SoilN~~All.Mycorhizae",
                        "SoilN~~L.4.L"))

# Explore individual model fits
sem.model.fits(mA)

coef.table2=sem.coefs(mA,data=byspecies,standardize="scale",
                      corr.errors=c("P.Exch.Conc~~P.Total.Conc",
                                    "P.Exch.Conc~~All.Mycorhizae",
                                    "P.Exch.Conc~~L.4.L",
                                    "P.Total.Conc~~All.Mycorhizae",
                                    "P.Total.Conc~~L.4.L",
                                    "All.Mycorhizae~~L.4.L",
                                    "SoilN~~P.Exch.Conc",
                                    "SoilN~~P.Total.Conc",
                                    "SoilN~~All.Mycorhizae",
                                    "SoilN~~L.4.L"))
sem.coefs(mA,data=byspecies)
summary(m1a)
summary(m2)
summary(m3)
summary(m4)
summary(m5)
summary(m6)

sem.plot(mA, byspecies,coef.table2, corr.errors = NULL,
         show.nonsig = TRUE, scaling = 15, alpha = 0.05)
View(coef.table2)
##Cornus####################################################
byspecies <- x[ which(x$Species=='Cornus'), ]

m=glmer(Survival.Binomial ~Dry.Year.gm2 + as.numeric(Deer) + (1 | Site), data = byspecies,
        family = binomial(link = "logit"), nAGQ = 100)
ma=glm(Survival.Binomial ~Dry.Year.gm2 + Deer, data = byspecies,
        family = binomial(link = "logit"))
m1 <- glmer(Survival.Binomial ~Dry.Year.gm2 + as.numeric(Deer) ++ SoilN + P.Total.Conc + P.Exch.Conc + (1 | Site), data = byspecies,
            family = binomial(link = "logit"), nAGQ = 100)
m1a <- glm(Survival.Binomial ~Dry.Year.gm2 + Deer+ SoilN + P.Total.Conc + P.Exch.Conc, data = byspecies,
            family = binomial(link = "logit"))
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m3=lme(L.4.L~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m4=lme(SoilN~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m5=lme(P.Total.Conc~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m6=lme(P.Exch.Conc~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)

mA = list(
  m1a,
  m2,
  m3,
  m4,
  m5,
  m6
)

# Run goodness-of-fit tests
sem.fit(mA, byspecies,conditional=T,
        corr.errors=c("P.Exch.Conc~~P.Total.Conc",
                      "P.Exch.Conc~~All.Mycorhizae",
                      "P.Exch.Conc~~L.4.L",
                      "P.Total.Conc~~All.Mycorhizae",
                      "P.Total.Conc~~L.4.L",
                      "All.Mycorhizae~~L.4.L",
                      "SoilN~~P.Exch.Conc",
                      "SoilN~~P.Total.Conc",
                      "SoilN~~All.Mycorhizae",
                      "SoilN~~L.4.L"))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("P.Exch.Conc~~P.Total.Conc",
                        "P.Exch.Conc~~All.Mycorhizae",
                        "P.Exch.Conc~~L.4.L",
                        "P.Total.Conc~~All.Mycorhizae",
                        "P.Total.Conc~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "SoilN~~P.Exch.Conc",
                        "SoilN~~P.Total.Conc",
                        "SoilN~~All.Mycorhizae",
                        "SoilN~~L.4.L"))

# Explore individual model fits
sem.model.fits(mA)

coef.table2=sem.coefs(mA,data=byspecies,standardize="scale")

sem.coefs(mA,data=byspecies)

summary(m1)
summary(m2)
summary(m3)
summary(m4)
summary(m5)
summary(m6)

sem.plot(mA, byspecies,coef.table2, corr.errors = NULL,
         show.nonsig = TRUE, scaling = 15, alpha = 0.05)

View(coef.table2)
##Prenanthes####################################################
byspecies <- x[ which(x$Species=='Prenanthes'), ]

m=glmer(Survival.Binomial ~Dry.Year.gm2 + Deer + (1 | Site), data = byspecies,
        family = binomial(link = "logit"), nAGQ = 100)
ma=glm(Survival.Binomial ~Dry.Year.gm2 + Deer, data = byspecies,
       family = binomial(link = "logit"))
m1 <- glmer(Survival.Binomial ~Dry.Year.gm2 + Deer + SoilN + P.Total.Conc + P.Exch.Conc + (1 | Site), data = byspecies,
            family = binomial(link = "logit"), nAGQ = 100)
m1a <- glm(Survival.Binomial ~Dry.Year.gm2 + Deer + SoilN + P.Total.Conc + P.Exch.Conc, data = byspecies,
           family = binomial(link = "logit"))
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m3=lme(L.4.L~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m4=lme(SoilN~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m5=lme(P.Total.Conc~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m6=lme(P.Exch.Conc~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)

mA = list(
  ma,
  m2,
  m3,
  m4,
  m5,
  m6
)

# Run goodness-of-fit tests
sem.fit(mA, byspecies,conditional=T,
        corr.errors=c("P.Exch.Conc~~P.Total.Conc",
                      "P.Exch.Conc~~All.Mycorhizae",
                      "P.Exch.Conc~~L.4.L",
                      "P.Total.Conc~~All.Mycorhizae",
                      "P.Total.Conc~~L.4.L",
                      "All.Mycorhizae~~L.4.L",
                      "SoilN~~P.Exch.Conc",
                      "SoilN~~P.Total.Conc",
                      "SoilN~~All.Mycorhizae",
                      "SoilN~~L.4.L"))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("P.Exch.Conc~~P.Total.Conc",
                        "P.Exch.Conc~~All.Mycorhizae",
                        "P.Exch.Conc~~L.4.L",
                        "P.Total.Conc~~All.Mycorhizae",
                        "P.Total.Conc~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "SoilN~~P.Exch.Conc",
                        "SoilN~~P.Total.Conc",
                        "SoilN~~All.Mycorhizae",
                        "SoilN~~L.4.L"))

# Explore individual model fits
sem.model.fits(mA)

coef.table2=sem.coefs(mA,data=byspecies,standardize="scale")
sem.coefs(mA,data=byspecies)

summary(ma)
summary(m2)
summary(m3)
summary(m4)
summary(m5)
summary(m6)

sem.plot(mA, byspecies,coef.table2, corr.errors = NULL,
         show.nonsig = TRUE, scaling = 15, alpha = 0.05)
View(coef.table2)
##Quercus####################################################
byspecies <- x[ which(x$Species=='Quercus'), ]

m=glmer(Survival.Binomial ~Dry.Year.gm2 + Deer + (1 | Site), data = byspecies,
        family = binomial(link = "logit"), nAGQ = 100)
ma=glm(Survival.Binomial ~Dry.Year.gm2 + Deer, data = byspecies,
       family = binomial(link = "logit"))
m1 <- glmer(Survival.Binomial ~Dry.Year.gm2 + Deer + SoilN + P.Total.Conc + P.Exch.Conc + (1 | Site), data = byspecies,
            family = binomial(link = "logit"), nAGQ = 100)
m1a <- glm(Survival.Binomial ~Dry.Year.gm2 + Deer + SoilN + P.Total.Conc + P.Exch.Conc, data = byspecies,
           family = binomial(link = "logit"))
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m3=lme(L.4.L~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m4=lme(SoilN~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m5=lme(P.Total.Conc~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m6=lme(P.Exch.Conc~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)

mA = list(
  m1a,
  m2,
  m3,
  m4,
  m5,
  m6
)

# Run goodness-of-fit tests
sem.fit(mA, byspecies,conditional=T,
        corr.errors=c("P.Exch.Conc~~P.Total.Conc",
                      "P.Exch.Conc~~All.Mycorhizae",
                      "P.Exch.Conc~~L.4.L",
                      "P.Total.Conc~~All.Mycorhizae",
                      "P.Total.Conc~~L.4.L",
                      "All.Mycorhizae~~L.4.L",
                      "SoilN~~P.Exch.Conc",
                      "SoilN~~P.Total.Conc",
                      "SoilN~~All.Mycorhizae",
                      "SoilN~~L.4.L"))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("P.Exch.Conc~~P.Total.Conc",
                        "P.Exch.Conc~~All.Mycorhizae",
                        "P.Exch.Conc~~L.4.L",
                        "P.Total.Conc~~All.Mycorhizae",
                        "P.Total.Conc~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "SoilN~~P.Exch.Conc",
                        "SoilN~~P.Total.Conc",
                        "SoilN~~All.Mycorhizae",
                        "SoilN~~L.4.L"))

# Explore individual model fits
sem.model.fits(mA)

coef.table2=sem.coefs(mA,data=byspecies,standardize="scale")
sem.coefs(mA,data=byspecies)

summary(m1a)
summary(m2)
summary(m3)
summary(m4)
summary(m5)
summary(m6)

sem.plot(mA, byspecies,coef.table2, corr.errors = NULL,
         show.nonsig = TRUE, scaling = 15, alpha = 0.05)
View(coef.table2)
#########Survival/Pools
########
byspecies <- x[ which(x$Species=='Actaea'), ]

m=glmer(Survival.Binomial~Dry.Year.gm2 + as.numeric(Deer) +1|Site, data = byspecies,na.action=na.omit)
m1=glmer(Survival.Binomial~Dry.Year.gm2 + as.numeric(Deer) + SoilN.20Pool + Standard.exchange.P.20 + Standard.Total.P.20 +1|Site, data = byspecies,na.action=na.omit)
m1a=glm(Survival.Binomial~Dry.Year.gm2 + as.numeric(Deer) + SoilN.20Pool + Standard.exchange.P.20 + Standard.Total.P.20, data = byspecies,na.action=na.omit)
m1b=glm(Survival.Binomial~Dry.Year.gm2 + as.numeric(Deer), data = byspecies,na.action=na.omit)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)

m3=lme(L.4.L~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m4=lme(SoilN.20Pool~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m5=lme(Standard.Total.P.20~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m6=lme(Standard.exchange.P.20~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)

mA = list(
  m1b,
  m2,
  m3,
  m4,
  m5,
  m6
)

# Run goodness-of-fit tests
sem.fit(mA, byspecies,conditional=T,
        corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                      "Standard.exchange.P.20~~All.Mycorhizae",
                      "Standard.exchange.P.20~~L.4.L",
                      "Standard.Total.P.20~~All.Mycorhizae",
                      "Standard.Total.P.20~~L.4.L",
                      "All.Mycorhizae~~L.4.L",
                      "SoilN.20Pool~~Standard.exchange.P.20",
                      "SoilN.20Pool~~Standard.Total.P.20",
                      "SoilN.20Pool~~All.Mycorhizae",
                      "SoilN.20Pool~~L.4.L"
        ))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                        "Standard.exchange.P.20~~All.Mycorhizae",
                        "Standard.exchange.P.20~~L.4.L",
                        "Standard.Total.P.20~~All.Mycorhizae",
                        "Standard.Total.P.20~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "SoilN.20Pool~~Standard.exchange.P.20",
                        "SoilN.20Pool~~Standard.Total.P.20",
                        "SoilN.20Pool~~All.Mycorhizae",
                        "SoilN.20Pool~~L.4.L"
          ))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                        "Standard.exchange.P.20~~All.Mycorhizae",
                        "Standard.exchange.P.20~~L.4.L",
                        "Standard.Total.P.20~~All.Mycorhizae",
                        "Standard.Total.P.20~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "SoilN.20Pool~~Standard.exchange.P.20",
                        "SoilN.20Pool~~Standard.Total.P.20",
                        "SoilN.20Pool~~All.Mycorhizae",
                        "SoilN.20Pool~~L.4.L"
          ))



coef.table2=sem.coefs(mA,data=byspecies,standardize="scale",
                      corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                                    "Standard.exchange.P.20~~All.Mycorhizae",
                                    "Standard.exchange.P.20~~L.4.L",
                                    "Standard.Total.P.20~~All.Mycorhizae",
                                    "Standard.Total.P.20~~L.4.L",
                                    "All.Mycorhizae~~L.4.L",
                                    "SoilN.20Pool~~Standard.exchange.P.20",
                                    "SoilN.20Pool~~Standard.Total.P.20",
                                    "SoilN.20Pool~~All.Mycorhizae",
                                    "SoilN.20Pool~~L.4.L"
                      ))
summary(m1b)
summary(m2)
summary(m3)
summary(m4)
summary(m5)
summary(m6)

sem.plot(mA, byspecies,coef.table2, corr.errors = NULL,
         show.nonsig = TRUE, scaling = 15, alpha = 0.05)
sem.model.fits(mA)
View(coef.table2)
##Aquilegia####################################################
byspecies <- x[ which(x$Species=='Aquilegia'), ]

m=glmer(Survival.Binomial~Dry.Year.gm2 + as.numeric(Deer) +1|Site, data = byspecies,na.action=na.omit)
m1=glmer(Survival.Binomial~Dry.Year.gm2 + as.numeric(Deer) + SoilN.20Pool + Standard.exchange.P.20 + Height + Standard.Total.P.20 +1|Site, data = byspecies,na.action=na.omit)
m1a=glm(Survival.Binomial~Dry.Year.gm2 + as.numeric(Deer)+ SoilN.20Pool + Standard.exchange.P.20 + Height + Standard.Total.P.20, data = byspecies,na.action=na.omit)
m1b=glm(Survival.Binomial~Dry.Year.gm2 + as.numeric(Deer), data = byspecies,na.action=na.omit)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m3=lme(L.4.L~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m4=lme(SoilN.20Pool~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m5=lme(Standard.Total.P.20~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m6=lme(Standard.exchange.P.20~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)

mA = list(
  m1b,
  m2,
  m3,
  m4,
  m5,
  m6
)

# Run goodness-of-fit tests
sem.fit(mA, byspecies,conditional=T,
        corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                      "Standard.exchange.P.20~~All.Mycorhizae",
                      "Standard.exchange.P.20~~L.4.L",
                      "Standard.Total.P.20~~All.Mycorhizae",
                      "Standard.Total.P.20~~L.4.L",
                      "All.Mycorhizae~~L.4.L",
                      "SoilN.20Pool~~Standard.exchange.P.20",
                      "SoilN.20Pool~~Standard.Total.P.20",
                      "SoilN.20Pool~~All.Mycorhizae",
                      "SoilN.20Pool~~L.4.L"))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                        "Standard.exchange.P.20~~All.Mycorhizae",
                        "Standard.exchange.P.20~~L.4.L",
                        "Standard.Total.P.20~~All.Mycorhizae",
                        "Standard.Total.P.20~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "SoilN.20Pool~~Standard.exchange.P.20",
                        "SoilN.20Pool~~Standard.Total.P.20",
                        "SoilN.20Pool~~All.Mycorhizae",
                        "SoilN.20Pool~~L.4.L"))

# Explore individual model fits
sem.model.fits(mA)

coef.table2=sem.coefs(mA,data=byspecies,standardize="scale",
                      corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                                    "Standard.exchange.P.20~~All.Mycorhizae",
                                    "Standard.exchange.P.20~~L.4.L",
                                    "Standard.Total.P.20~~All.Mycorhizae",
                                    "Standard.Total.P.20~~L.4.L",
                                    "All.Mycorhizae~~L.4.L",
                                    "SoilN.20Pool~~Standard.exchange.P.20",
                                    "SoilN.20Pool~~Standard.Total.P.20",
                                    "SoilN.20Pool~~All.Mycorhizae",
                                    "SoilN.20Pool~~L.4.L"))
sem.coefs(mA,data=byspecies)
summary(m1a)
summary(m2)
summary(m3)
summary(m4)
summary(m5)
summary(m6)

sem.plot(mA, byspecies,coef.table2, corr.errors = NULL,
         show.nonsig = TRUE, scaling = 15, alpha = 0.05)
View(coef.table2)
##Cornus####################################################
byspecies <- x[ which(x$Species=='Cornus'), ]

m=glmer(Survival.Binomial~Dry.Year.gm2+as.numeric(Deer) +1|Site, data = byspecies,na.action=na.omit)
m1=glmer(Survival.Binomial~Dry.Year.gm2+as.numeric(Deer) + SoilN.20Pool + Standard.exchange.P.20 + Height + Standard.Total.P.20 +1|Site, data = byspecies,na.action=na.omit)
m1a=glm(Survival.Binomial~Dry.Year.gm2 + as.numeric(Deer)+ SoilN.20Pool + Standard.exchange.P.20 + Height + Standard.Total.P.20, data = byspecies,na.action=na.omit)
m1b=glm(Survival.Binomial~Dry.Year.gm2 + as.numeric(Deer), data = byspecies,na.action=na.omit)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m3=lme(L.4.L~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m4=lme(SoilN.20Pool~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m5=lme(Standard.Total.P.20~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m6=lme(Standard.exchange.P.20~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)

mA = list(
  m1a,
  m2,
  m3,
  m4,
  m5,
  m6
)

# Run goodness-of-fit tests
sem.fit(mA, byspecies,conditional=T,
        corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                      "Standard.exchange.P.20~~All.Mycorhizae",
                      "Standard.exchange.P.20~~L.4.L",
                      "Standard.Total.P.20~~All.Mycorhizae",
                      "Standard.Total.P.20~~L.4.L",
                      "All.Mycorhizae~~L.4.L",
                      "SoilN.20Pool~~Standard.exchange.P.20",
                      "SoilN.20Pool~~Standard.Total.P.20",
                      "SoilN.20Pool~~All.Mycorhizae",
                      "SoilN.20Pool~~L.4.L"))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                        "Standard.exchange.P.20~~All.Mycorhizae",
                        "Standard.exchange.P.20~~L.4.L",
                        "Standard.Total.P.20~~All.Mycorhizae",
                        "Standard.Total.P.20~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "SoilN.20Pool~~Standard.exchange.P.20",
                        "SoilN.20Pool~~Standard.Total.P.20",
                        "SoilN.20Pool~~All.Mycorhizae",
                        "SoilN.20Pool~~L.4.L"))

# Explore individual model fits
sem.model.fits(mA)

coef.table2=sem.coefs(mA,data=byspecies,standardize="scale")
sem.coefs(mA,data=byspecies)

summary(m1)
summary(m2)
summary(m3)
summary(m4)
summary(m5)
summary(m6)

sem.plot(mA, byspecies,coef.table2, corr.errors = NULL,
         show.nonsig = TRUE, scaling = 15, alpha = 0.05)
View(coef.table2)
##Prenanthes####################################################
byspecies <- x[ which(x$Species=='Prenanthes'), ]

m=glmer(Survival.Binomial~Dry.Year.gm2 + as.numeric(Deer)+ Height+1|Site, data = byspecies,na.action=na.omit)
m1=glmer(Survival.Binomial~Dry.Year.gm2 + as.numeric(Deer) + SoilN.20Pool + Standard.exchange.P.20 + Height + Standard.Total.P.20 +1|Site, data = byspecies,na.action=na.omit)
m1a=glm(Survival.Binomial~Dry.Year.gm2 + as.numeric(Deer)+ SoilN.20Pool + Standard.exchange.P.20 + Height + Standard.Total.P.20, data = byspecies,na.action=na.omit)
m1b=glm(Survival.Binomial~Dry.Year.gm2 + as.numeric(Deer), data = byspecies,na.action=na.omit)

m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m3=lme(L.4.L~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m4=lme(SoilN.20Pool~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m5=lme(Standard.Total.P.20~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m6=lme(Standard.exchange.P.20~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)

mA = list(
  m1b,
  m2,
  m3,
  m4,
  m5,
  m6
)

# Run goodness-of-fit tests
sem.fit(mA, byspecies,conditional=T,
        corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                      "Standard.exchange.P.20~~All.Mycorhizae",
                      "Standard.exchange.P.20~~L.4.L",
                      "Standard.Total.P.20~~All.Mycorhizae",
                      "Standard.Total.P.20~~L.4.L",
                      "All.Mycorhizae~~L.4.L",
                      "SoilN.20Pool~~Standard.exchange.P.20",
                      "SoilN.20Pool~~Standard.Total.P.20",
                      "SoilN.20Pool~~All.Mycorhizae",
                      "SoilN.20Pool~~L.4.L"))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                        "Standard.exchange.P.20~~All.Mycorhizae",
                        "Standard.exchange.P.20~~L.4.L",
                        "Standard.Total.P.20~~All.Mycorhizae",
                        "Standard.Total.P.20~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "SoilN.20Pool~~Standard.exchange.P.20",
                        "SoilN.20Pool~~Standard.Total.P.20",
                        "SoilN.20Pool~~All.Mycorhizae",
                        "SoilN.20Pool~~L.4.L"))

# Explore individual model fits
sem.model.fits(mA)

coef.table2=sem.coefs(mA,data=byspecies,standardize="scale")
sem.coefs(mA,data=byspecies)

summary(m1)
summary(m2)
summary(m3)
summary(m4)
summary(m5)
summary(m6)

sem.plot(mA, byspecies,coef.table2, corr.errors = NULL,
         show.nonsig = TRUE, scaling = 15, alpha = 0.05)
View(coef.table2)
##Quercus####################################################
byspecies <- x[ which(x$Species=='Quercus'), ]

m=glmer(Survival.Binomial~Dry.Year.gm2 + as.numeric(Deer)  +1|Site, data = byspecies,na.action=na.omit)
m1=glmer(Survival.Binomial~Dry.Year.gm2 + as.numeric(Deer) + SoilN.20Pool + Standard.exchange.P.20 + Height + Standard.Total.P.20 +1|Site, data = byspecies,na.action=na.omit)
m1a=glm(Survival.Binomial~Dry.Year.gm2 + as.numeric(Deer)+ + SoilN.20Pool + Standard.exchange.P.20 + Height + Standard.Total.P.20, data = byspecies,na.action=na.omit)
m1b=glm(Survival.Binomial~Dry.Year.gm2 + as.numeric(Deer), data = byspecies,na.action=na.omit)

m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m3=lme(L.4.L~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m4=lme(SoilN.20Pool~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m5=lme(Standard.Total.P.20~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m6=lme(Standard.exchange.P.20~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)

mA = list(
  m1b,
  m2,
  m3,
  m4,
  m5,
  m6
)

# Run goodness-of-fit tests
sem.fit(mA, byspecies,conditional=T,
        corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                      "Standard.exchange.P.20~~All.Mycorhizae",
                      "Standard.exchange.P.20~~L.4.L",
                      "Standard.Total.P.20~~All.Mycorhizae",
                      "Standard.Total.P.20~~L.4.L",
                      "All.Mycorhizae~~L.4.L",
                      "SoilN.20Pool~~Standard.exchange.P.20",
                      "SoilN.20Pool~~Standard.Total.P.20",
                      "SoilN.20Pool~~All.Mycorhizae",
                      "SoilN.20Pool~~L.4.L"))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                        "Standard.exchange.P.20~~All.Mycorhizae",
                        "Standard.exchange.P.20~~L.4.L",
                        "Standard.Total.P.20~~All.Mycorhizae",
                        "Standard.Total.P.20~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "SoilN.20Pool~~Standard.exchange.P.20",
                        "SoilN.20Pool~~Standard.Total.P.20",
                        "SoilN.20Pool~~All.Mycorhizae",
                        "SoilN.20Pool~~L.4.L"))

# Explore individual model fits
sem.model.fits(mA)

coef.table2=sem.coefs(mA,data=byspecies,standardize="scale")
sem.coefs(mA,data=byspecies)

summary(m1)
summary(m2)
summary(m3)
summary(m4)
summary(m5)
summary(m6)

sem.plot(mA, byspecies,coef.table2, corr.errors = NULL,
         show.nonsig = TRUE, scaling = 15, alpha = 0.05)
View(coef.table2)
#######Biomass/Concentration

########
byspecies <- x[ which(x$Species=='Actaea'), ]

mfull=lme(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer) + L.4.L + P.Exch.Conc + All.Mycorhizae + Height + P.Total.Conc + SoilN, random = ~1|Site, data = byspecies,na.action=na.omit)
m=lme(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m1=lme(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer) + P.Exch.Conc, random = ~1|Site, data = byspecies,na.action=na.omit)
ma=lm(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer), data = byspecies,na.action=na.omit)
m1a=lm(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer)+ L.4.L, data = byspecies,na.action=na.omit)
summary(m1a)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m3=lme(L.4.L~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m4=lme(SoilN~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m5=lme(P.Total.Conc~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m6=lme(P.Exch.Conc~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)

mA = list(
  m1,
  m2,
  m3,
  m4,
  m5,
  m6
)

# Run goodness-of-fit tests
sem.fit(mA, byspecies,conditional=T,
        corr.errors=c("P.Exch.Conc~~P.Total.Conc",
                      "P.Exch.Conc~~All.Mycorhizae",
                      "P.Exch.Conc~~L.4.L",
                      "P.Total.Conc~~All.Mycorhizae",
                      "P.Total.Conc~~L.4.L",
                      "All.Mycorhizae~~L.4.L",
                      "SoilN~~P.Exch.Conc",
                      "SoilN~~P.Total.Conc",
                      "SoilN~~All.Mycorhizae",
                      "SoilN~~L.4.L"
        ))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("P.Exch.Conc~~P.Total.Conc",
                        "P.Exch.Conc~~All.Mycorhizae",
                        "P.Exch.Conc~~L.4.L",
                        "P.Total.Conc~~All.Mycorhizae",
                        "P.Total.Conc~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "SoilN~~P.Exch.Conc",
                        "SoilN~~P.Total.Conc",
                        "SoilN~~All.Mycorhizae",
                        "SoilN~~L.4.L"
          ))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("P.Exch.Conc~~P.Total.Conc",
                        "P.Exch.Conc~~All.Mycorhizae",
                        "P.Exch.Conc~~L.4.L",
                        "P.Total.Conc~~All.Mycorhizae",
                        "P.Total.Conc~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "SoilN~~P.Exch.Conc",
                        "SoilN~~P.Total.Conc",
                        "SoilN~~All.Mycorhizae",
                        "SoilN~~L.4.L"
          ))



coef.table2=sem.coefs(mA,data=byspecies,standardize="scale",
                      corr.errors=c("P.Exch.Conc~~P.Total.Conc",
                                    "P.Exch.Conc~~All.Mycorhizae",
                                    "P.Exch.Conc~~L.4.L",
                                    "P.Total.Conc~~All.Mycorhizae",
                                    "P.Total.Conc~~L.4.L",
                                    "All.Mycorhizae~~L.4.L",
                                    "SoilN~~P.Exch.Conc",
                                    "SoilN~~P.Total.Conc",
                                    "SoilN~~All.Mycorhizae",
                                    "SoilN~~L.4.L"
                      ))
summary(mfull)
summary(m2)
summary(m3)
summary(m4)
summary(m5)
summary(m6)

sem.plot(mA, byspecies,coef.table2, corr.errors = NULL,
         show.nonsig = TRUE, scaling = 15, alpha = 0.05)

##Aquilegia####################################################
byspecies <- x[ which(x$Species=='Aquilegia'), ]

mfull=lme(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer) + L.4.L + P.Exch.Conc + All.Mycorhizae + Height + P.Total.Conc + SoilN, random = ~1|Site, data = byspecies,na.action=na.omit)
m=lme(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m1=lme(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer) + L.4.L, random = ~1|Site, data = byspecies,na.action=na.omit)
ma=lm(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer), data = byspecies,na.action=na.omit)
m1a=lm(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer)+ L.4.L, data = byspecies,na.action=na.omit)
summary(m)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m3=lme(L.4.L~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m4=lme(SoilN~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m5=lme(P.Total.Conc~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m6=lme(P.Exch.Conc~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)

mA = list(
  mfull,
  m2,
  m3,
  m4,
  m5,
  m6
)

# Run goodness-of-fit tests
sem.fit(mA, byspecies,conditional=T,
        corr.errors=c("P.Exch.Conc~~P.Total.Conc",
                      "P.Exch.Conc~~All.Mycorhizae",
                      "P.Exch.Conc~~L.4.L",
                      "P.Total.Conc~~All.Mycorhizae",
                      "P.Total.Conc~~L.4.L",
                      "All.Mycorhizae~~L.4.L",
                      "SoilN~~P.Exch.Conc",
                      "SoilN~~P.Total.Conc",
                      "SoilN~~All.Mycorhizae",
                      "SoilN~~L.4.L"))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("P.Exch.Conc~~P.Total.Conc",
                        "P.Exch.Conc~~All.Mycorhizae",
                        "P.Exch.Conc~~L.4.L",
                        "P.Total.Conc~~All.Mycorhizae",
                        "P.Total.Conc~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "SoilN~~P.Exch.Conc",
                        "SoilN~~P.Total.Conc",
                        "SoilN~~All.Mycorhizae",
                        "SoilN~~L.4.L"))

# Explore individual model fits
sem.model.fits(mA)

coef.table2=sem.coefs(mA,data=byspecies,standardize="scale",
                      corr.errors=c("P.Exch.Conc~~P.Total.Conc",
                                    "P.Exch.Conc~~All.Mycorhizae",
                                    "P.Exch.Conc~~L.4.L",
                                    "P.Total.Conc~~All.Mycorhizae",
                                    "P.Total.Conc~~L.4.L",
                                    "All.Mycorhizae~~L.4.L",
                                    "SoilN~~P.Exch.Conc",
                                    "SoilN~~P.Total.Conc",
                                    "SoilN~~All.Mycorhizae",
                                    "SoilN~~L.4.L"))
sem.coefs(mA,data=byspecies)
summary(mfull)
summary(m2)
summary(m3)
summary(m4)
summary(m5)
summary(m6)

sem.plot(mA, byspecies,coef.table2, corr.errors = NULL,
         show.nonsig = TRUE, scaling = 15, alpha = 0.05)
##Cornus####################################################
byspecies <- x[ which(x$Species=='Cornus'), ]

mfull=lme(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer) + L.4.L + P.Exch.Conc + All.Mycorhizae + Height + P.Total.Conc + SoilN, random = ~1|Site, data = byspecies,na.action=na.omit)
m=lme(Dry.mass.above~Dry.Year.gm2+as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m1=lme(Dry.mass.above~Dry.Year.gm2+as.numeric(Deer) + L.4.L, random = ~1|Site, data = byspecies,na.action=na.omit)
m1a=lm(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer)+ L.4.L + P.Exch.Conc + All.Mycorhizae, data = byspecies,na.action=na.omit)
summary(ma)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m3=lme(L.4.L~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m4=lme(SoilN~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m5=lme(P.Total.Conc~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m6=lme(P.Exch.Conc~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)

mA = list(
  mfull,
  m2,
  m3,
  m4,
  m5,
  m6
)

# Run goodness-of-fit tests
sem.fit(mA, byspecies,conditional=T,
        corr.errors=c("P.Exch.Conc~~P.Total.Conc",
                      "P.Exch.Conc~~All.Mycorhizae",
                      "P.Exch.Conc~~L.4.L",
                      "P.Total.Conc~~All.Mycorhizae",
                      "P.Total.Conc~~L.4.L",
                      "All.Mycorhizae~~L.4.L",
                      "SoilN~~P.Exch.Conc",
                      "SoilN~~P.Total.Conc",
                      "SoilN~~All.Mycorhizae",
                      "SoilN~~L.4.L"))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("P.Exch.Conc~~P.Total.Conc",
                        "P.Exch.Conc~~All.Mycorhizae",
                        "P.Exch.Conc~~L.4.L",
                        "P.Total.Conc~~All.Mycorhizae",
                        "P.Total.Conc~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "SoilN~~P.Exch.Conc",
                        "SoilN~~P.Total.Conc",
                        "SoilN~~All.Mycorhizae",
                        "SoilN~~L.4.L"))

# Explore individual model fits
sem.model.fits(mA)

coef.table2=sem.coefs(mA,data=byspecies,standardize="scale")
sem.coefs(mA,data=byspecies)

summary(m1)
summary(m2)
summary(m3)
summary(m4)
summary(m5)
summary(m6)

sem.plot(mA, byspecies,coef.table2, corr.errors = NULL,
         show.nonsig = TRUE, scaling = 15, alpha = 0.05)

##Prenanthes####################################################
byspecies <- x[ which(x$Species=='Prenanthes'), ]

mfull=lme(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer) + L.4.L + P.Exch.Conc + All.Mycorhizae + Height + P.Total.Conc + SoilN, random = ~1|Site, data = byspecies,na.action=na.omit)
m=lme(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m1=lme(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer) + L.4.L + SoilN + Height, random = ~1|Site, data = byspecies,na.action=na.omit)
m1a=lm(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer)+ L.4.L + P.Exch.Conc + All.Mycorhizae, data = byspecies,na.action=na.omit)
summary(m1a)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m3=lme(L.4.L~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m4=lme(SoilN~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m5=lme(P.Total.Conc~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m6=lme(P.Exch.Conc~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)

mA = list(
  mfull,
  m2,
  m3,
  m4,
  m5,
  m6
)

# Run goodness-of-fit tests
sem.fit(mA, byspecies,conditional=T,
        corr.errors=c("P.Exch.Conc~~P.Total.Conc",
                      "P.Exch.Conc~~All.Mycorhizae",
                      "P.Exch.Conc~~L.4.L",
                      "P.Total.Conc~~All.Mycorhizae",
                      "P.Total.Conc~~L.4.L",
                      "All.Mycorhizae~~L.4.L",
                      "SoilN~~P.Exch.Conc",
                      "SoilN~~P.Total.Conc",
                      "SoilN~~All.Mycorhizae",
                      "SoilN~~L.4.L"))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("P.Exch.Conc~~P.Total.Conc",
                        "P.Exch.Conc~~All.Mycorhizae",
                        "P.Exch.Conc~~L.4.L",
                        "P.Total.Conc~~All.Mycorhizae",
                        "P.Total.Conc~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "SoilN~~P.Exch.Conc",
                        "SoilN~~P.Total.Conc",
                        "SoilN~~All.Mycorhizae",
                        "SoilN~~L.4.L"))

# Explore individual model fits
sem.model.fits(mA)

coef.table2=sem.coefs(mA,data=byspecies,standardize="scale")
sem.coefs(mA,data=byspecies)

summary(m1)
summary(m2)
summary(m3)
summary(m4)
summary(m5)
summary(m6)

sem.plot(mA, byspecies,coef.table2, corr.errors = NULL,
         show.nonsig = TRUE, scaling = 15, alpha = 0.05)

##Quercus####################################################
byspecies <- x[ which(x$Species=='Quercus'), ]

mfull=lme(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer) + L.4.L + P.Exch.Conc + All.Mycorhizae + Height + P.Total.Conc + SoilN, random = ~1|Site, data = byspecies,na.action=na.omit)
m=lme(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m1=lme(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer) + SoilN, random = ~1|Site, data = byspecies,na.action=na.omit)
m1a=lm(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer)+ L.4.L + P.Exch.Conc + All.Mycorhizae, data = byspecies,na.action=na.omit)
summary(m1a)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m3=lme(L.4.L~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m4=lme(SoilN~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m5=lme(P.Total.Conc~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m6=lme(P.Exch.Conc~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)

mA = list(
  mfull,
  m2,
  m3,
  m4,
  m5,
  m6
)

# Run goodness-of-fit tests
sem.fit(mA, byspecies,conditional=T,
        corr.errors=c("P.Exch.Conc~~P.Total.Conc",
                      "P.Exch.Conc~~All.Mycorhizae",
                      "P.Exch.Conc~~L.4.L",
                      "P.Total.Conc~~All.Mycorhizae",
                      "P.Total.Conc~~L.4.L",
                      "All.Mycorhizae~~L.4.L",
                      "SoilN~~P.Exch.Conc",
                      "SoilN~~P.Total.Conc",
                      "SoilN~~All.Mycorhizae",
                      "SoilN~~L.4.L"))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("P.Exch.Conc~~P.Total.Conc",
                        "P.Exch.Conc~~All.Mycorhizae",
                        "P.Exch.Conc~~L.4.L",
                        "P.Total.Conc~~All.Mycorhizae",
                        "P.Total.Conc~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "SoilN~~P.Exch.Conc",
                        "SoilN~~P.Total.Conc",
                        "SoilN~~All.Mycorhizae",
                        "SoilN~~L.4.L"))

# Explore individual model fits
sem.model.fits(mA)

coef.table2=sem.coefs(mA,data=byspecies,standardize="scale")
sem.coefs(mA,data=byspecies)

summary(m1)
summary(m2)
summary(m3)
summary(m4)
summary(m5)
summary(m6)

sem.plot(mA, byspecies,coef.table2, corr.errors = NULL,
         show.nonsig = TRUE, scaling = 15, alpha = 0.05)

#######Biomass/Pools
########
byspecies <- x[ which(x$Species=='Actaea'), ]

mfull=lme(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer) + L.4.L + Standard.exchange.P.20 + All.Mycorhizae + Height + Standard.Total.P.20 + SoilN.20Pool, random = ~1|Site, data = byspecies,na.action=na.omit)
m=lme(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m1=lme(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer) + L.4.L, random = ~1|Site, data = byspecies,na.action=na.omit)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m3=lme(L.4.L~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m4=lme(SoilN.20Pool~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m5=lme(Standard.Total.P.20~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m6=lme(Standard.exchange.P.20~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)

mA = list(
  mfull,
  m2,
  m3,
  m4,
  m5,
  m6
)

# Run goodness-of-fit tests
sem.fit(mA, byspecies,conditional=T,
        corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                      "Standard.exchange.P.20~~All.Mycorhizae",
                      "Standard.exchange.P.20~~L.4.L",
                      "Standard.Total.P.20~~All.Mycorhizae",
                      "Standard.Total.P.20~~L.4.L",
                      "All.Mycorhizae~~L.4.L",
                      "SoilN.20Pool~~Standard.exchange.P.20",
                      "SoilN.20Pool~~Standard.Total.P.20",
                      "SoilN.20Pool~~All.Mycorhizae",
                      "SoilN.20Pool~~L.4.L"
        ))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                        "Standard.exchange.P.20~~All.Mycorhizae",
                        "Standard.exchange.P.20~~L.4.L",
                        "Standard.Total.P.20~~All.Mycorhizae",
                        "Standard.Total.P.20~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "SoilN.20Pool~~Standard.exchange.P.20",
                        "SoilN.20Pool~~Standard.Total.P.20",
                        "SoilN.20Pool~~All.Mycorhizae",
                        "SoilN.20Pool~~L.4.L"
          ))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                        "Standard.exchange.P.20~~All.Mycorhizae",
                        "Standard.exchange.P.20~~L.4.L",
                        "Standard.Total.P.20~~All.Mycorhizae",
                        "Standard.Total.P.20~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "SoilN.20Pool~~Standard.exchange.P.20",
                        "SoilN.20Pool~~Standard.Total.P.20",
                        "SoilN.20Pool~~All.Mycorhizae",
                        "SoilN.20Pool~~L.4.L"
          ))



coef.table2=sem.coefs(mA,data=byspecies,standardize="scale",
                      corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                                    "Standard.exchange.P.20~~All.Mycorhizae",
                                    "Standard.exchange.P.20~~L.4.L",
                                    "Standard.Total.P.20~~All.Mycorhizae",
                                    "Standard.Total.P.20~~L.4.L",
                                    "All.Mycorhizae~~L.4.L",
                                    "SoilN.20Pool~~Standard.exchange.P.20",
                                    "SoilN.20Pool~~Standard.Total.P.20",
                                    "SoilN.20Pool~~All.Mycorhizae",
                                    "SoilN.20Pool~~L.4.L"
                      ))
summary(m1)
summary(m2)
summary(m3)
summary(m4)
summary(m5)
summary(m6)

sem.plot(mA, byspecies,coef.table2, corr.errors = NULL,
         show.nonsig = TRUE, scaling = 15, alpha = 0.05)

##Aquilegia####################################################
byspecies <- x[ which(x$Species=='Aquilegia'), ]

mfull=lme(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer) + L.4.L + Standard.exchange.P.20 + All.Mycorhizae + Height + Standard.Total.P.20 + SoilN.20Pool, random = ~1|Site, data = byspecies,na.action=na.omit)
m=lme(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m1=lme(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer) + L.4.L, random = ~1|Site, data = byspecies,na.action=na.omit)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m3=lme(L.4.L~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m4=lme(SoilN.20Pool~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m5=lme(Standard.Total.P.20~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m6=lme(Standard.exchange.P.20~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)

mA = list(
  mfull,
  m2,
  m3,
  m4,
  m5,
  m6
)

# Run goodness-of-fit tests
sem.fit(mA, byspecies,conditional=T,
        corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                      "Standard.exchange.P.20~~All.Mycorhizae",
                      "Standard.exchange.P.20~~L.4.L",
                      "Standard.Total.P.20~~All.Mycorhizae",
                      "Standard.Total.P.20~~L.4.L",
                      "All.Mycorhizae~~L.4.L",
                      "SoilN.20Pool~~Standard.exchange.P.20",
                      "SoilN.20Pool~~Standard.Total.P.20",
                      "SoilN.20Pool~~All.Mycorhizae",
                      "SoilN.20Pool~~L.4.L"))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                        "Standard.exchange.P.20~~All.Mycorhizae",
                        "Standard.exchange.P.20~~L.4.L",
                        "Standard.Total.P.20~~All.Mycorhizae",
                        "Standard.Total.P.20~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "SoilN.20Pool~~Standard.exchange.P.20",
                        "SoilN.20Pool~~Standard.Total.P.20",
                        "SoilN.20Pool~~All.Mycorhizae",
                        "SoilN.20Pool~~L.4.L"))

# Explore individual model fits
sem.model.fits(mA)

coef.table2=sem.coefs(mA,data=byspecies,standardize="scale",
                      corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                                    "Standard.exchange.P.20~~All.Mycorhizae",
                                    "Standard.exchange.P.20~~L.4.L",
                                    "Standard.Total.P.20~~All.Mycorhizae",
                                    "Standard.Total.P.20~~L.4.L",
                                    "All.Mycorhizae~~L.4.L",
                                    "SoilN.20Pool~~Standard.exchange.P.20",
                                    "SoilN.20Pool~~Standard.Total.P.20",
                                    "SoilN.20Pool~~All.Mycorhizae",
                                    "SoilN.20Pool~~L.4.L"))
sem.coefs(mA,data=byspecies)
summary(mfull)
summary(m2)
summary(m3)
summary(m4)
summary(m5)
summary(m6)

sem.plot(mA, byspecies,coef.table2, corr.errors = NULL,
         show.nonsig = TRUE, scaling = 15, alpha = 0.05)
##Cornus####################################################
byspecies <- x[ which(x$Species=='Cornus'), ]

mfull=lme(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer) + L.4.L + Standard.exchange.P.20 + All.Mycorhizae + Height + Standard.Total.P.20 + SoilN.20Pool, random = ~1|Site, data = byspecies,na.action=na.omit)
m=lme(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m1=lme(Dry.mass.above~Dry.Year.gm2+as.numeric(Deer) + L.4.L, random = ~1|Site, data = byspecies,na.action=na.omit)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m3=lme(L.4.L~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m4=lme(SoilN.20Pool~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m5=lme(Standard.Total.P.20~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m6=lme(Standard.exchange.P.20~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)

mA = list(
  mfull,
  m2,
  m3,
  m4,
  m5,
  m6
)

# Run goodness-of-fit tests
sem.fit(mA, byspecies,conditional=T,
        corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                      "Standard.exchange.P.20~~All.Mycorhizae",
                      "Standard.exchange.P.20~~L.4.L",
                      "Standard.Total.P.20~~All.Mycorhizae",
                      "Standard.Total.P.20~~L.4.L",
                      "All.Mycorhizae~~L.4.L",
                      "SoilN.20Pool~~Standard.exchange.P.20",
                      "SoilN.20Pool~~Standard.Total.P.20",
                      "SoilN.20Pool~~All.Mycorhizae",
                      "SoilN.20Pool~~L.4.L"))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                        "Standard.exchange.P.20~~All.Mycorhizae",
                        "Standard.exchange.P.20~~L.4.L",
                        "Standard.Total.P.20~~All.Mycorhizae",
                        "Standard.Total.P.20~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "SoilN.20Pool~~Standard.exchange.P.20",
                        "SoilN.20Pool~~Standard.Total.P.20",
                        "SoilN.20Pool~~All.Mycorhizae",
                        "SoilN.20Pool~~L.4.L"))

# Explore individual model fits
sem.model.fits(mA)

coef.table2=sem.coefs(mA,data=byspecies,standardize="scale")
sem.coefs(mA,data=byspecies)

summary(m1)
summary(m2)
summary(m3)
summary(m4)
summary(m5)
summary(m6)

sem.plot(mA, byspecies,coef.table2, corr.errors = NULL,
         show.nonsig = TRUE, scaling = 15, alpha = 0.05)

##Prenanthes####################################################
byspecies <- x[ which(x$Species=='Prenanthes'), ]

mfull=lme(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer) + L.4.L + Standard.exchange.P.20 + All.Mycorhizae + Height + Standard.Total.P.20 + SoilN.20Pool, random = ~1|Site, data = byspecies,na.action=na.omit)
m=lme(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m1=lme(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer) + L.4.L + Standard.Total.P.20 , random = ~1|Site, data = byspecies,na.action=na.omit)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m3=lme(L.4.L~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m4=lme(SoilN.20Pool~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m5=lme(Standard.Total.P.20~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m6=lme(Standard.exchange.P.20~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)

mA = list(
  mfull,
  m2,
  m3,
  m4,
  m5,
  m6
)

# Run goodness-of-fit tests
sem.fit(mA, byspecies,conditional=T,
        corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                      "Standard.exchange.P.20~~All.Mycorhizae",
                      "Standard.exchange.P.20~~L.4.L",
                      "Standard.Total.P.20~~All.Mycorhizae",
                      "Standard.Total.P.20~~L.4.L",
                      "All.Mycorhizae~~L.4.L",
                      "SoilN.20Pool~~Standard.exchange.P.20",
                      "SoilN.20Pool~~Standard.Total.P.20",
                      "SoilN.20Pool~~All.Mycorhizae",
                      "SoilN.20Pool~~L.4.L"))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                        "Standard.exchange.P.20~~All.Mycorhizae",
                        "Standard.exchange.P.20~~L.4.L",
                        "Standard.Total.P.20~~All.Mycorhizae",
                        "Standard.Total.P.20~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "SoilN.20Pool~~Standard.exchange.P.20",
                        "SoilN.20Pool~~Standard.Total.P.20",
                        "SoilN.20Pool~~All.Mycorhizae",
                        "SoilN.20Pool~~L.4.L"))

# Explore individual model fits
sem.model.fits(mA)

coef.table2=sem.coefs(mA,data=byspecies,standardize="scale")
sem.coefs(mA,data=byspecies)

summary(m1)
summary(m2)
summary(m3)
summary(m4)
summary(m5)
summary(m6)

sem.plot(mA, byspecies,coef.table2, corr.errors = NULL,
         show.nonsig = TRUE, scaling = 15, alpha = 0.05)

##Quercus####################################################
byspecies <- x[ which(x$Species=='Quercus'), ]

mfull=lme(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer) + L.4.L + Standard.exchange.P.20 + All.Mycorhizae + Height + Standard.Total.P.20 + SoilN.20Pool, random = ~1|Site, data = byspecies,na.action=na.omit)
m=lme(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m1=lme(Dry.mass.above~Dry.Year.gm2 + as.numeric(Deer) + Standard.exchange.P.20, random = ~1|Site, data = byspecies,na.action=na.omit)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m3=lme(L.4.L~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m4=lme(SoilN.20Pool~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m5=lme(Standard.Total.P.20~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)
m6=lme(Standard.exchange.P.20~Dry.Year.gm2 + as.numeric(Deer), random = ~1|Site, data = byspecies,na.action=na.omit)

mA = list(
  mfull,
  m2,
  m3,
  m4,
  m5,
  m6
)

# Run goodness-of-fit tests
sem.fit(mA, byspecies,conditional=T,
        corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                      "Standard.exchange.P.20~~All.Mycorhizae",
                      "Standard.exchange.P.20~~L.4.L",
                      "Standard.Total.P.20~~All.Mycorhizae",
                      "Standard.Total.P.20~~L.4.L",
                      "All.Mycorhizae~~L.4.L",
                      "SoilN.20Pool~~Standard.exchange.P.20",
                      "SoilN.20Pool~~Standard.Total.P.20",
                      "SoilN.20Pool~~All.Mycorhizae",
                      "SoilN.20Pool~~L.4.L"))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                        "Standard.exchange.P.20~~All.Mycorhizae",
                        "Standard.exchange.P.20~~L.4.L",
                        "Standard.Total.P.20~~All.Mycorhizae",
                        "Standard.Total.P.20~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "SoilN.20Pool~~Standard.exchange.P.20",
                        "SoilN.20Pool~~Standard.Total.P.20",
                        "SoilN.20Pool~~All.Mycorhizae",
                        "SoilN.20Pool~~L.4.L"))

# Explore individual model fits
sem.model.fits(mA)

coef.table2=sem.coefs(mA,data=byspecies,standardize="scale")
sem.coefs(mA,data=byspecies)

summary(m1)
summary(m2)
summary(m3)
summary(m4)
summary(m5)
summary(m6)

sem.plot(mA, byspecies,coef.table2, corr.errors = NULL,
         show.nonsig = TRUE, scaling = 15, alpha = 0.05)

####DSE###########
m1=lmer(DSE.binomial ~ Dry.Year.gm2 + as.numeric(Deer) + All.Mycorhizae+ P.Total.Conc + (1|Site), data = x)
m2=lmer(All.Mycorhizae ~ Dry.Year.gm2 + as.numeric(Deer)+ (1|Site) + Species,data = x)
m4=lmer(P.Total.Conc ~ Dry.Year.gm2 + as.numeric(Deer)+ (1|Site) + Species,data = x)
m5=lmer(P.Exch.Conc ~ Dry.Year.gm2 + as.numeric(Deer)+ (1|Site) + Species,data = x)
m6=lmer(SoilN ~ Dry.Year.gm2 + as.numeric(Deer)+ (1|Site) + Species,data = x)

mA = list(
  m1,
  m2,
  m4,
  m5,
  m6
)

# Run goodness-of-fit tests
sem.fit(mA,x1,conditional=T,
        corr.errors=c("P.Exch.Conc~~P.Total.Conc",
                      "P.Exch.Conc~~All.Mycorhizae",
                      "P.Total.Conc~~All.Mycorhizae",
                      "SoilN~~P.Exch.Conc",
                      "SoilN~~P.Total.Conc",
                      "SoilN~~All.Mycorhizae"))

# Obtain standardized regression coefficients
sem.coefs(mA,x1, standardize = "scale",
          corr.errors=c("P.Exch.Conc~~P.Total.Conc",
                        "P.Exch.Conc~~All.Mycorhizae",
                        "P.Total.Conc~~All.Mycorhizae",
                        "SoilN~~P.Exch.Conc",
                        "SoilN~~P.Total.Conc",
                        "SoilN~~All.Mycorhizae"))

# Explore individual model fits
sem.model.fits(mA)

coef.table2=sem.coefs(mA,data=x1,standardize="scale",
                      corr.errors=c("P.Exch.Conc~~P.Total.Conc",
                                    "P.Exch.Conc~~All.Mycorhizae",
                                    "P.Total.Conc~~All.Mycorhizae",
                                    "SoilN~~P.Exch.Conc",
                                    "SoilN~~P.Total.Conc",
                                    "SoilN~~All.Mycorhizae"))
sem.coefs(mA,data=x1)
summary(m1)
summary(m2)
summary(m4)
summary(m5)
summary(m6)

sem.plot(mA, x1,coef.table2, corr.errors = NULL,
         show.nonsig = TRUE, scaling = 15, alpha = 0.05)

m1=lmer(DSE ~ Dry.Year.gm2 + as.numeric(Deer) + 1|Site, data = x1)
summary(m1)

####DSE################
x1=read.csv("RootScans.June281.csv",header=T)

m1=lmer(DSE ~ Dry.Year.gm2 + as.numeric(Deer) + All.Mycorhizae+ P.Total.Conc + (1|Site), data = x1)
m2=lmer(All.Mycorhizae ~ Dry.Year.gm2 + as.numeric(Deer)+ (1|Site) + Species,data = x1)
m4=lmer(P.Total.Conc ~ Dry.Year.gm2 + as.numeric(Deer)+ (1|Site) + Species,data = x1)
m5=lmer(P.Exch.Conc ~ Dry.Year.gm2 + as.numeric(Deer)+ (1|Site) + Species,data = x1)
m6=lmer(SoilN ~ Dry.Year.gm2 + as.numeric(Deer)+ (1|Site) + Species,data = x1)

m=glmer(Zero ~ Dry.Year.gm2 + as.numeric(Deer)+ (1|Site) + (1|Species),family=binomial(link = "logit"), data = x1)
mA = list(
  m1,
  m2,
  m4,
  m5,
  m6
)

# Run goodness-of-fit tests
sem.fit(mA,x1,conditional=T,
        corr.errors=c("P.Exch.Conc~~P.Total.Conc",
                      "P.Exch.Conc~~All.Mycorhizae",
                      "P.Total.Conc~~All.Mycorhizae",
                      "SoilN~~P.Exch.Conc",
                      "SoilN~~P.Total.Conc",
                      "SoilN~~All.Mycorhizae"))

# Obtain standardized regression coefficients
sem.coefs(mA,x1, standardize = "scale",
          corr.errors=c("P.Exch.Conc~~P.Total.Conc",
                        "P.Exch.Conc~~All.Mycorhizae",
                        "P.Total.Conc~~All.Mycorhizae",
                        "SoilN~~P.Exch.Conc",
                        "SoilN~~P.Total.Conc",
                        "SoilN~~All.Mycorhizae"))

# Explore individual model fits
sem.model.fits(mA)

coef.table2=sem.coefs(mA,data=x1,standardize="scale",
                      corr.errors=c("P.Exch.Conc~~P.Total.Conc",
                                    "P.Exch.Conc~~All.Mycorhizae",
                                    "P.Total.Conc~~All.Mycorhizae",
                                    "SoilN~~P.Exch.Conc",
                                    "SoilN~~P.Total.Conc",
                                    "SoilN~~All.Mycorhizae"))
sem.coefs(mA,data=x1)
summary(m1)
summary(m2)
summary(m4)
summary(m5)
summary(m6)

sem.plot(mA, x1,coef.table2, corr.errors = NULL,
         show.nonsig = TRUE, scaling = 15, alpha = 0.05)

m1=lmer(DSE ~ Dry.Year.gm2 + as.numeric(Deer) + 1|Site, data = x1)
summary(m1)
