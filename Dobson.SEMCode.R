library(lavaan)
library(qgraph)
library(semPlot)
library(lme4)
library(devtools)
install_github("jslefche/piecewiseSEM")
library(piecewiseSEM)
library(ape)
library(caper)
library(nlme)
library(glmm)
library(pscl)

x=read.csv("Dobson.Data.csv",header=T)

############Survival/Concentration models
########
byspecies <- x[ which(x$Species=='Actaea'), ]

m1=lme(Survival~Dry.Year.gm2 + as.factor(Deer) + L.4.L + P.Exch.Conc + All.Mycorhizae, random = ~1|Site, na.action = na.omit, data = byspecies)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m3=lme(L.4.L~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m4=lme(SoilN~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m5=lme(P.Total.Conc~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m6=lme(P.Exch.Conc~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)

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

m1=lme(Dry.mass.above~Dry.Year.gm2 + as.factor(Deer) + L.4.L + SoilN + P.Exch.Conc + All.Mycorhizae + Height, random = ~1|Site, na.action = na.omit, data = byspecies)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m3=lme(L.4.L~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m4=lme(SoilN~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m5=lme(P.Total.Conc~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m6=lme(P.Exch.Conc~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)

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
summary(m1)
summary(m2)
summary(m3)
summary(m4)
summary(m5)
summary(m6)

sem.plot(mA, byspecies,coef.table2, corr.errors = NULL,
         show.nonsig = TRUE, scaling = 15, alpha = 0.05)
##Cornus####################################################
byspecies <- x[ which(x$Species=='Cornus'), ]

m1=lme(Survival~Dry.Year.gm2*as.factor(Deer) + L.4.L + SoilN + P.Exch.Conc + All.Mycorhizae + Height, random = ~1|Site, na.action = na.omit, data = byspecies)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m3=lme(L.4.L~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m4=lme(SoilN~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m5=lme(P.Total.Conc~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m6=lme(P.Exch.Conc~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)

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

m1=lme(Survival~Dry.Year.gm2 + as.factor(Deer) + L.4.L + SoilN + P.Exch.Conc + All.Mycorhizae + Height + P.Total.Conc, random = ~1|Site, na.action = na.omit, data = byspecies)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m3=lme(L.4.L~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m4=lme(SoilN~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m5=lme(P.Total.Conc~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m6=lme(P.Exch.Conc~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)

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

m1=lme(Survival~Dry.Year.gm2 + as.factor(Deer) + All.Mycorhizae + L.4.L+ P.Exch.Conc + SoilN + P.Total.Conc + Height, random = ~1|Site, na.action = na.omit, data = byspecies)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m3=lme(L.4.L~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m4=lme(SoilN~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m5=lme(P.Total.Conc~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m6=lme(P.Exch.Conc~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)

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

#########Survival/Pools
########
byspecies <- x[ which(x$Species=='Actaea'), ]

m1=lme(Survival~Dry.Year.gm2 + as.factor(Deer) + L.4.L + Standard.exchange.P.20 + All.Mycorhizae, random = ~1|Site, na.action = na.omit, data = byspecies)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m3=lme(L.4.L~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m4=lme(Soil.N.by.exch.20~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m5=lme(Standard.Total.P.20~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m6=lme(Standard.exchange.P.20~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)

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
        corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                      "Standard.exchange.P.20~~All.Mycorhizae",
                      "Standard.exchange.P.20~~L.4.L",
                      "Standard.Total.P.20~~All.Mycorhizae",
                      "Standard.Total.P.20~~L.4.L",
                      "All.Mycorhizae~~L.4.L",
                      "Soil.N.by.exch.20~~Standard.exchange.P.20",
                      "Soil.N.by.exch.20~~Standard.Total.P.20",
                      "Soil.N.by.exch.20~~All.Mycorhizae",
                      "Soil.N.by.exch.20~~L.4.L"
        ))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                        "Standard.exchange.P.20~~All.Mycorhizae",
                        "Standard.exchange.P.20~~L.4.L",
                        "Standard.Total.P.20~~All.Mycorhizae",
                        "Standard.Total.P.20~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "Soil.N.by.exch.20~~Standard.exchange.P.20",
                        "Soil.N.by.exch.20~~Standard.Total.P.20",
                        "Soil.N.by.exch.20~~All.Mycorhizae",
                        "Soil.N.by.exch.20~~L.4.L"
          ))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                        "Standard.exchange.P.20~~All.Mycorhizae",
                        "Standard.exchange.P.20~~L.4.L",
                        "Standard.Total.P.20~~All.Mycorhizae",
                        "Standard.Total.P.20~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "Soil.N.by.exch.20~~Standard.exchange.P.20",
                        "Soil.N.by.exch.20~~Standard.Total.P.20",
                        "Soil.N.by.exch.20~~All.Mycorhizae",
                        "Soil.N.by.exch.20~~L.4.L"
          ))



coef.table2=sem.coefs(mA,data=byspecies,standardize="scale",
                      corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                                    "Standard.exchange.P.20~~All.Mycorhizae",
                                    "Standard.exchange.P.20~~L.4.L",
                                    "Standard.Total.P.20~~All.Mycorhizae",
                                    "Standard.Total.P.20~~L.4.L",
                                    "All.Mycorhizae~~L.4.L",
                                    "Soil.N.by.exch.20~~Standard.exchange.P.20",
                                    "Soil.N.by.exch.20~~Standard.Total.P.20",
                                    "Soil.N.by.exch.20~~All.Mycorhizae",
                                    "Soil.N.by.exch.20~~L.4.L"
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

m1=lme(Dry.mass.above~Dry.Year.gm2 + as.factor(Deer) + L.4.L + Soil.N.by.exch.20 + Standard.exchange.P.20 + All.Mycorhizae + Height, random = ~1|Site, na.action = na.omit, data = byspecies)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m3=lme(L.4.L~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m4=lme(Soil.N.by.exch.20~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m5=lme(Standard.Total.P.20~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m6=lme(Standard.exchange.P.20~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)

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
        corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                      "Standard.exchange.P.20~~All.Mycorhizae",
                      "Standard.exchange.P.20~~L.4.L",
                      "Standard.Total.P.20~~All.Mycorhizae",
                      "Standard.Total.P.20~~L.4.L",
                      "All.Mycorhizae~~L.4.L",
                      "Soil.N.by.exch.20~~Standard.exchange.P.20",
                      "Soil.N.by.exch.20~~Standard.Total.P.20",
                      "Soil.N.by.exch.20~~All.Mycorhizae",
                      "Soil.N.by.exch.20~~L.4.L"))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                        "Standard.exchange.P.20~~All.Mycorhizae",
                        "Standard.exchange.P.20~~L.4.L",
                        "Standard.Total.P.20~~All.Mycorhizae",
                        "Standard.Total.P.20~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "Soil.N.by.exch.20~~Standard.exchange.P.20",
                        "Soil.N.by.exch.20~~Standard.Total.P.20",
                        "Soil.N.by.exch.20~~All.Mycorhizae",
                        "Soil.N.by.exch.20~~L.4.L"))

# Explore individual model fits
sem.model.fits(mA)

coef.table2=sem.coefs(mA,data=byspecies,standardize="scale",
                      corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                                    "Standard.exchange.P.20~~All.Mycorhizae",
                                    "Standard.exchange.P.20~~L.4.L",
                                    "Standard.Total.P.20~~All.Mycorhizae",
                                    "Standard.Total.P.20~~L.4.L",
                                    "All.Mycorhizae~~L.4.L",
                                    "Soil.N.by.exch.20~~Standard.exchange.P.20",
                                    "Soil.N.by.exch.20~~Standard.Total.P.20",
                                    "Soil.N.by.exch.20~~All.Mycorhizae",
                                    "Soil.N.by.exch.20~~L.4.L"))
sem.coefs(mA,data=byspecies)
summary(m1)
summary(m2)
summary(m3)
summary(m4)
summary(m5)
summary(m6)

sem.plot(mA, byspecies,coef.table2, corr.errors = NULL,
         show.nonsig = TRUE, scaling = 15, alpha = 0.05)
##Cornus####################################################
byspecies <- x[ which(x$Species=='Cornus'), ]

m1=lme(Survival~Dry.Year.gm2*as.factor(Deer) + L.4.L + Soil.N.by.exch.20 + Standard.exchange.P.20 + All.Mycorhizae + Height, random = ~1|Site, na.action = na.omit, data = byspecies)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m3=lme(L.4.L~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m4=lme(Soil.N.by.exch.20~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m5=lme(Standard.Total.P.20~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m6=lme(Standard.exchange.P.20~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)

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
        corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                      "Standard.exchange.P.20~~All.Mycorhizae",
                      "Standard.exchange.P.20~~L.4.L",
                      "Standard.Total.P.20~~All.Mycorhizae",
                      "Standard.Total.P.20~~L.4.L",
                      "All.Mycorhizae~~L.4.L",
                      "Soil.N.by.exch.20~~Standard.exchange.P.20",
                      "Soil.N.by.exch.20~~Standard.Total.P.20",
                      "Soil.N.by.exch.20~~All.Mycorhizae",
                      "Soil.N.by.exch.20~~L.4.L"))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                        "Standard.exchange.P.20~~All.Mycorhizae",
                        "Standard.exchange.P.20~~L.4.L",
                        "Standard.Total.P.20~~All.Mycorhizae",
                        "Standard.Total.P.20~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "Soil.N.by.exch.20~~Standard.exchange.P.20",
                        "Soil.N.by.exch.20~~Standard.Total.P.20",
                        "Soil.N.by.exch.20~~All.Mycorhizae",
                        "Soil.N.by.exch.20~~L.4.L"))

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

m1=lme(Survival~Dry.Year.gm2 + as.factor(Deer) + L.4.L + Soil.N.by.exch.20 + Standard.exchange.P.20 + All.Mycorhizae + Height + Standard.Total.P.20, random = ~1|Site, na.action = na.omit, data = byspecies)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m3=lme(L.4.L~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m4=lme(Soil.N.by.exch.20~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m5=lme(Standard.Total.P.20~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m6=lme(Standard.exchange.P.20~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)

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
        corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                      "Standard.exchange.P.20~~All.Mycorhizae",
                      "Standard.exchange.P.20~~L.4.L",
                      "Standard.Total.P.20~~All.Mycorhizae",
                      "Standard.Total.P.20~~L.4.L",
                      "All.Mycorhizae~~L.4.L",
                      "Soil.N.by.exch.20~~Standard.exchange.P.20",
                      "Soil.N.by.exch.20~~Standard.Total.P.20",
                      "Soil.N.by.exch.20~~All.Mycorhizae",
                      "Soil.N.by.exch.20~~L.4.L"))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                        "Standard.exchange.P.20~~All.Mycorhizae",
                        "Standard.exchange.P.20~~L.4.L",
                        "Standard.Total.P.20~~All.Mycorhizae",
                        "Standard.Total.P.20~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "Soil.N.by.exch.20~~Standard.exchange.P.20",
                        "Soil.N.by.exch.20~~Standard.Total.P.20",
                        "Soil.N.by.exch.20~~All.Mycorhizae",
                        "Soil.N.by.exch.20~~L.4.L"))

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

m1=lme(Survival~Dry.Year.gm2 + as.factor(Deer) + All.Mycorhizae + L.4.L+ Standard.exchange.P.20 + Soil.N.by.exch.20 + Standard.Total.P.20 + Height, random = ~1|Site, na.action = na.omit, data = byspecies)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m3=lme(L.4.L~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m4=lme(Soil.N.by.exch.20~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m5=lme(Standard.Total.P.20~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m6=lme(Standard.exchange.P.20~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)

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
        corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                      "Standard.exchange.P.20~~All.Mycorhizae",
                      "Standard.exchange.P.20~~L.4.L",
                      "Standard.Total.P.20~~All.Mycorhizae",
                      "Standard.Total.P.20~~L.4.L",
                      "All.Mycorhizae~~L.4.L",
                      "Soil.N.by.exch.20~~Standard.exchange.P.20",
                      "Soil.N.by.exch.20~~Standard.Total.P.20",
                      "Soil.N.by.exch.20~~All.Mycorhizae",
                      "Soil.N.by.exch.20~~L.4.L"))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                        "Standard.exchange.P.20~~All.Mycorhizae",
                        "Standard.exchange.P.20~~L.4.L",
                        "Standard.Total.P.20~~All.Mycorhizae",
                        "Standard.Total.P.20~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "Soil.N.by.exch.20~~Standard.exchange.P.20",
                        "Soil.N.by.exch.20~~Standard.Total.P.20",
                        "Soil.N.by.exch.20~~All.Mycorhizae",
                        "Soil.N.by.exch.20~~L.4.L"))

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

#######Biomass/Concentration
setwd("~/Box Sync/amd338@cornell.edu - Annise Dobson's Files and Folders/Root architecture")
library(lavaan)
library(qgraph)
library(semPlot)
library(lme4)
library(devtools)
install_github("jslefche/piecewiseSEM")
library(piecewiseSEM)
library(ape)
library(caper)
library(nlme)
library(glmm)
library(pscl)

x=read.csv("RootScans.June28.csv",header=T)

####DSE
x1=read.csv("RootScans.June281.csv",header=T)

m1=lmer(DSE ~ Dry.Year.gm2 + as.factor(Deer) + All.Mycorhizae+ P.Total.Conc + (1|Site), data = x1)
m2=lmer(All.Mycorhizae ~ Dry.Year.gm2 + as.factor(Deer)+ (1|Site) + Species,data = x1)
m4=lmer(P.Total.Conc ~ Dry.Year.gm2 + as.factor(Deer)+ (1|Site) + Species,data = x1)
m5=lmer(P.Exch.Conc ~ Dry.Year.gm2 + as.factor(Deer)+ (1|Site) + Species,data = x1)
m6=lmer(SoilN ~ Dry.Year.gm2 + as.factor(Deer)+ (1|Site) + Species,data = x1)

m=glmer(Zero ~ Dry.Year.gm2 + as.factor(Deer)+ (1|Site) + (1|Species),family=binomial(link = "logit"), data = x1)
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

m1=lmer(DSE ~ Dry.Year.gm2 + as.factor(Deer) + 1|Site, na.action = na.omit, data = x1)
summary(m1)

########
byspecies <- x[ which(x$Species=='Actaea'), ]

m1=lme(Dry.mass.above~Dry.Year.gm2 + as.factor(Deer) + L.4.L + P.Exch.Conc + All.Mycorhizae, random = ~1|Site, na.action = na.omit, data = byspecies)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m3=lme(L.4.L~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m4=lme(SoilN~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m5=lme(P.Total.Conc~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m6=lme(P.Exch.Conc~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)

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

m1=lme(Dry.mass.above~Dry.Year.gm2 + as.factor(Deer) + L.4.L + SoilN + P.Exch.Conc + All.Mycorhizae + Height, random = ~1|Site, na.action = na.omit, data = byspecies)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m3=lme(L.4.L~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m4=lme(SoilN~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m5=lme(P.Total.Conc~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m6=lme(P.Exch.Conc~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)

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
summary(m1)
summary(m2)
summary(m3)
summary(m4)
summary(m5)
summary(m6)

sem.plot(mA, byspecies,coef.table2, corr.errors = NULL,
         show.nonsig = TRUE, scaling = 15, alpha = 0.05)
##Cornus####################################################
byspecies <- x[ which(x$Species=='Cornus'), ]

m1=lme(Dry.mass.above~Dry.Year.gm2*as.factor(Deer) + L.4.L + SoilN + P.Exch.Conc + All.Mycorhizae + Height, random = ~1|Site, na.action = na.omit, data = byspecies)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m3=lme(L.4.L~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m4=lme(SoilN~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m5=lme(P.Total.Conc~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m6=lme(P.Exch.Conc~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)

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

m1=lme(Dry.mass.above~Dry.Year.gm2 + as.factor(Deer) + L.4.L + SoilN + P.Exch.Conc + All.Mycorhizae + Height + P.Total.Conc, random = ~1|Site, na.action = na.omit, data = byspecies)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m3=lme(L.4.L~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m4=lme(SoilN~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m5=lme(P.Total.Conc~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m6=lme(P.Exch.Conc~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)

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

m1=lme(Dry.mass.above~Dry.Year.gm2 + as.factor(Deer) + All.Mycorhizae + L.4.L+ P.Exch.Conc + SoilN + P.Total.Conc + Height, random = ~1|Site, na.action = na.omit, data = byspecies)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m3=lme(L.4.L~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m4=lme(SoilN~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m5=lme(P.Total.Conc~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m6=lme(P.Exch.Conc~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)

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

m1=lme(Dry.mass.above~Dry.Year.gm2 + as.factor(Deer) + L.4.L + Standard.exchange.P.20 + All.Mycorhizae, random = ~1|Site, na.action = na.omit, data = byspecies)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m3=lme(L.4.L~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m4=lme(Soil.N.by.exch.20~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m5=lme(Standard.Total.P.20~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m6=lme(Standard.exchange.P.20~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)

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
        corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                      "Standard.exchange.P.20~~All.Mycorhizae",
                      "Standard.exchange.P.20~~L.4.L",
                      "Standard.Total.P.20~~All.Mycorhizae",
                      "Standard.Total.P.20~~L.4.L",
                      "All.Mycorhizae~~L.4.L",
                      "Soil.N.by.exch.20~~Standard.exchange.P.20",
                      "Soil.N.by.exch.20~~Standard.Total.P.20",
                      "Soil.N.by.exch.20~~All.Mycorhizae",
                      "Soil.N.by.exch.20~~L.4.L"
        ))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                        "Standard.exchange.P.20~~All.Mycorhizae",
                        "Standard.exchange.P.20~~L.4.L",
                        "Standard.Total.P.20~~All.Mycorhizae",
                        "Standard.Total.P.20~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "Soil.N.by.exch.20~~Standard.exchange.P.20",
                        "Soil.N.by.exch.20~~Standard.Total.P.20",
                        "Soil.N.by.exch.20~~All.Mycorhizae",
                        "Soil.N.by.exch.20~~L.4.L"
          ))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                        "Standard.exchange.P.20~~All.Mycorhizae",
                        "Standard.exchange.P.20~~L.4.L",
                        "Standard.Total.P.20~~All.Mycorhizae",
                        "Standard.Total.P.20~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "Soil.N.by.exch.20~~Standard.exchange.P.20",
                        "Soil.N.by.exch.20~~Standard.Total.P.20",
                        "Soil.N.by.exch.20~~All.Mycorhizae",
                        "Soil.N.by.exch.20~~L.4.L"
          ))



coef.table2=sem.coefs(mA,data=byspecies,standardize="scale",
                      corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                                    "Standard.exchange.P.20~~All.Mycorhizae",
                                    "Standard.exchange.P.20~~L.4.L",
                                    "Standard.Total.P.20~~All.Mycorhizae",
                                    "Standard.Total.P.20~~L.4.L",
                                    "All.Mycorhizae~~L.4.L",
                                    "Soil.N.by.exch.20~~Standard.exchange.P.20",
                                    "Soil.N.by.exch.20~~Standard.Total.P.20",
                                    "Soil.N.by.exch.20~~All.Mycorhizae",
                                    "Soil.N.by.exch.20~~L.4.L"
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

m1=lme(Dry.mass.above~Dry.Year.gm2 + as.factor(Deer) + L.4.L + Soil.N.by.exch.20 + Standard.exchange.P.20 + All.Mycorhizae + Height, random = ~1|Site, na.action = na.omit, data = byspecies)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m3=lme(L.4.L~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m4=lme(Soil.N.by.exch.20~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m5=lme(Standard.Total.P.20~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m6=lme(Standard.exchange.P.20~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)

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
        corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                      "Standard.exchange.P.20~~All.Mycorhizae",
                      "Standard.exchange.P.20~~L.4.L",
                      "Standard.Total.P.20~~All.Mycorhizae",
                      "Standard.Total.P.20~~L.4.L",
                      "All.Mycorhizae~~L.4.L",
                      "Soil.N.by.exch.20~~Standard.exchange.P.20",
                      "Soil.N.by.exch.20~~Standard.Total.P.20",
                      "Soil.N.by.exch.20~~All.Mycorhizae",
                      "Soil.N.by.exch.20~~L.4.L"))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                        "Standard.exchange.P.20~~All.Mycorhizae",
                        "Standard.exchange.P.20~~L.4.L",
                        "Standard.Total.P.20~~All.Mycorhizae",
                        "Standard.Total.P.20~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "Soil.N.by.exch.20~~Standard.exchange.P.20",
                        "Soil.N.by.exch.20~~Standard.Total.P.20",
                        "Soil.N.by.exch.20~~All.Mycorhizae",
                        "Soil.N.by.exch.20~~L.4.L"))

# Explore individual model fits
sem.model.fits(mA)

coef.table2=sem.coefs(mA,data=byspecies,standardize="scale",
                      corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                                    "Standard.exchange.P.20~~All.Mycorhizae",
                                    "Standard.exchange.P.20~~L.4.L",
                                    "Standard.Total.P.20~~All.Mycorhizae",
                                    "Standard.Total.P.20~~L.4.L",
                                    "All.Mycorhizae~~L.4.L",
                                    "Soil.N.by.exch.20~~Standard.exchange.P.20",
                                    "Soil.N.by.exch.20~~Standard.Total.P.20",
                                    "Soil.N.by.exch.20~~All.Mycorhizae",
                                    "Soil.N.by.exch.20~~L.4.L"))
sem.coefs(mA,data=byspecies)
summary(m1)
summary(m2)
summary(m3)
summary(m4)
summary(m5)
summary(m6)

sem.plot(mA, byspecies,coef.table2, corr.errors = NULL,
         show.nonsig = TRUE, scaling = 15, alpha = 0.05)
##Cornus####################################################
byspecies <- x[ which(x$Species=='Cornus'), ]

m1=lme(Dry.mass.above~Dry.Year.gm2*as.factor(Deer) + L.4.L + Soil.N.by.exch.20 + Standard.exchange.P.20 + All.Mycorhizae + Height, random = ~1|Site, na.action = na.omit, data = byspecies)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m3=lme(L.4.L~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m4=lme(Soil.N.by.exch.20~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m5=lme(Standard.Total.P.20~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m6=lme(Standard.exchange.P.20~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)

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
        corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                      "Standard.exchange.P.20~~All.Mycorhizae",
                      "Standard.exchange.P.20~~L.4.L",
                      "Standard.Total.P.20~~All.Mycorhizae",
                      "Standard.Total.P.20~~L.4.L",
                      "All.Mycorhizae~~L.4.L",
                      "Soil.N.by.exch.20~~Standard.exchange.P.20",
                      "Soil.N.by.exch.20~~Standard.Total.P.20",
                      "Soil.N.by.exch.20~~All.Mycorhizae",
                      "Soil.N.by.exch.20~~L.4.L"))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                        "Standard.exchange.P.20~~All.Mycorhizae",
                        "Standard.exchange.P.20~~L.4.L",
                        "Standard.Total.P.20~~All.Mycorhizae",
                        "Standard.Total.P.20~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "Soil.N.by.exch.20~~Standard.exchange.P.20",
                        "Soil.N.by.exch.20~~Standard.Total.P.20",
                        "Soil.N.by.exch.20~~All.Mycorhizae",
                        "Soil.N.by.exch.20~~L.4.L"))

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

m1=lme(Dry.mass.above~Dry.Year.gm2 + as.factor(Deer) + L.4.L + Soil.N.by.exch.20 + Standard.exchange.P.20 + All.Mycorhizae + Height + Standard.Total.P.20, random = ~1|Site, na.action = na.omit, data = byspecies)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m3=lme(L.4.L~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m4=lme(Soil.N.by.exch.20~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m5=lme(Standard.Total.P.20~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m6=lme(Standard.exchange.P.20~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)

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
        corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                      "Standard.exchange.P.20~~All.Mycorhizae",
                      "Standard.exchange.P.20~~L.4.L",
                      "Standard.Total.P.20~~All.Mycorhizae",
                      "Standard.Total.P.20~~L.4.L",
                      "All.Mycorhizae~~L.4.L",
                      "Soil.N.by.exch.20~~Standard.exchange.P.20",
                      "Soil.N.by.exch.20~~Standard.Total.P.20",
                      "Soil.N.by.exch.20~~All.Mycorhizae",
                      "Soil.N.by.exch.20~~L.4.L"))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                        "Standard.exchange.P.20~~All.Mycorhizae",
                        "Standard.exchange.P.20~~L.4.L",
                        "Standard.Total.P.20~~All.Mycorhizae",
                        "Standard.Total.P.20~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "Soil.N.by.exch.20~~Standard.exchange.P.20",
                        "Soil.N.by.exch.20~~Standard.Total.P.20",
                        "Soil.N.by.exch.20~~All.Mycorhizae",
                        "Soil.N.by.exch.20~~L.4.L"))

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

m1=lme(Dry.mass.above~Dry.Year.gm2 + as.factor(Deer) + All.Mycorhizae + L.4.L+ Standard.exchange.P.20 + Soil.N.by.exch.20 + Standard.Total.P.20 + Height, random = ~1|Site, na.action = na.omit, data = byspecies)
m2=lme(All.Mycorhizae~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m3=lme(L.4.L~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m4=lme(Soil.N.by.exch.20~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m5=lme(Standard.Total.P.20~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)
m6=lme(Standard.exchange.P.20~Dry.Year.gm2 + as.factor(Deer), random = ~1|Site, na.action = na.omit, data = byspecies)

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
        corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                      "Standard.exchange.P.20~~All.Mycorhizae",
                      "Standard.exchange.P.20~~L.4.L",
                      "Standard.Total.P.20~~All.Mycorhizae",
                      "Standard.Total.P.20~~L.4.L",
                      "All.Mycorhizae~~L.4.L",
                      "Soil.N.by.exch.20~~Standard.exchange.P.20",
                      "Soil.N.by.exch.20~~Standard.Total.P.20",
                      "Soil.N.by.exch.20~~All.Mycorhizae",
                      "Soil.N.by.exch.20~~L.4.L"))

# Obtain standardized regression coefficients
sem.coefs(mA,byspecies, standardize = "scale",
          corr.errors=c("Standard.exchange.P.20~~Standard.Total.P.20",
                        "Standard.exchange.P.20~~All.Mycorhizae",
                        "Standard.exchange.P.20~~L.4.L",
                        "Standard.Total.P.20~~All.Mycorhizae",
                        "Standard.Total.P.20~~L.4.L",
                        "All.Mycorhizae~~L.4.L",
                        "Soil.N.by.exch.20~~Standard.exchange.P.20",
                        "Soil.N.by.exch.20~~Standard.Total.P.20",
                        "Soil.N.by.exch.20~~All.Mycorhizae",
                        "Soil.N.by.exch.20~~L.4.L"))

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
m1=lmer(DSE.binomial ~ Dry.Year.gm2 + as.factor(Deer) + All.Mycorhizae+ P.Total.Conc + (1|Site), data = x)
m2=lmer(All.Mycorhizae ~ Dry.Year.gm2 + as.factor(Deer)+ (1|Site) + Species,data = x)
m4=lmer(P.Total.Conc ~ Dry.Year.gm2 + as.factor(Deer)+ (1|Site) + Species,data = x)
m5=lmer(P.Exch.Conc ~ Dry.Year.gm2 + as.factor(Deer)+ (1|Site) + Species,data = x)
m6=lmer(SoilN ~ Dry.Year.gm2 + as.factor(Deer)+ (1|Site) + Species,data = x)

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

m1=lmer(DSE ~ Dry.Year.gm2 + as.factor(Deer) + 1|Site, na.action = na.omit, data = x1)
summary(m1)