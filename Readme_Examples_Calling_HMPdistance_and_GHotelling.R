rm(list=ls())
setwd("~/Desktop/Mitch/Bioinformatics/Final_Script")  ##User should modify to directory which hold all required files.

source('HMPdistance_function.r')

##Caculating HMP distances for reference stool sample and nasal samples.

d.l2.hmp <- HMPdistance(tax.level ='l2.phylum',
            d.new.filename = "HMP_L2_030317_phylum.csv", # this is what the user would specify as the filename of their dataset
            measure = 'bc',
            d.new.ix.col.not.rel.abu = 1:4,
            print.details = F)
head(d.l2.hmp)

##Caculating HMP reference distances for independent baxter data

d.l2.baxter <- HMPdistance(tax.level ='l2.phylum',
                    d.new.filename = "Baxter_L2_rel_abundance.csv", # this is what the user would specify as the filename of their dataset
                    d.new.ix.col.not.rel.abu = 1:4,
                    measure = 'corr',
                    print.details = F)
head(d.l2.baxter)

##Caculating Hoteling T-squared test
source("GHotelling_function.r")
load("test-mx.rdata")
data=st.na
nBoot=10000
dim(data)
GHotelling_function(data,nBoot)   
