#-------------------------------------------------------------------------#
# Effect of diet restriction on variance in longevity: Analysis
# of Sinead and Uller data
# Authors: AM Senior & D Noble
# Description: Processing data, running some random effects and
# multi-level (meta-regression) models using lnRR, lnCVR, lnVr
# and lnSD with lnMean in a random slopes model.
#---------------------------------------------------------------------------#

# 1. Packages and data processing
#--------------------------------------------------------#
    # Load packages
    library(metafor)
    library(MCMCglmm)

    # Source functions. Working directory assumed to be the project folder main directory.
    source("./Analyses/EffectSizeFunctions.3.0.R")

    # Load data. Working directory assumed to be the project folder main directory.
    data<-read.csv("./data/Clean_Data.csv")
    head(data)

    # Define variables
    data$StudyNo<-as.factor(data$StudyNo)
    data$EffectID<-as.factor(data$EffectID)

    # Re-classify "PreNatal-PostNatal" to "Prenatal" as was done in original analysis
    data$ExptLifeStage[data$ExptLifeStage=="PreNatal-PostNatal"] <- "Prenatal" 
    data$ExptLifeStage <- droplevels(data$ExptLifeStage)

    # Sorting out some shared control IDs
    data$Control.ID<-0

    studies<-unique(data$StudyNo)
    for(i in 1:length(studies)){
    	locations<-which(data$StudyNo == studies[i])
    	data$Control.ID[locations]<-seq(1, length(locations), 1)
    }

    data$Control.ID<-paste(data$StudyNo, data$Control.ID, sep="_")
    data$Control.ID[which(data$ShTrt > 0)]<-paste(data$StudyNo[which(data$ShTrt > 0)], data$ShTrt[which(data$ShTrt > 0)], sep="_")

    # A list of the unique control ids - useful later
    controls<-unique(data$Control.ID)

# 2. Process data in long format for Random regression with lnSD
#-------------------------------------------------------------------------------#

    # Creating a long format dataframe for the random-slopes model
    data.C<-data[,-c(16, 18, 20)]  
    data.C$Trt<-"C"
    data.E<-data[,-c(15, 17, 19)]
    data.E$Trt<-"E"
    names(data.C)[c(15,16,17)]<-c("N", "Mean", "SD")
    names(data.E)[c(15,16,17)]<-c("N", "Mean", "SD")

    # Drop the control repeats because control values were re-used multiple times and we want to create a stacked dataset
    data.C<-data.C[match(controls, data.C$Control.ID),]

    # Paste the two together the treatment and control means, errors etc
    data.long<-rbind(data.C, data.E)

    # Factorise the treatment variable
    data.long$Trt<-as.factor(data.long$Trt)

    # Calculate some effect sizes for later
    data.long$lnSD         <-Calc.lnSD(SD = data.long$SD, N = data.long$N)
    data.long$V.lnSD     <-Calc.var.lnSD(N = data.long$N)
    data.long$lnMean     <-log(data.long$Mean)
    data.long$V.lnMean <-(data.long$SD^2) / (data.long$N * data.long$Mean^2)

# 3. Explore mean-variance relationship in each group
#-----------------------------------------------------------------#

    par(mfrow=c(1, 1))
    plot(data.long$lnMean, data.long$lnSD, col=data.long$Trt)

# 4. Analysis of lnCVR
#-----------------------------------------------------------------#

    # Calculate lnCVRs
    data$lnCVR<-Calc.lnCVR(CMean = data$MeanC, EMean = data$MeanE, CSD=data$SD_C, ESD=data$SD_E, CN=data$NStartControl, EN=data$NStartExp)

    data$V.lnCVR<-Calc.var.lnCVR(CMean = data$MeanC, EMean = data$MeanE, CSD=data$SD_C, ESD=data$SD_E, CN=data$NStartControl, EN=data$NStartExp, repeated.control = T, Control.IDs = controls)

    # Create square matrix to estimate covariance for shared controls
    VlnCVR<-matrix(0,nrow = dim(data)[1],ncol = dim(data)[1])

    Data.ID<-seq(1, dim(data)[1], 1)
    rownames(VlnCVR) <- Data.ID
    colnames(VlnCVR) <- Data.ID

    Shared.Control<-data$Control.ID
    shared_coord <- which(Shared.Control%in%Shared.Control[duplicated(Shared.Control)]==TRUE)

    # matrix of combinations of coordinates for each experiment with shared control
    combinations <- do.call("rbind", tapply(shared_coord, data[shared_coord,"Control.ID"], function(x) t(combn(x,2))))
    combinations

    # Get correlation for control data alone without repeats
    mv.corr<-cor(log(data$MeanC[match(controls, data$Control.ID)]), log(data$SD_C[match(controls, data$Control.ID)]))

    # Function for getting covariance. @Noble - This is interesting. Not what I had expected to get the covariance for shared control. You'll need to explain this to me to clear up my understanding. I sort of expected this to be derived from the variance of the two VlnCV values for the effect sizes, but this is basically the variance of the control only. (eqn 12 - control).
    Calc.cov.lnCVR<-function(CMean, CSD, CN, mvcorr){
    	Cov<- CSD^2 / (CN * (CMean^2)) + 1 / (2 * (CN - 1)) - 2 * mvcorr * sqrt((CSD^2 / (CN * (CMean^2))) * (1 / (2 * (CN - 1)))) 
    	return(Cov)
    }

    # Calculate covariance values between lnCVR values at the positions in shared_list and place them on the matrix
    for (i in 1:dim(combinations)[1]){
      p1 <- combinations[i,1]
      p2 <- combinations[i,2]
      p1_p2_cov <- Calc.cov.lnCVR(CMean = data[p1,"MeanC"], CSD = data[p1,"SD_C"], CN = data[p1,"NStartControl"], mvcorr = mv.corr)
      VlnCVR[p1,p2] <- p1_p2_cov
      VlnCVR[p2,p1] <- p1_p2_cov

    }

# add in variance
diag(VlnCVR)<-data$V.lnCVR

# do REMA and MLMA in metafor. Note: EffectID needed for the V matrix.
REMA<-rma.mv(yi = lnCVR, V = VlnCVR, random=~1|EffectID, data=data)
MLMA<-rma.mv(yi = lnCVR, V = VlnCVR, random=~1|StudyNo/EffectID, data=data)
MLMA<-rma.mv(yi = lnCVR, V = VlnCVR, random=list(~1|StudyNo, ~1|EffectID), data=data) # Notation I would have used.

anova(REMA, MLMA)

summary(REMA)
summary(MLMA)

# Analyse using MLMA in MCMCglmm
mat<-solve(VlnCVR)
AnivGlnCVR<-as(mat, "dgCMatrix")

rownames(AnivGlnCVR) <- colnames(AnivGlnCVR) <- data$EffectID


prior<-list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002, alpha.mu=0, alpha.V=1000), G2=list(V=1, fix=1)))
nitts<-1000000
burn<-500000
thins<-(nitts - burn) / 1000

#model<-MCMCglmm(lnCVR ~ 1, random=~StudyNo + Map, ginverse=list(Map = AnivG), data=data, prior=prior, nitt=nitts, burnin=burn, thin=thins, verbose=F, pr=T)
#summary(model)

# Check some moderators

MLMR<-rma.mv(yi = lnCVR, V = V, mods=~AdultDiet + ExptLifeStage + ManipType + Sex + CatchUp + Phylum, random=~1|StudyNo/EffectID, data=data)
summary(MLMR)

# NOTES (apply to below also):
# 1) expt. stage has three levels in the dataset, but only 2 in the paper
# 2) not sure about global model approach e.g this has lower AICc - MLMR<-rma.mv(yi = lnCVR, V = V, mods=~AdultDiet + ExptLifeStage + Sex, random=~1|StudyNo/EffectID, data=data)

############ Repeat with lnVR

data$lnVR<-Calc.lnVR(CSD=data$SD_C, ESD=data$SD_E, CN=data$NStartControl, EN=data$NStartExp)

data$V.lnVR<-Calc.var.lnVR(CN=data$NStartControl, EN=data$NStartExp)

# Create square matrix
V<-matrix(0,nrow = dim(data)[1],ncol = dim(data)[1])

Data.ID<-seq(1, dim(data)[1], 1)
rownames(V) <- Data.ID
colnames(V) <- Data.ID

Shared.Control<-data$Control.ID
shared_coord <- which(Shared.Control%in%Shared.Control[duplicated(Shared.Control)]==TRUE)
shared_coord

# matrix of combinations of coordinates for each experiment with shared control
combinations <- do.call("rbind", tapply(shared_coord, data[shared_coord,"Control.ID"], function(x) t(combn(x,2))))
combinations


mv.corr<-cor(log(data$MeanC[match(controls, data$Control.ID)]), log(data$SD_C[match(controls, data$Control.ID)]))

Calc.cov.lnVR<-function(CN){	
	Cov<-(1 / (2 * (CN - 1))) 
	return(Cov)
}


for (i in 1:dim(combinations)[1]){

  p1 <- combinations[i,1]
  p2 <- combinations[i,2]
  p1_p2_cov <- Calc.cov.lnVR(CN = data[p1,"NStartControl"])
  V[p1,p2] <- p1_p2_cov
  V[p2,p1] <- p1_p2_cov

}

diag(V)<-data$V.lnVR

# Analyse using REMA and MLMA in metafor
REMA<-rma.mv(yi = lnVR, V = V, random=~1|EffectID, data=data)
MLMA<-rma.mv(yi = lnVR, V = V, random=~1|StudyNo/EffectID, data=data)

anova(REMA, MLMA)

summary(REMA)
summary(MLMA)

# Analyse using MLMA in MCMCglmm

mat<-solve(V)
AnivG<-as(mat, "dgCMatrix")

rownames(AnivG) = colnames(AnivG) <- data$EffectID
data$Map<-row.names(AnivG)

prior<-list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002, alpha.mu=0, alpha.V=1000), G2=list(V=1, fix=1)))
nitts<-1000000
burn<-500000
thins<-(nitts - burn) / 1000

#model<-MCMCglmm(lnVR ~ 1, random=~StudyNo + Map, ginverse=list(Map = AnivG), data=data, prior=prior, nitt=nitts, burnin=burn, thin=thins, verbose=F, pr=T)

summary(model)

# Check some moderators

MLMR<-rma.mv(yi = lnVR, V = V, mods=~AdultDiet + ExptLifeStage + ManipType + Sex + CatchUp + Phylum, random=~1|StudyNo/EffectID, data=data)
summary(MLMR)

############ lnSD

MLMR.Uncor<-rma.mv(yi=lnSD, V=V.lnSD, mods=~Trt + lnMean, random=~1|StudyNo/EffectID, data=data.long)
summary(MLMR)

MLMR.Cor<-rma.mv(yi=lnSD, V=V.lnSD, mods=~Trt + lnMean, random=~Trt|StudyNo/EffectID, data=data.long)
summary(MLMR)

prior<-list(R=list(V=1, nu=0.002), G=list(G1=list(V=diag(2), nu=0.002, alpha.mu=c(0, 0), alpha.V=1000*diag(2))))
nitts<-1000000
burn<-500000
thins<-(nitts - burn) / 1000

#model<-MCMCglmm(lnSD ~ Trt + lnMean, random=~us(1 + Trt):StudyNo, mev=data.long$V.lnSD, data=data.long, prior=prior, nitt=nitts, burnin=burn, thin=thins, verbose=F, pr=T)

summary(model)

plot(model$VCV)

#################

par(mfrow=c(1,2))

plot(data$lnVR, 1/sqrt(data$V.lnVR))
plot(data$lnCVR, 1/sqrt(data$V.lnCVR))

# New thought - if SD is affected then d could be biased so should repeat their averages analysis with lnRR which is less biased by variance

# Repeat with lnRR

data$lnRR<-Calc.lnRR(CMean = data$MeanC, EMean = data$MeanE)

data$V.lnRR<-Calc.var.lnRR(CN=data$NStartControl, EN=data$NStartExp, CMean = data$MeanC, EMean = data$MeanE, CSD=data$SD_C, ESD=data$SD_E)

# Create square matrix
V<-matrix(0,nrow = dim(data)[1],ncol = dim(data)[1])

Data.ID<-seq(1, dim(data)[1], 1)
rownames(V) <- Data.ID
colnames(V) <- Data.ID

Shared.Control<-data$Control.ID
shared_coord <- which(Shared.Control%in%Shared.Control[duplicated(Shared.Control)]==TRUE)

# matrix of combinations of coordinates for each experiment with shared control
combinations <- do.call("rbind", tapply(shared_coord, data[shared_coord,"Control.ID"], function(x) t(combn(x,2))))

Calc.cov.lnRR<-function(CN, CSD, CMean){
	Cov<-(CSD^2) / (CN * CMean^2)
	return(Cov)
}

for (i in 1:dim(combinations)[1]){

  p1 <- combinations[i,1]
  p2 <- combinations[i,2]
  p1_p2_cov <- Calc.cov.lnRR(CN = data[p1,"NStartControl"], CSD = data[p1,"SD_C"], CMean = data[p1,"MeanC"])
  V[p1,p2] <- p1_p2_cov
  V[p2,p1] <- p1_p2_cov

}

diag(V)<-data$V.lnRR

# Analyse using REMA and MLMA in metafor
REMA<-rma.mv(yi = lnRR, V = V, random=~1|EffectID, data=data)
MLMA<-rma.mv(yi = lnRR, V = V, random=~1|StudyNo/EffectID, data=data)

anova(REMA, MLMA)

summary(REMA)
summary(MLMA)

MLMR<-rma.mv(yi = lnRR, V = V, mods=~AdultDiet + ExptLifeStage + ManipType + Sex + CatchUp + Phylum, random=~1|StudyNo/EffectID, data=data)
summary(MLMR)

# Reuslts seem robust to d and lnRR

