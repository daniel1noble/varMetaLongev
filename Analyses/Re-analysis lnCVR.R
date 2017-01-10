#--------------------------------------------------------------------------------------#
# Effect of diet restriction on variance in longevity
# Authors: AM Senior & D Noble
# Description: Processing data and create (co)variance matrices. Check
# priors, mixing and convergence of chains.
#--------------------------------------------------------------------------------------#

# 1. Packages and data processing
#--------------------------------------------------------------------------------------#
    # Clear working space
    rm(list = ls())

    # Load packages
    library(metafor)
    library(MCMCglmm)

    # Source functions. Working directory assumed to be the project folder main directory.
    source("./Analyses/EffectSizeFunctions.3.0.R")

    # Load data. Working directory assumed to be the project folder main directory.
    data<-read.csv("./data/Clean_Data_AS.csv")
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

# 2. Process data in long format for random regression with lnSD
#--------------------------------------------------------------------------------------#
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
#--------------------------------------------------------------------------------------#
    par(mfrow=c(1, 1))
    plot(data.long$lnMean, data.long$lnSD, col=data.long$Trt)

# 4. Create lnCVR matrix
#--------------------------------------------------------------------------------------#

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

    # Get correlation for control data alone without repeats
    mv.corr<-cor(log(data$MeanC[match(controls, data$Control.ID)]), log(data$SD_C[match(controls, data$Control.ID)]))

    # Calculate covariance values between lnCVR values at the positions in shared_list and place them on the matrix
        for (i in 1:dim(combinations)[1]){
          p1 <- combinations[i,1]
          p2 <- combinations[i,2]
          p1_p2_cov <- Calc.cov.lnCVR(CMean = data[p1,"MeanC"], CSD = data[p1,"SD_C"], CN = data[p1,"NStartControl"], mvcorr = mv.corr)
          VlnCVR[p1,p2] <- p1_p2_cov
          VlnCVR[p2,p1] <- p1_p2_cov

        }

    # Add in variance
        diag(VlnCVR)<-data$V.lnCVR

   # Solve matrix for MCMCglmm
        AnivGlnCVR <- as(solve(VlnCVR), "dgCMatrix")
    
# 5. Create lnVR matrix
#--------------------------------------------------------------------------------------#

    data$lnVR<-Calc.lnVR(CSD=data$SD_C, ESD=data$SD_E, CN=data$NStartControl, EN=data$NStartExp)

    data$V.lnVR<-Calc.var.lnVR(CN=data$NStartControl, EN=data$NStartExp)

    # Create square matrix
    VlnVR<-matrix(0,nrow = dim(data)[1],ncol = dim(data)[1])

    Data.ID<-seq(1, dim(data)[1], 1)
    rownames(VlnVR) <- Data.ID
    colnames(VlnVR) <- Data.ID

    Shared.Control<-data$Control.ID
    shared_coord <- which(Shared.Control%in%Shared.Control[duplicated(Shared.Control)]==TRUE)
    shared_coord

    # matrix of combinations of coordinates for each experiment with shared control
    combinations <- do.call("rbind", tapply(shared_coord, data[shared_coord,"Control.ID"], function(x) t(combn(x,2))))

    # Correlation between mean and sd
    mv.corr<-cor(log(data$MeanC[match(controls, data$Control.ID)]), log(data$SD_C[match(controls, data$Control.ID)]))

    # Use combinations of shared controls to calculated the covariance
    for (i in 1:dim(combinations)[1]){
      p1 <- combinations[i,1]
      p2 <- combinations[i,2]
      p1_p2_cov <- Calc.cov.lnVR(CN = data[p1,"NStartControl"])
      VlnVR[p1,p2] <- p1_p2_cov
      VlnVR[p2,p1] <- p1_p2_cov
    }

    # Add variance along diagonal
    diag(VlnVR)<-data$V.lnVR

    # Analyse using MLMA in MCMCglmm
    AnivGlnVR<-as(solve(VlnVR), "dgCMatrix")
    
# 6. Create lnRR matrix
#--------------------------------------------------------------------------------------#

    # New thought - if SD is affected then d could be biased so should repeat their averages analysis with lnRR which is less biased by variance

    data$lnRR<-Calc.lnRR(CMean = data$MeanC, EMean = data$MeanE)

    data$V.lnRR<-Calc.var.lnRR(CN=data$NStartControl, EN=data$NStartExp, CMean = data$MeanC, EMean = data$MeanE, CSD=data$SD_C, ESD=data$SD_E)

    # Create square matrix
    VlnRR<-matrix(0,nrow = dim(data)[1],ncol = dim(data)[1])

    Data.ID<-seq(1, dim(data)[1], 1)
    rownames(VlnRR) <- Data.ID
    colnames(VlnRR) <- Data.ID

    Shared.Control<-data$Control.ID
    shared_coord <- which(Shared.Control%in%Shared.Control[duplicated(Shared.Control)]==TRUE)

    # matrix of combinations of coordinates for each experiment with shared control
    combinations <- do.call("rbind", tapply(shared_coord, data[shared_coord,"Control.ID"], function(x) t(combn(x,2))))
    
    # Use combinations of shared controls to calculated the covariance
    for (i in 1:dim(combinations)[1]){
      p1 <- combinations[i,1]
      p2 <- combinations[i,2]
      p1_p2_cov <- Calc.cov.lnRR(CN = data[p1,"NStartControl"], CSD = data[p1,"SD_C"], CMean = data[p1,"MeanC"])
      VlnRR[p1,p2] <- p1_p2_cov
      VlnRR[p2,p1] <- p1_p2_cov
    }

    # Add in variance for lnRR along diagonal
    diag(VlnRR)<-data$V.lnRR

    # Solve matrix for MCMCglmm
    AinvlnRR <- as(solve(VlnRR), "dgCMatrix")
    data$Map <-row.names(AinvlnRR)
    
# 7. Save all relevant objects into a single list object for future use
#--------------------------------------------------------------------------------------#

    datObjects <- list(data = data, data.long = data.long, AnivGlnCVR=AnivGlnCVR, AnivGlnVR=AnivGlnVR, AinvlnRR=AinvlnRR, VlnCVR = VlnCVR, VlnVR = VlnVR, VlnRR = VlnRR)
    saveRDS(datObjects, file = "./data/datObjects")

# 8. Model parameters and priors. Test that they are sufficient
#--------------------------------------------------------------------------------------#
    
    # Specify the number of iterations, burnin period and thinning interval for MCMCglmm
    itts  <-5000000
    burn<-3000000
    thins<-(itts - burn) / 1000

    # Specify the priors for MCMCglmm; note that the final random factor of the covariance matrices has been fixed - meta-analysis because it is the known sampling (co)variance matrix.
    prior<-list(R=list(V=1, nu=0.002),
                      G=list(
            G1=list(V=1, nu=0.002, alpha.mu=0, alpha.V=1000),   
            G2=list(V=1, fix=1)
              )
    )

# 9. Check model diagnostics of full model prior to model averaging
#     to test that iterations specified and, burn in etc is OK.
#--------------------------------------------------------------------------------------#
    # Fit saturated models and do diagnostic checks. Run three chains.
    modSatlnCVR <- list()
    seed <- round(runif(min = 1,max = 100, 3))
    
    # Tests with lnCVR
    for(i in 1:3){
        set.seed(seed[i])
        modSatlnCVR[[i]] <- MCMCglmm(lnCVR ~ ExptLifeStage + ManipType + CatchUp + Sex + AdultDiet + Phylum, random = ~StudyNo + Map, ginverse = list(Map = datObjects$AnivGlnCVR), data = datObjects$data, prior = prior, nitt = itts, burnin = burn, thin = thins)
    }
    
    saveRDS(modSatlnCVR, file = "./output/modSatlnCVR")
    
    # Tests with lnVR
    modSatlnVR <- list()
    seed <- round(runif(min = 1,max = 100, 3))

    for(i in 1:3){
        set.seed(seed[i])
        modSatlnVR[[i]] <- MCMCglmm(lnVR ~ AdultDiet + ExptLifeStage + ManipType + Sex + CatchUp + Phylum, random = ~StudyNo + Map, ginverse = list(Map = datObjects$AnivGlnVR), data = datObjects$data, prior = prior, nitt = itts, burnin = burn, thin = thins)
    }
     
    saveRDS(modSatlnVR, file = "./output/modSatlnVR")

    # Tests with lnRR
    modSatlnRR <- list()
    seed <- round(runif(min = 1,max = 100, 3))
    
    for(i in 1:3){
        set.seed(seed[i])
        modSatlnRR[[i]] <- MCMCglmm(lnRR ~ ExptLifeStage + ManipType + CatchUp + Sex + AdultDiet + Phylum, random = ~StudyNo + Map, ginverse = list(Map = datObjects$AinvlnRR), data = datObjects$data, prior = prior, nitt = itts, burnin = burn, thin = thins)
    }

    saveRDS(modSatlnRR, file = "./output/modSatlnRR")

# 10. Checking models prior to model averaging
#--------------------------------------------------------------------------------------#
    # Read in model objects. LnCVR
        modSatlnCVR <- readRDS("./output/modSatlnCVR")
        chainsCVR      <- MCMC.chains("./output/modSatlnCVR", ginv = "Map")
    # Explore chains separately
        plot(modSatlnCVR[[2]])
    # Diagnostics on chains. Note 3 chains will be pooled for the autocorrelation and heidel diagnostics. 
        MCMC.diag(chainsCVR, cols=1:2) #cols refers to VCV matrix.
    
    # Read in model objects. LnCVR
        modSatlnVR <- readRDS("./output/modSatlnVR")
        chainsVR      <- MCMC.chains(path = "./output/modSatlnVR", ginv = "Map")
    # Explore chains separately
        plot(modSatlnVR[[2]])
    # Diagnostics on chains. Note 3 chains will be pooled for the autocorrelation and heidel diagnostics. 
        MCMC.diag(chainsVR, cols=1:2) #cols refers to VCV matrix.

    # Read in model objects. LnRR
        modSatlnRR <- readRDS("./output/modSatlnRR")
        chainsRR      <- MCMC.chains("./output/modSatlnRR", ginv = "Map")
    # Explore chains separately
        plot(modSatlnRR[[3]])
    # Diagnostics on chains. Note 3 chains will be pooled for the autocorrelation and heidel diagnostics. 
        MCMC.diag(chainsRR, cols=1:2) #cols refers to VCV matrix.



