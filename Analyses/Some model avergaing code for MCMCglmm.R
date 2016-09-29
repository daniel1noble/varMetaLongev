#--------------------------------------------------------------------------------------------#
# Effect of diet restriction on variance in longevity: Analysis
# of Sinead and Uller data
# Authors: AM Senior & D Noble
# Description: Will model average MCMC model output using
# existing moderator variables, with known covariance structures
#--------------------------------------------------------------------------------------------#

# 1. Packages and workspace prep, functions
#--------------------------------------------------------------------------------------------#
	
	# Remove any old objects from the R environment
	rm(list=ls())

	# Load the packages
	library("metafor")
	library("MCMCglmm")

	# Load the list of model objects
	datObjects <- readRDS("./data/datObjects")

# 2. Model parameters and priors
#--------------------------------------------------------------------------------------------#
	
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

# 3. Extract variables create combinations of predictors for model averaging
#--------------------------------------------------------------------------------------------#

	# Specify the predictors and responses
	variables <-c("ExptLifeStage", "ManipType", "CatchUp", "Sex", "AdultDiet", "Phylum")
	responses<-c("lnRR", "lnCVR", "lnVR")

	# Create the combinations of predictor variables
	models<-"1"

	for(i in 1:length(variables)){
		models<-c(models, apply(combn(variables, i), 2, paste, collapse=" + "))
	}

	# Create the Model.Fits table to hold the DIC values for all models
	Model.Fits              <-as.data.frame(array(NA, c(length(models), length(responses) + 1)))
	names(Model.Fits) <-c("Model", responses)
	Model.Fits[,1]        <-models

# 4. Model averaging of "lnRR", "lnCVR", "lnVr". 
#---------------------------------------------------------------------------------------------#

	# Run through the table running each model, extracting the DIC values and saving the model (these models can be reloaded for individual interpretation); The script also saves each model to be used later
	for(i in 1:length(Model.Fits[,1])){
			
			# Run the models for lnRR
			model.form<-paste(responses[1], Model.Fits[i,1], sep=" ~ ")
			model<-MCMCglmm(as.formula(model.form),  random = ~StudyNo + Map, ginverse = list(Map = datObjects$AinvlnRR), data=datObjects$data, nitt=itts, burnin=burn, thin=thins, verbose=FALSE, pr=TRUE, prior=prior)
			Model.Fits[i,2]<-model$DIC
			
			# Save the model with a names that corresponds to the model formula
		 	name<-paste(responses[1], Model.Fits[i,1], "Rdata", sep=".")
			save(model, file=paste0("./output/models/", name))
			
			# Run the models for lnCVR
			model.form<-paste(responses[2], Model.Fits[i,1], sep=" ~ ")
			model<-MCMCglmm(as.formula(model.form), random = ~StudyNo + Map, ginverse = list(Map = datObjects$AnivGlnCVR), data = datObjects$data, nitt=itts, burnin=burn, thin=thins, verbose=FALSE, pr=TRUE, prior=prior)
			Model.Fits[i,3]<-model$DIC
			
			# Save the model with a names that corresponds to the model formula
			name<-paste(responses[2], Model.Fits[i,1], "Rdata", sep=".")
			save(model, file=paste0("./output/models/", name))

			# Run models for lnVR
			model.form<-paste(responses[3], Model.Fits[i,1], sep=" ~ ")
			model<-MCMCglmm(as.formula(model.form), random = ~StudyNo + Map, ginverse = list(Map = datObjects$AnivGlnVR), data = datObjects$data, nitt=itts, burnin=burn, thin=thins, verbose=FALSE, pr=TRUE, prior=prior)
			Model.Fits[i,4]<-model$DIC

			# Save the model with a names that corresponds to the model formula
			name<-paste(responses[3], Model.Fits[i,1], "Rdata", sep=".")
			save(model, file=paste0("./output/models/", name))	
			
			# Print the progress of the script
			print(i/length(Model.Fits[,1]) * 100)

	}

	# See the output
	Model.Fits <- readRDS("./output/Model.Fits")

	# Calculate the delta DIC values (delta.DIC)
	Model.Fits$lnRR.delta.DIC   <-Model.Fits$lnRR - min(Model.Fits$lnRR)
	Model.Fits$lnCVR.delta.DIC<-Model.Fits$lnCVR - min(Model.Fits$lnCVR)
	Model.Fits$lnVR.delta.DIC    <-Model.Fits$lnVR - min(Model.Fits$lnVR)

	# Select models Based on DIC < 3
	lnRR.Set   <-Model.Fits[which(Model.Fits$lnRR.delta.DIC < 3), c(1,2,5)]
	lnCVR.Set<-Model.Fits[which(Model.Fits$lnCVR.delta.DIC < 3), c(1,3,6)]
	lnVR.Set    <-Model.Fits[which(Model.Fits$lnVR.delta.DIC < 3), c(1,4,7)]

	# Calculate weights for top model set
	lnRR.Set$wi   <-exp(-0.5 * lnRR.Set$lnRR.delta.DIC) / sum(exp(-0.5 * lnRR.Set$lnRR.delta.DIC))
	lnCVR.Set$wi<-exp(-0.5 * lnCVR.Set$lnCVR.delta.DIC) / sum(exp(-0.5 * lnCVR.Set$lnCVR.delta.DIC))
	lnVR.Set$wi    <-exp(-0.5 * lnVR.Set$lnVR.delta.DIC) / sum(exp(-0.5 * lnVR.Set$lnVR.delta.DIC))

	# Sort model set based on DDIC
	lnRR.Set   <-lnRR.Set[order(lnRR.Set$lnRR.delta.DIC),]
	lnCVR.Set<-lnCVR.Set[order(lnCVR.Set$lnCVR.delta.DIC),]
	lnVR.Set   <-lnVR.Set[order(lnVR.Set$lnVR.delta.DIC),]

# 7. Processing and calculating model averages coefficients
#---------------------------------------------------------------------------------------------#
	# See the top model sets and the model weights
	lnRR.Set    <- round_df(lnRR.Set, digits = 3)
	lnCVR.Set <- round_df(lnCVR.Set, digits = 3)
	lnVR.Set    <- round_df(lnVR.Set, digits = 3)

	# Write these tables
	write.csv(lnRR.Set, file = "./output/tables/lnRR.Set.csv")
	write.csv(lnCVR.Set, file = "./output/tables/lnCVR.Set.csv")
	write.csv(lnVR.Set, file = "./output/tables/lnVR.Set.csv")

	# Load the model set for lnRR using the function above and a vector of file names
	File.NameslnRR<-paste0("./output/models/","lnRR.", lnRR.Set$Model)
	modelslnRR<-loadModels(File.NameslnRR)

	# Grab the unique parameters in the global model set. This can be used for all models. Anything not containing parameters should just be zero.
	global <- paste0("./output/models/","lnRR.", Model.Fits[64, "Model"])
	globalMod<-loadModels(global)

	# Create a table to hold the model averaged coefficients - take the coefficients names from the global model in the models list (if the global model did not come up, you will need to adjust this code) 
	ma.coefslnRR<-as.data.frame(array(NA, c(length(row.names(summary(globalMod[[1]])$solutions)), 4)))
	names(ma.coefslnRR)<-c("Parameter", "Mode", "LCI", "UCI")

	# Add these as row names, quick solution to getting global estimates
	ma.coefslnRR$Parameter<-row.names(summary(globalMod[[1]])$solutions)

	# Do the model averaging using the function above
	for(i in 1:length(ma.coefslnRR$Parameter)){
		ma.coefslnRR[i,c(2:4)]<-averageParameter(parameter = ma.coefslnRR$Parameter[i], weight=lnRR.Set$wi, models=modelslmRR)
	}

	# Here are the model averaged paramters; note that these are averaged using a method equivalent to the 'zero' method of model averaging
	ma.coefslnRR <- round_df(ma.coefslnRR, digits = 3)

	# Repeat for lnCVR

	# Load the model set for lnRR using the function above and a vector of file names
	File.NameslnCVR<-paste0("./output/models/","lnCVR.", lnCVR.Set$Model)
	modelslnCVR<-loadModels(File.NameslnCVR)

	ma.coefslnCVR<-as.data.frame(array(NA, c(length(row.names(summary(globalMod[[1]])$solutions)), 4)))
	names(ma.coefslnCVR)<-c("Parameter", "Mode", "LCI", "UCI")

	# Add these as row names, quick solution to getting global estimates
	ma.coefslnCVR$Parameter<-row.names(summary(globalMod[[1]])$solutions)

	# Do the model averaging using the function above
	for(i in 1:length(ma.coefslnCVR$Parameter)){
		ma.coefslnCVR[i,c(2:4)]<-averageParameter(parameter = ma.coefslnCVR$Parameter[i], weight=lnCVR.Set$wi, models=modelslnCVR)
	}

	# Here are the model averaged paramters; note that these are averaged using a method equivalent to the 'zero' method of model averaging
	ma.coefslnCVR <- round_df(ma.coefslnCVR, digits = 3)

	# Repeat for lnVR

	# Load the model set for lnRR using the function above and a vector of file names
	File.NameslnVR<-paste0("./output/models/","lnVR.", lnVR.Set$Model)
	modelslnVR<-loadModels(File.NameslnVR)

	# Create a table to hold the model averaged coefficents - take the coefficients names from the global model in the models list (if the global model did not come up, you will need to adjust this code) 
	ma.coefslnVR<-as.data.frame(array(NA, c(length(row.names(summary(globalMod[[1]])$solutions)), 4)))
	names(ma.coefslnVR)<-c("Parameter", "Mode", "LCI", "UCI")

	# Add these as row names, quick solution to getting global estimates
	ma.coefslnVR$Parameter<-row.names(summary(globalMod[[1]])$solutions)

	# Do the model averaging using the function above
	for(i in 1:length(ma.coefslnVR$Parameter)){
		ma.coefslnVR[i,c(2:4)]<-averageParameter(parameter = ma.coefslnVR$Parameter[i], weight=lnVR.Set$wi, models=modelslnVR)
	}

	# Here are the model averaged paramters; note that these are averaged using a method equivalent to the 'zero' method of model averaging
	ma.coefslnVR <- round_df(ma.coefslnVR, digits = 3)

	# cbind Tables together into a single one. 

	modAvgTable <- cBind(ma.coefslnRR, ma.coefslnCVR[-1], ma.coefslnVR[-1])	
	write.csv(modAvgTable, file = "./output/tables/modAvgTable.csv", row.names = FALSE)

# 8. Generate unconditional estimates in each of the levels of the categorical predictors
#------------------------------------------------------------------------------------------------------------
	# lnRR
	solPostlnRR <- matrix(nrow = 1000, ncol = length(ma.coefslnRR$Parameter))
			for(i in 1:length(ma.coefslnRR$Parameter)){
				solPostlnRR[,i] <-averageParamDist(parameter = ma.coefslnRR$Parameter[i], weight=lnRR.Set$wi, models=modelslnRR)
			}
			colnames(solPostlnRR) <- ma.coefslnRR$Parameter

			# Check that this gives the same results to Table 1
			posterior.mode(solPostlnRR)

	# lnCVR
	solPostlnCVR <- matrix(nrow = 1000, ncol = length(ma.coefslnCVR$Parameter))
			for(i in 1:length(ma.coefslnCVR$Parameter)){
				solPostlnCVR[,i] <-averageParamDist(parameter = ma.coefslnCVR$Parameter[i], weight=lnCVR.Set$wi, models=modelslnCVR)
			}
			colnames(solPostlnCVR) <- ma.coefslnCVR$Parameter

			# Check that this gives the same results to Table 1
			posterior.mode(solPostlnCVR)

	# lnVR
		solPostlnVR <- matrix(nrow = 1000, ncol = length(ma.coefslnVR$Parameter))
		for(i in 1:length(ma.coefslnVR$Parameter)){
			solPostlnVR[,i] <-averageParamDist(parameter = ma.coefslnVR$Parameter[i], weight=lnVR.Set$wi, models=modelslnVR)
		}
		colnames(solPostlnVR) <- ma.coefslnVR$Parameter

		# Check that this gives the same results to Table 1
		posterior.mode(solPostlnVR)

    # 9. Analysis with lnSD
#-------------------------------------------------------------------------------#

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

# 8. Analysis with metafor
#---------------------------------------------------------------------------------------------#
	# Analyse using REMA and MLMA in metafor
	        REMA<-rma.mv(yi = lnVR, V = VlnVR, random=~1|EffectID, data=data)
	        MLMA<-rma.mv(yi = lnVR, V = VlnVR, random=~1|StudyNo/EffectID, data=data)

	        anova(REMA, MLMA)

	        summary(REMA)
	        summary(MLMA)


	    # Check some moderators

	    MLMR<-rma.mv(yi = lnVR, V = V, mods=~AdultDiet + ExptLifeStage + ManipType + Sex + CatchUp + Phylum, random=~1|StudyNo/EffectID, data=data)
	    summary(MLMR)
	    
	    #model<-MCMCglmm(lnVR ~ AdultDiet + ExptLifeStage + ManipType + Sex + CatchUp + Phylum, random=~StudyNo + Map, ginverse=list(Map = AnivGlnVR), data=data, prior=prior, nitt=nitts, burnin=burn, thin=thins, verbose=T, pr=T)
	    #summary(model)
		# Analyse using REMA and MLMA in metafor
	    REMA<-rma.mv(yi = lnRR, V = VlnRR, random=~1|Map, data=data)
	    MLMA<-rma.mv(yi = lnRR, V = VlnRR, random=~1|StudyNo/Map, data=data)

	    anova(REMA, MLMA)

	    summary(REMA)
	    summary(MLMA)

	    MLMR<-rma.mv(yi = lnRR, V = VlnRR, mods=~AdultDiet + ExptLifeStage + ManipType + Sex + CatchUp + Phylum, random=~1|StudyNo/EffectID, data=data)
	    summary(MLMR)

	    # Results seem robust to d and lnRR
