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

	# Load functions
	source("./Analyses/EffectSizeFunctions.3.0.R")

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

	# Write table 
	write.csv(round_df(Model.Fits, digits = 3), "./output/tables/export/Model.Fits.csv")

# 4. Model averaging of "lnRR", "lnCVR", "lnVr". 
#--------------------------------------------------------------------------------------------#

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
#--------------------------------------------------------------------------------------------#
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
		ma.coefslnRR[i,c(2:4)]<-averageParameter(parameter = ma.coefslnRR$Parameter[i], weight=lnRR.Set$wi, models=modelslnRR)
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

# 8. Generate unconditional estimates for levels of the categorical predictors
#--------------------------------------------------------------------------------------------#
	# Sample sizes for respective categories and their levels
	nExptLifeStage <- table(datObjects$data$ExptLifeStage)
	nManipType     <- table(datObjects$data$ManipType)
	nCatchup           <- table(datObjects$data$CatchUp)
	nSex                  <- table(datObjects$data$Sex)
	nAdultDiet        <- table(datObjects$data$AdultDiet)
	nPhyl                 <- table(datObjects$data$Phylum)

	nVals <- c(dim(datObjects$data)[1], nExptLifeStage, nAdultDiet, nPhyl, nManipType, nSex, nCatchup)
	
	# lnRR
	#---------------------------------------------------------------------------------------------------#
		solPostlnRR <- matrix(nrow = 1000, ncol = length(ma.coefslnRR$Parameter))
			for(i in 1:length(ma.coefslnRR$Parameter)){
				solPostlnRR[,i] <-averageParamDist(parameter = ma.coefslnRR$Parameter[i], weight=lnRR.Set$wi, models=modelslnRR)
			}
			colnames(solPostlnRR) <- ma.coefslnRR$Parameter

			# Check that this gives the same results to Table 1
			posterior.mode(solPostlnRR)

			# Get unconditional estimates for: Postnatal
			PostnatalsolPostlnRR <- as.mcmc(
				solPostlnRR[, "(Intercept)"] + 
				solPostlnRR[, "AdultDietR"] * (nAdultDiet[2] / sum(nAdultDiet)) + 
				solPostlnRR[, "PhylumVertebrate"] * (nPhyl[2]/sum(nPhyl)) + 
				solPostlnRR[, "ManipTypeQuantity"] * (nManipType[2] / sum(nManipType)) + 
				solPostlnRR[, "CatchUp2.No"] * (nCatchup[2] / sum(nCatchup)) + 
				solPostlnRR[, "CatchUp3.Yes"] * (nCatchup[3] / sum(nCatchup)) + 
				solPostlnRR[, "SexF"] * (nSex[2] / sum(nSex)) + 
				solPostlnRR[, "SexM"] * (nSex[3] / sum(nSex)) 
				)

			# second step is easy, now we have unconditional estimate

			PrenatalsolPostlnRR <- as.mcmc(PostnatalsolPostlnRR + solPostlnRR[, "ExptLifeStagePrenatal"])

			PrePostsolPostlnRR <- rbind(cbind(mean(PostnatalsolPostlnRR), HPDinterval(PostnatalsolPostlnRR)), cbind(mean(PrenatalsolPostlnRR), HPDinterval(PrenatalsolPostlnRR)))
			row.names(PrePostsolPostlnRR) <- c("Postnatal", "Prenatal")

		# Manipulation Type
			ManipQualitysolPostlnRR <- as.mcmc(
				solPostlnRR[, "(Intercept)"] + 
				solPostlnRR[, "AdultDietR"] * (nAdultDiet[2] / sum(nAdultDiet)) + 
				solPostlnRR[, "PhylumVertebrate"] * (nPhyl[2]/sum(nPhyl)) + 
				solPostlnRR[, "ExptLifeStagePrenatal"] * (nExptLifeStage[2] / sum(nExptLifeStage)) + 
				solPostlnRR[, "CatchUp2.No"] * (nCatchup[2] / sum(nCatchup)) + 
				solPostlnRR[, "CatchUp3.Yes"] * (nCatchup[3] / sum(nCatchup)) + 
				solPostlnRR[, "SexF"] * (nSex[2] / sum(nSex)) + 
				solPostlnRR[, "SexM"] * (nSex[3] / sum(nSex)) 
				)

			# second step is easy, now we have unconditional estimate

			ManipQuantitysolPostlnRR <- as.mcmc(ManipQualitysolPostlnRR + solPostlnRR[, "ManipTypeQuantity"])

			ManipsolPostlnRR <- rbind(cbind(mean(ManipQualitysolPostlnRR), HPDinterval(ManipQualitysolPostlnRR)), cbind(mean(ManipQuantitysolPostlnRR), HPDinterval(ManipQuantitysolPostlnRR)))
			row.names(ManipsolPostlnRR) <- c("Quality", "Quantity")

		# Adult Diet
			AdUnrestrsolPostlnRR <- as.mcmc(
				solPostlnRR[, "(Intercept)"] + 
				solPostlnRR[, "ExptLifeStagePrenatal"] * (nExptLifeStage[2] / sum(nExptLifeStage)) + 
				solPostlnRR[, "PhylumVertebrate"] * (nPhyl[2]/sum(nPhyl)) + 
				solPostlnRR[, "ManipTypeQuantity"] * (nManipType[2] / sum(nManipType)) + 
				solPostlnRR[, "CatchUp2.No"] * (nCatchup[2] / sum(nCatchup)) + 
				solPostlnRR[, "CatchUp3.Yes"] * (nCatchup[3] / sum(nCatchup)) + 
				solPostlnRR[, "SexF"] * (nSex[2] / sum(nSex)) + 
				solPostlnRR[, "SexM"] * (nSex[3] / sum(nSex)) 
				)

			# second step is easy, now we have unconditional estimate

			RestrictsolPostlnRR <- as.mcmc(AdUnrestrsolPostlnRR + solPostlnRR[, "AdultDietR"])

			AdultDietsolPostlnRR <- rbind(cbind(mean(AdUnrestrsolPostlnRR), HPDinterval(AdUnrestrsolPostlnRR)), cbind(mean(RestrictsolPostlnRR), HPDinterval(RestrictsolPostlnRR)))
			row.names(AdultDietsolPostlnRR) <- c("Unrestricted", "Restricted")
	
		# Catchup
			Catchup1solPostlnRR <- as.mcmc(
				solPostlnRR[, "(Intercept)"] + 
				solPostlnRR[, "ExptLifeStagePrenatal"] * (nExptLifeStage[2] / sum(nExptLifeStage)) + 
				solPostlnRR[, "PhylumVertebrate"] * (nPhyl[2]/sum(nPhyl)) + 
				solPostlnRR[, "ManipTypeQuantity"] * (nManipType[2] / sum(nManipType)) + 
				solPostlnRR[, "AdultDietR"] * (nAdultDiet[2] / sum(nAdultDiet)) + 
				solPostlnRR[, "SexF"] * (nSex[2] / sum(nSex)) + 
				solPostlnRR[, "SexM"] * (nSex[3] / sum(nSex)) 
				)

			# second step is easy, now we have unconditional estimate
				Catchup2solPostlnRR <- Catchup1solPostlnRR + solPostlnRR[, "CatchUp2.No"] 
				Catchup3solPostlnRR <- Catchup1solPostlnRR + solPostlnRR[, "CatchUp3.Yes"]
			
			CatchsolPostlnRR <- rbind(cbind(mean(Catchup1solPostlnRR), HPDinterval(Catchup1solPostlnRR)), cbind(mean(Catchup2solPostlnRR), HPDinterval(Catchup2solPostlnRR)), cbind(mean(Catchup3solPostlnRR), HPDinterval(Catchup3solPostlnRR)))
			row.names(CatchsolPostlnRR) <- c("Unknown", "No", "Yes")
		
		#Sex
			SexBothsolPostlnRR <- as.mcmc(
				solPostlnRR[, "(Intercept)"] + 
				solPostlnRR[, "ExptLifeStagePrenatal"] * (nExptLifeStage[2] / sum(nExptLifeStage)) + 
				solPostlnRR[, "PhylumVertebrate"] * (nPhyl[2]/sum(nPhyl)) + 
				solPostlnRR[, "ManipTypeQuantity"] * (nManipType[2] / sum(nManipType)) + 
				solPostlnRR[, "AdultDietR"] * (nAdultDiet[2] / sum(nAdultDiet)) + 
				solPostlnRR[, "CatchUp2.No"] * (nCatchup[2] / sum(nCatchup)) + 
				solPostlnRR[, "CatchUp3.Yes"] * (nCatchup[3] / sum(nCatchup))
				)

			# second step is easy, now we have unconditional estimate
				SexMsolPostlnRR <- SexBothsolPostlnRR + solPostlnRR[, "SexM"] 
				SexFsolPostlnRR <- SexBothsolPostlnRR + solPostlnRR[, "SexF"]
			 
			SexsolPostlnRR <- rbind(cbind(mean(SexBothsolPostlnRR), HPDinterval(SexBothsolPostlnRR)), cbind(mean(SexMsolPostlnRR), HPDinterval(SexMsolPostlnRR)), cbind(posterior.mode(SexFsolPostlnRR), HPDinterval(SexFsolPostlnRR)))
			row.names(SexsolPostlnRR) <- c("Both", "Male", "Female")
	
		#Phylum
			InvertsolPostlnRR <- as.mcmc(
				solPostlnRR[, "(Intercept)"] + 
				solPostlnRR[, "AdultDietR"] * (nAdultDiet[2] / sum(nAdultDiet)) + 
				solPostlnRR[, "ManipTypeQuantity"] * (nManipType[2]/sum(nManipType)) + 
				solPostlnRR[, "ExptLifeStagePrenatal"] * (nExptLifeStage[2] / sum(nExptLifeStage)) + 
				solPostlnRR[, "CatchUp2.No"] * (nCatchup[2] / sum(nCatchup)) + 
				solPostlnRR[, "CatchUp3.Yes"] * (nCatchup[3] / sum(nCatchup)) + 
				solPostlnRR[, "SexF"] * (nSex[2] / sum(nSex)) + 
				solPostlnRR[, "SexM"] * (nSex[3] / sum(nSex)) 
				)

			# second step is easy, now we have unconditional estimate

			VertsolPostlnRR <- as.mcmc(InvertsolPostlnRR + solPostlnRR[, "PhylumVertebrate"])

			PhylumsolPostlnRR <- rbind(cbind(mean(InvertsolPostlnRR), HPDinterval(InvertsolPostlnRR)), cbind(mean(VertsolPostlnRR), HPDinterval(VertsolPostlnRR)))
			row.names(PhylumsolPostlnRR) <- c("Invertebrate", "Vertebrate")	
	# lnCVR
	#---------------------------------------------------------------------------------------------------#
		solPostlnCVR <- matrix(nrow = 1000, ncol = length(ma.coefslnCVR$Parameter))
			for(i in 1:length(ma.coefslnCVR$Parameter)){
				solPostlnCVR[,i] <-averageParamDist(parameter = ma.coefslnCVR$Parameter[i], weight=lnCVR.Set$wi, models=modelslnCVR)
			}
			colnames(solPostlnCVR) <- ma.coefslnCVR$Parameter

			# Check that this gives the same results to Table 1
			posterior.mode(solPostlnCVR)

		# Get unconditional estimates for: Postnatal
			PostnatalCVR <- as.mcmc(
				solPostlnCVR[, "(Intercept)"] + 
				solPostlnCVR[, "AdultDietR"] * (nAdultDiet[2] / sum(nAdultDiet)) + 
				solPostlnCVR[, "PhylumVertebrate"] * (nPhyl[2]/sum(nPhyl)) + 
				solPostlnCVR[, "ManipTypeQuantity"] * (nManipType[2] / sum(nManipType)) + 
				solPostlnCVR[, "CatchUp2.No"] * (nCatchup[2] / sum(nCatchup)) + 
				solPostlnCVR[, "CatchUp3.Yes"] * (nCatchup[3] / sum(nCatchup)) + 
				solPostlnCVR[, "SexF"] * (nSex[2] / sum(nSex)) + 
				solPostlnCVR[, "SexM"] * (nSex[3] / sum(nSex)) 
				)

			# second step is easy, now we have unconditional estimate

			PrenatalCVR <- as.mcmc(PostnatalCVR + solPostlnCVR[, "ExptLifeStagePrenatal"])

			PrePostCVR <- rbind(cbind(mean(PostnatalCVR), HPDinterval(PostnatalCVR)), cbind(mean(PrenatalCVR), HPDinterval(PrenatalCVR)))
			row.names(PrePostCVR) <- c("Postnatal", "Prenatal")

		# Manipulation Type
			ManipQualityCVR <- as.mcmc(
				solPostlnCVR[, "(Intercept)"] + 
				solPostlnCVR[, "AdultDietR"] * (nAdultDiet[2] / sum(nAdultDiet)) + 
				solPostlnCVR[, "PhylumVertebrate"] * (nPhyl[2]/sum(nPhyl)) + 
				solPostlnCVR[, "ExptLifeStagePrenatal"] * (nExptLifeStage[2] / sum(nExptLifeStage)) + 
				solPostlnCVR[, "CatchUp2.No"] * (nCatchup[2] / sum(nCatchup)) + 
				solPostlnCVR[, "CatchUp3.Yes"] * (nCatchup[3] / sum(nCatchup)) + 
				solPostlnCVR[, "SexF"] * (nSex[2] / sum(nSex)) + 
				solPostlnCVR[, "SexM"] * (nSex[3] / sum(nSex)) 
				)

			# second step is easy, now we have unconditional estimate

			ManipQuantityCVR <- as.mcmc(ManipQualityCVR + solPostlnCVR[, "ManipTypeQuantity"])

			ManipCVR <- rbind(cbind(mean(ManipQualityCVR), HPDinterval(ManipQualityCVR)), cbind(mean(ManipQuantityCVR), HPDinterval(ManipQuantityCVR)))
			row.names(ManipCVR) <- c("Quality", "Quantity")

		# Adult Diet
			AdUnrestrCVR <- as.mcmc(
				solPostlnCVR[, "(Intercept)"] + 
				solPostlnCVR[, "ExptLifeStagePrenatal"] * (nExptLifeStage[2] / sum(nExptLifeStage)) + 
				solPostlnCVR[, "PhylumVertebrate"] * (nPhyl[2]/sum(nPhyl)) + 
				solPostlnCVR[, "ManipTypeQuantity"] * (nManipType[2] / sum(nManipType)) + 
				solPostlnCVR[, "CatchUp2.No"] * (nCatchup[2] / sum(nCatchup)) + 
				solPostlnCVR[, "CatchUp3.Yes"] * (nCatchup[3] / sum(nCatchup)) + 
				solPostlnCVR[, "SexF"] * (nSex[2] / sum(nSex)) + 
				solPostlnCVR[, "SexM"] * (nSex[3] / sum(nSex)) 
				)

			# second step is easy, now we have unconditional estimate

			RestrictCVR <- as.mcmc(AdUnrestrCVR + solPostlnCVR[, "AdultDietR"])

			AdultDietCVR <- rbind(cbind(mean(AdUnrestrCVR), HPDinterval(AdUnrestrCVR)), cbind(mean(RestrictCVR), HPDinterval(RestrictCVR)))
			row.names(AdultDietCVR) <- c("Unrestricted", "Restricted")
	
		# Catchup
			Catchup1CVR <- as.mcmc(
				solPostlnCVR[, "(Intercept)"] + 
				solPostlnCVR[, "ExptLifeStagePrenatal"] * (nExptLifeStage[2] / sum(nExptLifeStage)) + 
				solPostlnCVR[, "PhylumVertebrate"] * (nPhyl[2]/sum(nPhyl)) + 
				solPostlnCVR[, "ManipTypeQuantity"] * (nManipType[2] / sum(nManipType)) + 
				solPostlnCVR[, "AdultDietR"] * (nAdultDiet[2] / sum(nAdultDiet)) + 
				solPostlnCVR[, "SexF"] * (nSex[2] / sum(nSex)) + 
				solPostlnCVR[, "SexM"] * (nSex[3] / sum(nSex)) 
				)

			# second step is easy, now we have unconditional estimate
				Catchup2CVR <- Catchup1CVR + solPostlnCVR[, "CatchUp2.No"] 
				Catchup3CVR <- Catchup1CVR + solPostlnCVR[, "CatchUp3.Yes"]
			
			CatchCVR <- rbind(cbind(mean(Catchup1CVR), HPDinterval(Catchup1CVR)), cbind(mean(Catchup2CVR), HPDinterval(Catchup2CVR)), cbind(mean(Catchup3CVR), HPDinterval(Catchup3CVR)))
			row.names(CatchCVR) <- c("Unknown", "No", "Yes")
		
		#Sex
			SexBothCVR <- as.mcmc(
				solPostlnCVR[, "(Intercept)"] + 
				solPostlnCVR[, "ExptLifeStagePrenatal"] * (nExptLifeStage[2] / sum(nExptLifeStage)) + 
				solPostlnCVR[, "PhylumVertebrate"] * (nPhyl[2]/sum(nPhyl)) + 
				solPostlnCVR[, "ManipTypeQuantity"] * (nManipType[2] / sum(nManipType)) + 
				solPostlnCVR[, "AdultDietR"] * (nAdultDiet[2] / sum(nAdultDiet)) + 
				solPostlnCVR[, "CatchUp2.No"] * (nCatchup[2] / sum(nCatchup)) + 
				solPostlnCVR[, "CatchUp3.Yes"] * (nCatchup[3] / sum(nCatchup))
				)

			# second step is easy, now we have unconditional estimate
				SexMCVR <- SexBothCVR + solPostlnCVR[, "SexM"] 
				SexFCVR <- SexBothCVR + solPostlnCVR[, "SexF"]
			 
			SexCVR <- rbind(cbind(mean(SexBothCVR), HPDinterval(SexBothCVR)), cbind(mean(SexMCVR), HPDinterval(SexMCVR)), cbind(posterior.mode(SexFCVR), HPDinterval(SexFCVR)))
			row.names(SexCVR) <- c("Both", "Male", "Female")
	
		#Phylum
			InvertCVR <- as.mcmc(
				solPostlnCVR[, "(Intercept)"] + 
				solPostlnCVR[, "AdultDietR"] * (nAdultDiet[2] / sum(nAdultDiet)) + 
				solPostlnCVR[, "ManipTypeQuantity"] * (nManipType[2]/sum(nManipType)) + 
				solPostlnCVR[, "ExptLifeStagePrenatal"] * (nExptLifeStage[2] / sum(nExptLifeStage)) + 
				solPostlnCVR[, "CatchUp2.No"] * (nCatchup[2] / sum(nCatchup)) + 
				solPostlnCVR[, "CatchUp3.Yes"] * (nCatchup[3] / sum(nCatchup)) + 
				solPostlnCVR[, "SexF"] * (nSex[2] / sum(nSex)) + 
				solPostlnCVR[, "SexM"] * (nSex[3] / sum(nSex)) 
				)

			# second step is easy, now we have unconditional estimate

			VertCVR <- as.mcmc(InvertCVR + solPostlnCVR[, "PhylumVertebrate"])

			PhylumCVR <- rbind(cbind(mean(InvertCVR), HPDinterval(InvertCVR)), cbind(mean(VertCVR), HPDinterval(VertCVR)))
			row.names(PhylumCVR) <- c("Invertebrate", "Vertebrate")	
	# lnVR
	#---------------------------------------------------------------------------------------------------#
		solPostlnVR <- matrix(nrow = 1000, ncol = length(ma.coefslnVR$Parameter))
		for(i in 1:length(ma.coefslnVR$Parameter)){
			solPostlnVR[,i] <-averageParamDist(parameter = ma.coefslnVR$Parameter[i], weight=lnVR.Set$wi, models=modelslnVR)
		}
		colnames(solPostlnVR) <- ma.coefslnVR$Parameter

		# Check that this gives the same results to Table 1
		posterior.mode(solPostlnVR)

		# Get unconditional estimates for Postnatal
		
		Postnatal <- as.mcmc(solPostlnVR[, "(Intercept)"] + solPostlnVR[, "AdultDietR"] * (nAdultDiet[2] / sum(nAdultDiet)) + solPostlnVR[, "PhylumVertebrate"] * (nPhyl[2]/sum(nPhyl)))

		# second step is easy, now we have unconditional estimate

		Prenatal <- as.mcmc(Postnatal + solPostlnVR[, "ExptLifeStagePrenatal"])

		PrePost <- rbind(cbind(mean(Postnatal), HPDinterval(Postnatal)), cbind(mean(Prenatal), HPDinterval(Prenatal)))
		row.names(PrePost) <- c("Postnatal", "Prenatal")
		
		# A second example for Adult diet
		AdultUnrestricted <- as.mcmc(solPostlnVR[, "(Intercept)"] + solPostlnVR[, "ExptLifeStagePrenatal"] * (nExptLifeStage[2] / sum(nExptLifeStage)) + solPostlnVR[, "PhylumVertebrate"] * (nPhyl[2]/sum(nPhyl)))

		# second step is easy, now we have unconditional estimate
		AdultRestricted <- as.mcmc(AdultUnrestricted + solPostlnVR[, "AdultDietR"])

		AdultDiet <- rbind(cbind(mean(AdultUnrestricted), HPDinterval(AdultUnrestricted)), cbind(mean(AdultRestricted), HPDinterval(AdultRestricted)))
		row.names(AdultDiet) <- c("Unrestricted", "Restricted")
		
		# A second example for Adult diet
		Invertebrate <- as.mcmc(solPostlnVR[, "(Intercept)"] + solPostlnVR[, "ExptLifeStagePrenatal"] * (nExptLifeStage[2] / sum(nExptLifeStage)) + solPostlnVR[, "AdultDietR"] * (nAdultDiet[2] / sum(nAdultDiet)))

		# second step is easy, now we have unconditional estimate
		Vertebrate <- as.mcmc(Invertebrate + solPostlnVR[, "PhylumVertebrate"])

		Phylum <- rbind(cbind(mean(Invertebrate), HPDinterval(Invertebrate)), cbind(mean(Vertebrate), HPDinterval(Vertebrate)))
		row.names(Phylum) <- c("Invertebrate", "Vertebrate")


	# Make table with estimates
	lnRRCondEst <- rbind(CatchsolPostlnRR, SexsolPostlnRR, ManipsolPostlnRR, PhylumsolPostlnRR, AdultDietsolPostlnRR, PrePostsolPostlnRR)


	lnCVRCondEst <- rbind(CatchCVR, SexCVR, ManipCVR, PhylumCVR, AdultDietCVR, PrePostCVR)

	lnVRCondEst <- rbind(rep(NA, 3), rep(NA, 3), rep(NA, 3), rep(NA, 3), rep(NA, 3), rep(NA, 3), rep(NA, 3), rep(NA, 3), Phylum, AdultDiet, PrePost)

	

	#Overall 
			load("./output/models/lnRR.1.Rdata")
			TotlnRR <- cbind(as.numeric(mean(model$Sol[,"(Intercept)"])), HPDinterval(model$Sol[,"(Intercept)"]))

			load("./output/models/lnCVR.1.Rdata")
			TotlnCVR <- cbind(as.numeric(mean(model$Sol[,"(Intercept)"])), HPDinterval(model$Sol[,"(Intercept)"]))

			load("./output/models/lnVR.1.Rdata")
			TotalnVR <- cbind(as.numeric(mean(model$Sol[,"(Intercept)"])), HPDinterval(model$Sol[,"(Intercept)"]))
			Overall <- cbind(TotlnRR, TotlnCVR, TotalnVR)
	ESTTable <- round_df(cbind(lnRRCondEst, lnCVRCondEst, lnVRCondEst), digits = 3)
	ESTTable <- round_df(rbind(ESTTable, Overall), digits = 3)
	write.csv(ESTTable, "./output/tables/export/ESTTableFig1.csv")
		
# 10. Heterogeneity
#--------------------------------------------------------------------------------------------#
	#lnVR
		load("./output/models/lnVR.1.Rdata")
		modVR <- model
		summary(modVR)

		#weights and sampling error variance
		wi <- 1/datObjects$data$V.lnVR
		wVR  <- sum((wi) * (length(wi) - 1))  / (((sum(wi)^2) - sum((wi)^2)))

		styVR <- modVR$VCV[,"StudyNo"]
		rVR   <- modVR$VCV[,"units"]
		
		#Total heterogeneity 
		I2t      <-  (styVR + rVR) / (styVR + rVR + wVR)
		totlnVr <- cbind(posterior.mode(I2t), HPDinterval(I2t))

		#Study heterogeneity
		I2stdy <-  styVR / (styVR + rVR + wVR)
		studyVarVR <- posterior.mode(styVR + rVR) # Needed for Jensen's inequality
		StdylnVr <- cbind(posterior.mode(I2stdy), HPDinterval(I2stdy))

		# Increase Prenatal (pre) and adult restriction. Estimates from Table S1. Accounts for Jensen's Inequality. Estimates need to be exponentiated.
		prelnVR    <- exp(0.228 + 0.5*studyVarVR)
		adultlnVR <- exp(0.285 + 0.5*studyVarVR)

	#lnCVR
		load("./output/models/lnCVR.1.Rdata")
		modCVR <- model
		summary(modCVR)

		#weights and sampling error variance
		wi <- 1/datObjects$data$V.lnCVR
		wCVR  <- sum((wi) * (length(wi) - 1))  / (((sum(wi)^2) - sum((wi)^2)))

		styCVR <- modCVR$VCV[,"StudyNo"]
		rCVR    <- modCVR$VCV[,"units"]
		
		#Total heterogeneity 
		I2t      <-  (styCVR + rCVR) / (styCVR + rCVR + wCVR)
		totlnCVR <- cbind(posterior.mode(I2t), HPDinterval(I2t))

		#Study heterogeneity
		I2stdy <-  styCVR / (styCVR + rCVR + wCVR)
		studyVarCVR <- posterior.mode(styCVR + rCVR) # Needed for Jensen's inequality
		StdylnCVR <- cbind(posterior.mode(I2stdy), HPDinterval(I2stdy))

		# Increase Prenatal (pre) and adult restriction. Estimates from Table S1. Accounts for Jensen's Inequality.Estimates need to be exponentiated.
		prelnCVR   <- exp(0.404 + 0.5*studyVarCVR)
		adultlnCVR <- exp(0.344 + 0.5*studyVarCVR)
	
	#lnRR
		load("./output/models/lnRR.1.Rdata")
		modRR <- model
		summary(modRR)

		#weights and sampling error variance
		wiRR <- 1/datObjects$data$V.lnRR
		wRR  <- sum((wiRR) * (length(wiRR) - 1))  / (((sum(wiRR)^2) - sum((wiRR)^2)))

		styRR <- modRR$VCV[,"StudyNo"]
		rRR    <- modRR$VCV[,"units"]
		
		#Total heterogeneity 
		I2t      <-  (styRR + rRR) / (styRR + rRR + wRR)
		totlnRR <- cbind(posterior.mode(I2t), HPDinterval(I2t))

		#Study heterogeneity
		I2stdy <-  styRR / (styRR + rRR + wRR)
		StdylnRR <- cbind(posterior.mode(I2stdy), HPDinterval(I2stdy))
		stdyVarRR <- posterior.mode(styRR + rRR) #Jensen's Inequality

		# Increase Prenatal (pre) and adult restriction. Estimates from Table S1. Accounts for Jensen's Inequality.Estimates need to be exponentiated.
		prelnRR   <- exp(-0.224 + 0.5*stdyVarRR)

	# Table of heterogenetity 
		HeterTable <- cbind(rbind(StdylnVr, totlnVr), rbind(StdylnCVR, totlnCVR), rbind(StdylnRR, totlnRR))
		row.names(HeterTable) <- c("Study", "Total")
		colnames(HeterTable) <- rep(c("Est.", "LCI", "UCI"), 3)

		HeterTable <- round_df(HeterTable, digits = 3)
		write.csv(HeterTable, "./output/tables/export/HeterTable.csv")

# 11. Publication bias
#--------------------------------------------------------------------------------------------#
	#lnVR
		load("./output/models/lnVR.1.Rdata")
		modVR <- model
		summary(modVR)
		  res <- datObjects$data$lnVR - predict(modVR, marginal = ~Map)
		prec <- 1/datObjects$data$V.lnVR
		   W <- res*prec

		#Egger regression
		EggerVR <- lm(W ~ prec)
		summary(EggerVR)   # Intercept just border line non-significant

		# Trim and fill
		modlnVR <- rma(yi = res, vi = datObjects$data$V.lnVR)
		tflnVR      <- trimfill(modlnVR, estimator = "L0")
		tflnVR.R      <- trimfill(modlnVR, estimator = "R0")

	#lnCVR
		load("./output/models/lnCVR.1.Rdata")
		modCVR <- model
		summary(modCVR)
		  res <- datObjects$data$lnCVR - predict(modCVR, marginal = ~Map)
		prec <- 1/datObjects$data$V.lnCVR
		   W <- res*prec

		#Egger regression
		EggerCVR <- lm(W ~ prec)
		summary(EggerCVR)   # Intercept non-significant

		# Trim and fill
		modlnCVR <- rma(yi = res, vi = datObjects$data$V.lnCVR)
		tflnCVR      <- trimfill(modlnCVR)
	
	#lnRR
		load("./output/models/lnRR.1.Rdata")
		modRR <- model
		summary(modRR)
		  res <- datObjects$data$lnRR - predict(modRR, marginal = ~Map)
		prec <- 1/datObjects$data$V.lnRR
		   W <- res*prec

		#Egger regression
		EggerRR <- lm(W ~ prec)
		summary(EggerRR)   # Intercept non-significant

		# Trim and fill
		modlnRR <- rma(yi = res, vi = datObjects$data$V.lnRR)
		tflnRR      <- trimfill(modlnRR)

# 12. Analysis with lnSD
#-------------------------------------------------------------------------------#
    prior<-list(R=list(V=1, nu=0.002), G=list(G1=list(V=diag(2), nu=0.002, alpha.mu=c(0, 0), alpha.V=1000*diag(2))))
   	nitts  <-5000000
	burn<-3000000
	thins<-(itts - burn) / 1000

    modellnSD<-MCMCglmm(lnSD ~ Trt + lnMean, random=~us(1 + Trt):StudyNo, mev=datObjects$data.long$V.lnSD, data=datObjects$data.long, prior=prior, nitt=nitts, burnin=burn, thin=thins, verbose=F, pr=T)
    summary(modellnSD)


	fe <- round_df(rbind(cbind(mean(modellnSD$Sol[,1]), HPDinterval(modellnSD$Sol[,1])), cbind(mean(modellnSD$Sol[,2]), HPDinterval(modellnSD$Sol[,2])), cbind(mean(modellnSD$Sol[,3]), HPDinterval(modellnSD$Sol[,3]))), digits = 3)

	re <- round_df(cbind(posterior.mode(modellnSD$VCV), HPDinterval(modellnSD$VCV)), digits =3)
   
   write.csv(fe, "./output/tables/export/lnSDfe.csv")
   write.csv(re, "./output/tables/export/lnSDre.csv")

   # Some comparable models with metafor

    MLMR.Uncor<-rma.mv(yi=lnSD, V=V.lnSD, mods=~Trt + lnMean, random=list(~1|StudyNo, ~1|EffectID), data=datObjects$data.long)
    summary(MLMR.Uncor)

    MLMR.Cor<-rma.mv(yi=lnSD, V=V.lnSD, mods=~Trt + lnMean, random=list(~Trt|StudyNo), data=datObjects$data.long)
    summary(MLMR.Cor)