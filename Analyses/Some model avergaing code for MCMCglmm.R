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

# 2. Model parameters and priors
#--------------------------------------------------------------------------------------------#
	
	# Specify the number of iterations, burnin period and thinning interval for MCMCglmm
	itts  <-5000000
	burn<-3000000
	thins<-(itts - burn) / 1000

	# Specify the priors for MCMCglmm; note that the final random factor of the covariance matrices has been fixed - meta-analysis because it is the known sampling (co)variance matrix.
	prior<-list(R=list(V=1, nu=0),
	                  G=list(
			G1=list(V=1, nu=0.002, alpha.mu=0, alpha.V=1000),  	
			G2=list(V=1, fix=1)
			  )
	)

	prior2<-list(R=list(V=1, nu=0),
	                  G=list(G1 = list(V=1, nu=2, alpha.mu=0, alpha.V=1000)  	
			  )
	)

# 3. Extract variables create combinations of predictors for model averaging
#--------------------------------------------------------------------------------------------#

	# Specify the predictors and responses
	variables <-c("ExptLifeStage", "ManipType", "CatchUp", "Sex", "AdultDiet", "Phylum")
	responses<-c("lnRR", "lnCVR", "lnVr")

	# Create the combinations of predictor variables
	models<-"1"

	for(i in 1:length(variables)){
		models<-c(models, apply(combn(variables, i), 2, paste, collapse=" + "))
	}

	# Create the Model.Fits table to hold the DIC values for all models
	Model.Fits              <-as.data.frame(array(NA, c(length(models), length(responses) + 1)))
	names(Model.Fits) <-c("Model", responses)
	Model.Fits[,1]        <-models

# 4. Check model diagnostics of full model prior to model averaging
#     to test that iterations specified and, burnin etc is OK.
#--------------------------------------------------------------------------------------------#
	# Fit saturated models and do diagnostic checks. Run three chains.
	modSatlnCVR <- list()
	seed <- round(runif(min = 1,max = 100, 3))
	
	# Tests with lnCVR
	for(i in 1:3){
		set.seed(seed[i])
		modSatlnCVR[[i]] <- MCMCglmm(lnCVR ~ ExptLifeStage + ManipType + CatchUp + Sex + AdultDiet + Phylum, random = ~StudyNo + Map, ginverse = list(Map = AnivGlnCVR), data = data, prior = prior, nitt = itts, burnin = burn, thin = thins)
	}
	
	saveRDS(modSatlnCVR, file = "./output/modSatlnCVR")
	
	# Tests with lnVR
	modSatlnVR <- list()
	seed <- round(runif(min = 1,max = 100, 3))

	for(i in 1:3){
		set.seed(seed[i])
		modSatlnVR[[i]] <- MCMCglmm(lnVR ~ AdultDiet + ExptLifeStage + ManipType + Sex + CatchUp + Phylum, random = ~StudyNo + Map, ginverse = list(Map = AnivGlnVR), data = data, prior = prior, nitt = itts, burnin = burn, thin = thins)
	}
	 
	for(i in 1:3){
		set.seed(seed[i])
		modSatlnVR[[i]] <- MCMCglmm(lnVR ~ AdultDiet + ExptLifeStage + ManipType + Sex + CatchUp + Phylum, random = ~StudyNo, data = data, prior = prior2, nitt = itts, burnin = burn, thin = thins)
	}

	saveRDS(modSatlnVR, file = "./output/modSatlnVR")

	# Tests with lnRR
	modSatlnRR <- list()
	seed <- round(runif(min = 1,max = 100, 3))
	
	for(i in 1:3){
		set.seed(seed[i])
		modSatlnRR[[i]] <- MCMCglmm(lnRR ~ ExptLifeStage + ManipType + CatchUp + Sex + AdultDiet + Phylum, random = ~StudyNo + Map, ginverse = list(Map = AinvlnRR), data = data, prior = prior, nitt = itts, burnin = burn, thin = thins)
	}

	saveRDS(modSatlnRR, file = "./output/modSatlnRR")

# 4.1 Checking models
#---------------------------------------------------------------------------------------------#
	# Read in model objects. LnCVR
		modSatlnCVR <- readRDS("./output/modSatlnCVR")
		chainsCVR <- MCMC.chains("./output/modSatlnCVR")
	# Explore chains separately
		plot(modSatlnCVR[[2]])
	# Diagnostics on chains. Note 3 chains will be pooled for the autocorrelation and heidel diagnostics. 
	MCMC.diag(chainsCVR, cols=1:2) #cols refers to VCV matrix.
	

	# Read in model objects. LnRR
		modSatlnRR <- readRDS("./output/modSatlnRR")
		chainsRR <- MCMC.chains("./output/modSatlnRR")
	# Explore chains separately
		plot(modSatlnRR[[3]])
	# Diagnostics on chains. Note 3 chains will be pooled for the autocorrelation and heidel diagnostics. 
		MCMC.diag(chainsRR, cols=1:2) #cols refers to VCV matrix.

# 5. Model averaging of "lnRR", "lnCVR", "lnVr". 
#---------------------------------------------------------------------------------------------#

	# Run through the table running each model, extracting the DIC values and saving the model (these models can be reloaded for individual interpretation); The script also saves each model to be used later
	# NOTE: THIS WILL TAKE A LONG TIME
	for(i in 1:length(Model.Fits[,1])){
			
			# Run the models for lnRR
			model.form<-paste(responses[1], Model.Fits[i,1], sep=" ~ ")
			model<-MCMCglmm(as.formula(model.form), random = ~Article.ID + Consumer.Sp + Data.ID, data=data, ginverse=list(Data.ID = lnRR.AnivG), nitt=itts, burnin=burn, thin=thins, verbose=F, pr=T, prior=prior)
			Model.Fits[i,2]<-model$DIC
			
			# Save the model with a names that corresponds to the model formula
		 	name<-paste(responses[1], Model.Fits[i,1], "Rdata", sep=".")
			save(model, file=name)
			
			# Run the models for lnCVR
			model.form<-paste(responses[2], Model.Fits[i,1], sep=" ~ ")
			model<-MCMCglmm(as.formula(model.form), random = ~Article.ID + Consumer.Sp + Data.ID, data=data, ginverse=list(Data.ID = lnCVR.AnivG), nitt=itts, burnin=burn, thin=thins, verbose=F, pr=T, prior=prior)
			Model.Fits[i,3]<-model$DIC
			
			# Save the model with a names that corresponds to the model formula
			name<-paste(responses[2], Model.Fits[i,1], "Rdata", sep=".")
			save(model, file=name)	
			
			# Print the progress of the script
			print(i/length(Model.Fits[,1]) * 100)

	}

	# See the output
	Model.Fits

	# Calculate the delta DIC values (delta.DIC)
	Model.Fits$lnRR.delta.DIC   <-Model.Fits$lnRR - min(Model.Fits$lnRR)
	Model.Fits$lnCVR.delta.DIC<-Model.Fits$lnCVR - min(Model.Fits$lnCVR)
	Model.Fits$lnVr.delta.DIC    <-Model.Fits$lnVr - min(Model.Fits$lnVr)

	# Select models Based on DIC < 3
	lnRR.Set   <-Model.Fits[which(Model.Fits$lnRR.delta.DIC < 3), c(1,2,5)]
	lnCVR.Set<-Model.Fits[which(Model.Fits$lnCVR.delta.DIC < 3), c(1,3,6)]
	lnVr.Set    <-Model.Fits[which(Model.Fits$lnVr.delta.DIC < 3), c(1,4,7)]

	# Calculate weights for top model set
	lnRR.Set$wi   <-exp(-0.5 * lnRR.Set$lnRR.delta.DIC) / sum(exp(-0.5 * lnRR.Set$lnRR.delta.DIC))
	lnCVR.Set$wi<-exp(-0.5 * lnCVR.Set$lnCVR.delta.DIC) / sum(exp(-0.5 * lnCVR.Set$lnCVR.delta.DIC))
	lnVr.Set$wi    <-exp(-0.5 * lnVr.Set$lnVr.delta.DIC) / sum(exp(-0.5 * lnVr.Set$lnVr.delta.DIC))

	# Sort model set based on DDIC
	lnRR.Set<-lnRR.Set[order(lnRR.Set$lnRR.delta.DIC),]
	lnCVR.Set<-lnCVR.Set[order(lnCVR.Set$lnCVR.delta.DIC),]

# See the top model sets and the model weights
lnRR.Set
lnCVR.Set

# Load the model set for lnRR using the function above and a vector of file names
File.Names<-paste("lnRR", lnRR.Set$Model, sep=".")
models<-loadModels(File.Names)

# Create a table to hold the model averaged coefficients - take the coefficients names from the global model in the models list (if the global model did not come up, you will need to adjust this code) 
ma.coefs<-as.data.frame(array(NA, c(7, 4)))
names(ma.coefs)<-c("Parameter", "Mode", "LCI", "UCI")
ma.coefs$Parameter<-row.names(summary(models[[which(lnRR.Set == "Trait.Type + Habitat + Trophic.Level + Defence + Z.Mixed.Breadth")]])$solutions)

# Do the model averaging using the function above
for(i in 1:length(ma.coefs$Parameter)){
ma.coefs[i,c(2:4)]<-averageParameter(parameter = ma.coefs$Parameter[i], weight=lnRR.Set$wi, models=models)
}

# Here are the model averaged paramters; note that these are averaged using a method equivalent to the 'zero' method of model averaging
ma.coefs

# Repeat for lnCVR

# Load the model set for lnRR using the function above and a vector of file names
File.Names<-paste("lnCVR", lnCVR.Set$Model, sep=".")
models<-loadModels(File.Names)

# Create a table to hold the model averaged coefficents - take the coefficients names from the global model in the models list (if the global model did not come up, you will need to adjust this code) 
ma.coefs<-as.data.frame(array(NA, c(7, 4)))
names(ma.coefs)<-c("Parameter", "Mode", "LCI", "UCI")
ma.coefs$Parameter<-row.names(summary(models[[which(lnCVR.Set == "Trait.Type + Habitat + Trophic.Level + Defence + Z.Mixed.Breadth")]])$solutions)

# Do the model averaging using the function above
for(i in 1:length(ma.coefs$Parameter)){
ma.coefs[i,c(2:4)]<-averageParameter(parameter = ma.coefs$Parameter[i], weight=lnCVR.Set$wi, models=models)
}

# Here are the model averaged paramters; note that these are averaged using a method equivalent to the 'zero' method of model averaging
ma.coefs

