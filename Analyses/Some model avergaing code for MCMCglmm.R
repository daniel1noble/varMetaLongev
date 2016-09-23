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

	# Function to reload in models that will be averaged and return them in a list: these functions are available as loadable functions including comments in a more flexible format from the author
	loadModels<-function(model.names){
		
	# Create a list for models
		models<-list()
		
	# loop through the model names, loading the model and storing it in the list called models
		for(i in 1:length(model.names)){	
			name<-paste(model.names[i], ".Rdata", sep="")
			load(name)
			models[[i]]<-model
		}
		return(models)
	}

	# A function to average one of the parameters in the model set using an equivalent of the zero method. The function takes a paramter name, a list of models for avergaing and a vector of weights for each model 
	averageParameter<-function(parameter, weight, models){
		val<-0
		for(i in 1:length(models)){
			if(is.na(models[[i]]$Sol[, which(row.names(summary(models[[i]])$solutions) == parameter)][1]) == F){
			val<-val + models[[i]]$Sol[, which(row.names(summary(models[[i]])$solutions) == parameter)] * weight[i]
			}
		}
		return(c(posterior.mode(val), HPDinterval(val)))
	}	


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
					G2=list(V=1, nu=0.002, alpha.mu=0, alpha.V=1000), 	
					G3=list(V=1, fix=1)
					   )
				)

# 3. Extract variables create combinations of predictors for model averaging
#--------------------------------------------------------------------------------------------#

	# Specify the predictors and responses
	variables <-c("ExptLifeStage", "ManipType", "CatchUp", "Sex", "StageType", "AdultDiet", "Phylum")
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
	modSat <- list()
	for(i in 1:3){
		set.seed(seed[i])
		modSat[[i]] <- MCMCglmm(lnCVR ~ ExptLifeStage + ManipType + CatchUp + Sex + StageType + AdultDiet + Phylum, random = ~EffectID + StudyNo, ginverse = list(EffectID = AnivG), data = data, prior = prior, nitt = itts, burnin = burn, thin = thins)
	}
	
# 5. Model averaging of "lnRR", "lnCVR", "lnVr". 
#---------------------------------------------------------------------------------------------#
# MLMRs for all other variables; I have written a script to automate this process
# The script creates a table called Model.Fits which contains the DIC of each model
# The script also saves each model to be used later

# NOTE: the script will take some time to complete a

# Run through the table running each model, extracting the DIC values and saving the model (these models can be reloaded for individual interpretation)
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

