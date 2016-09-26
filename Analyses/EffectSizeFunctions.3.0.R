
# A Collection of functions for calculating effect size 

# Functions were created by A. M Senior at various times (details below)
# and were combined in to this script on 13/06/2014 
# @ The University of Sydney. 

# The script was adapted on 29/04/2016, by Alistair Senior to include a function to simulate meta-analytic datasets for d-type (and lnRR/lnCVR) data.

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

	 # Functions for getting covariance for various effect sizes. @Noble - This is interesting. Not what I had expected to get the covariance for shared control. You'll need to explain this to me to clear up my understanding. I sort of expected this to be derived from the variance of the two VlnCV values for the effect sizes, but this is basically the variance of the control only. (eqn 12 - control).
	    Calc.cov.lnCVR<-function(CMean, CSD, CN, mvcorr){
	    	Cov<- CSD^2 / (CN * (CMean^2)) + 1 / (2 * (CN - 1)) - 2 * mvcorr * sqrt((CSD^2 / (CN * (CMean^2))) * (1 / (2 * (CN - 1)))) 
	    	return(Cov)
	    }


	Calc.cov.lnVR<-function(CN){	
	    	Cov<-(1 / (2 * (CN - 1))) 
	    	return(Cov)
	    }


	    Calc.cov.lnRR<-function(CN, CSD, CMean){
	    	Cov<-(CSD^2) / (CN * CMean^2)
	    	return(Cov)
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



#### Created by D Noble @ UNSW 23/09/2016
#Load an MCMC model object of multiple chains. 
MCMC.chains <- function(path, ginv = "EffectID"){
    imp <- readRDS(path)
	            #Import chains
		 VCV.list <- mcmc.list(lapply(imp, function(x) x$VCV[, -which(ginv == colnames(x$VCV))]))
		     Sol.list <- mcmc.list(lapply(imp, function(x) x$Sol))
		         DIC <- plyr::ldply(lapply(imp, function(x) x$DIC))
		#Combine chains
		     VCV.mode <- MCMCglmm::posterior.mode(as.mcmc(plyr::ldply(VCV.list)))
		      VCV.HPD <- coda::HPDinterval(as.mcmc(plyr::ldply(VCV.list)))

		       Sol.mode <- MCMCglmm::posterior.mode(as.mcmc(plyr::ldply(Sol.list)))
		      Sol.HPD <- coda::HPDinterval(as.mcmc(plyr::ldply(Sol.list)))
		
	return(list(solVCVlist = list(VCV = VCV.list, Sol = Sol.list ), solVCVchain = list(VCV = cbind(VCV.mode, VCV.HPD), Sol = cbind( Sol.mode, Sol.HPD)), DIC = DIC))
}

#### Created by D Noble @ UNSW 23/09/2016
## MCMC diagnostics function on a model object with 3 chains. 
MCMC.diag <- function(MCMC.Chains, cols = c(1:3)){
	gel.fixed <- gelman.diag(MCMC.Chains$solVCVlist$Sol)
	gel.VCV <- gelman.diag(MCMC.Chains$solVCVlist$VCV[,cols])
	autocor.fixed <- autocorr.diag(mcmc.list(MCMC.Chains$solVCVlist$Sol))
	autocor.VCV <- autocorr.diag(mcmc.list(MCMC.Chains$solVCVlist$VCV))
	heidel.fixed <- heidel.diag(mcmc.list(MCMC.Chains$solVCVlist$Sol))
	heidel.VCV <- heidel.diag(mcmc.list(MCMC.Chains$solVCVlist$VCV)[,cols])
  return(list(gel.fixed = gel.fixed, gel.VCV = gel.VCV, autocor.fixed = autocor.fixed, autocor.VCV= autocor.VCV, heidel.fixed = heidel.fixed, heidel.VCV= heidel.VCV))
}


#### Created by A M Senior @ The University of Otago 07/01/2014

##### Two functions to calculate d effect sizes for meta-analysis and standard error for d. Equations are taken from Nakagawa and Cuthill 2007 Biol. Reviews.
##### Calc.d calculates Hedges'd (sometimes g), which is 'unbiased for small sample size'. Cohen's d can be returned using adjusted = F.
##### Calc.SE.d calculates for d values given the necessery n's


Calc.d<-function(CMean, CSD, CN, EMean, ESD, EN, adjusted=T){
	
	
	varianceoftreatment<-(ESD)^2
	varianceofcontrols<-(CSD)^2

	sPooled<-sqrt(((EN - 1)*varianceoftreatment + (CN - 1)*varianceofcontrols) / (EN + CN - 2))	
	
	d<- (EMean - CMean) / sPooled 
	
	H.d<-d * (1 - (3 / (4 * (EN + CN - 2) - 1)))
	
	if(adjusted==T){return(H.d)}
	if(adjusted==F){return(d)}
}


Calc.SE.d<-function(CN, EN, d){
	
	SE<-sqrt(( (CN + EN) / (CN * EN) ) + ( (d^2) / (2 * (CN + EN - 2) ) ) )
	
	return(SE)
	
}




#### Created by A M Senior @ The University of Otago 07/01/2014

##### Two functions to calculate lnRR effect sizes for meta-analysis and var for lnRR.
##### Means for the two groups, sample sizes and standard deviations are required. 


Calc.lnRR<-function(CMean, EMean){
	
	lnRR<-log(EMean / CMean)
	
	return(lnRR)
	
}





Calc.var.lnRR<-function(CMean, CN, CSD, EMean, EN, ESD){
	
	EVar<-ESD^2
	CVar<-CSD^2
	
	V<- (CVar / (CN * (CMean^2))) + (EVar / (EN * (EMean^2)))
	
	return(V)
	
}


#### Created by A M Senior @ the University of Otago NZ 03/01/2014

#### Below are funcitons for calculating effect sizes for meta-analysis of variance - (modified on 16/05/2016 to account for repeated control data in calculation of correlation) 
#### Both functions take the mean, sd and n from the control and experimental groups.

#### The first function, Cal.lnCVR, calculates the the log repsonse-ratio of the coefficient of variance (lnCVR) - see Nakagawa et al in prep.

#### The second function calculates the measuremnt error variance for lnCVR. As well as the aforementioned parameters, this function also takes
#### Equal.E.C.Corr (default = T), which must be True or False. If true, the funciton assumes that the correlaiton between mean and sd (Taylor's Law) 
#### is equal for the mean and control groups, and, thus these data are pooled. If False the mean-SD correlation for the experimental and control groups
#### are calculated seperatley from one another.

#### A modification was added on 16/05/2016, by AM Senior at the university of Sydney, to the calc.var.lnCVR function, which allows for repeated control data.
#### By Default, this is ignored, but can be invoked by repeated.control = TRUE, in which case you need to pass a vector of IDs corresponding to those rows that contain the shared control.

Calc.lnCVR<-function(CMean, CSD, CN, EMean, ESD, EN){
	
	ES<-(log(ESD) - log(EMean) + 1 / (2*(EN - 1))) - (log(CSD) - log(CMean) + 1 / (2*(CN - 1)))
	
	return(ES)
	
}



Calc.var.lnCVR<-function(CMean, CSD, CN, EMean, ESD, EN, Equal.E.C.Corr=T, repeated.control = F, Control.IDs){
	

	
	if(repeated.control == T){
		
			mean.control.for.cor<-CMean[match(unique(Control.IDs), Control.IDs)]
			sd.control.for.cor<-CSD[match(unique(Control.IDs), Control.IDs)]
		
		}
		else{
		
			mean.control.for.cor<-CMean
			sd.control.for.cor<-CSD
			
		}
		
		if(Equal.E.C.Corr==T){
	
			mvcorr<-cor.test(log(c(mean.control.for.cor, EMean)), log(c(sd.control.for.cor, ESD)))$estimate
			
			S2<- CSD^2 / (CN * (CMean^2)) + 1 / (2 * (CN - 1)) - 2 * mvcorr * sqrt((CSD^2 / (CN * (CMean^2))) * (1 / (2 * (CN - 1)))) + ESD^2 / (EN * (EMean^2)) + 1 / (2 * (EN - 1)) - 2 * mvcorr * sqrt((ESD^2 / (EN * (EMean^2))) * (1 / (2 * (EN - 1))))
	
		}
		else{
		
			Cmvcorr<-cor.test(log(mean.control.for.cor), log(sd.control.for.cor))$estimate
			Emvcorr<-cor.test(log(EMean), (ESD))$estimate
	
			S2<- CSD^2 / (CN * (CMean^2)) + 1 / (2 * (CN - 1)) - 2 * Cmvcorr * sqrt((CSD^2 / (CN * (CMean^2))) * (1 / (2 * (CN - 1)))) + ESD^2 / (EN * (EMean^2)) + 1 / (2 * (EN - 1)) - 2 * Emvcorr * sqrt((ESD^2 / (EN * (EMean^2))) * (1 / (2 * (EN - 1))))		
		
		
		}
	
	return(S2)
	
}




#### Created by A M Senior @ The University of Sydney, Australia, 18/05/2014

### Below are functions for caluculating log variance estimates of variance and it's associated measurement error.

### The first function requires the estimate of standard deviation and a sample size, the second just a sample size


Calc.lnSD<-function(SD, N){
	
	lnSD <- log(SD) + (1 / (2 * (N - 1)))
	
	return(lnSD)
}



Calc.var.lnSD<-function(N){
	
	var.lnSD <- (1 / (2 * (N - 1)))
	
	return(var.lnSD)
	
}



#### Created by A M Senior @ The University of Sydney, Australia, 29/01/2016

### Below are functions for caluculating log variance ratio estimates of variance and it's associated measurement error.

### The first function requires the estimate of standard deviation and a sample size for each group, the second just a sample size for each group


Calc.lnVR<-function(ESD, EN, CSD, CN){
	
	lnVR <- log(ESD / CSD) + (1 / (2 * (EN - 1))) - (1 / (2 * (CN - 1)))
	
	return(lnVR)
}



Calc.var.lnVR<-function(EN, CN){
	
	var.lnVR <- (1 / (2 * (CN - 1))) + (1 / (2 * (EN - 1)))
	
	return(var.lnVR)
	
}

################## Created by A M Senior on the 29/04/2016 @ the University of Sydney

### Below is a function for simulating meta-analytic datasets, for a d/lnRR/lnCVR type meta-analysis, based on a normally distributed response.

### The function is flexible, and allows for the dataset to include heterogeneity and/or heteroscedasiticty, as well as simulaiton of publication bias.
### The minimum arguments that must be specified are control.mean, treatment.mean, SD and n.
### These correspond to the mean of the control and treatment groups, the SD (of both groups - homogeneity of variance) and n of both groups.
### By default 50 effect sizes are simulated, although this can be manipulated by altering k.
### By specifying control.SD and treatment.SD one can manipulate the SD in each group.
### var.n induces variance in the sample sizes by adding or subtracting (with equal probability) an integer drawn from a random poisson distribution, with lamda of the specified value - all smaple sizes have a lower bound of 2.
### control.n and treatment.n allow you to specify different average sample sizes in both groups.
### b_study.t.M abd b_study.c.M allow you to control heterogenetiy in group means, by adding a value to study specific group means, drawn from a random-normal distribution with mean = 0 and SD = these values (default is 0).
### b_study.t.SD abd b_study.c.SD allow you to control heterogenetiy in group SDs, by adding a value to study specific group means, drawn from a random folded-normal distribution with mean = 0 and SD = these values (default is 0).
### PB incorporates aspects of publication bias. If true, k significant effects are generated. All generated non-significant effects are also included, but for each effect size is given its significance and p.value.
### Subsequent algorithms can then be applied to subset the dataset based on some probability oif publicaiton criteria that is a function of significance.
### return.true.effects, allows the study specific 'true population statistics' to be returned - note if heterogeneity is absent these will be equivelent to those specified in the arguments - by default these data are not given.
### return.d specifies whether you want hedges d and its SE back as well


sim.d.data<-function(control.mean, SD, treatment.mean, n, k = 50, control.SD = SD, treatment.SD = SD, var.n = 0, control.n = n, treatment.n = n, b_study.t.M = 0, b_study.c.M = 0, b_study.t.SD = 0, b_study.c.SD = 0, PB = FALSE, return.true.effects = FALSE, return.d = TRUE){

	# Objects to store data
	ds<-NA
	SE.ds<-NA
	means.c<-NA
	SD.c<-NA
	n.c<-NA
	means.e<-NA
	SD.e<-NA
	n.e<-NA
	tests<-NA
	data.count<-0
	true.controls<-NA
	true.treats<-NA
	true.controls.SD<-NA
	true.treats.SD<-NA
	
	# Simulate k datasets
	while(data.count < k){
		
		# Get the sample size
		Study.nc<-ifelse(rbinom(1, 1, 0.5) == 0, control.n + rpois(1, var.n), control.n - rpois(1, var.n))
		Study.ne<-ifelse(rbinom(1, 1, 0.5) == 0, treatment.n + rpois(1, var.n), treatment.n - rpois(1, var.n))
		
		# must be at least 2
		Study.nc<-max(2, Study.nc)
		Study.ne<-max(2, Study.ne)
		
		# Find the study means
		study.control<-control.mean + rnorm(1, 0, b_study.c.M)
		study.treat<-treatment.mean + rnorm(1, 0, b_study.t.M)
		
		true.controls<-c(true.controls, study.control)
		true.treats<-c(true.treats, study.treat)
		
		study.control.SD<-abs(control.SD + rnorm(1, 0, b_study.c.SD))
		study.treatment.SD<-abs(treatment.SD + rnorm(1, 0, b_study.t.SD))
		
		true.controls.SD<-c(true.controls.SD, study.control.SD)
		true.treats.SD<-c(true.treats.SD, study.treatment.SD)
		
		# Simulate the data
		control.sample<-rnorm(Study.nc, study.control, study.control.SD)
		treatment.sample<-rnorm(Study.ne,  study.treat, study.treatment.SD)
		
		p.val<-t.test(control.sample, treatment.sample)$p.value
		tests<-c(tests, p.val)
				
		means.c<-c(means.c, mean(control.sample))
		SD.c<-c(SD.c, sd(control.sample))
		n.c<-c(n.c, Study.nc)

		means.e<-c(means.e, mean(treatment.sample))
		SD.e<-c(SD.e, sd(treatment.sample))
		n.e<-c(n.e, Study.ne)
			
		# Calculate d and SE.d
		d<-Calc.d(CMean = mean(control.sample), CSD = sd(control.sample), CN = length(control.sample), EMean = mean(treatment.sample), ESD = sd(treatment.sample), EN = length(treatment.sample))
		SE.d<-Calc.SE.d(CN = length(control.sample), EN = length(treatment.sample), d=d)

		ds<-c(ds, d)
		SE.ds<-c(SE.ds, SE.d)

		
		ifelse(PB == TRUE, data.count<-data.count + as.numeric(p.val < 0.05), data.count<-data.count+1)
	}
	
	# Remove NAs
	ds<-ds[-1]
	SE.ds<-SE.ds[-1]
	means.c<-means.c[-1]
	means.e<-means.e[-1]
	SD.c<-SD.c[-1]
	SD.e<-SD.e[-1]
	n.c<-n.c[-1]
	n.e<-n.e[-1]
	tests<-tests[-1]
	
	true.controls<-true.controls[-1]
	true.treats<-true.treats[-1]
	true.controls.SD<-true.controls.SD[-1]
	true.treats.SD<-true.treats.SD[-1]

	# Create dataframe and return it
	output<-data.frame(C.Mean = means.c, E.Mean = means.e)
	
	if(return.true.effects == TRUE){
		output$True.C.Mean<-true.controls
		output$True.E.Mean<-true.treats 
	}
	
	output$C.SD<-SD.c
	output$E.SD<-SD.e
		
	if(return.true.effects == TRUE){
		output$True.C.SD<-true.controls.SD
		output$True.E.SD<-true.treats.SD 
	}
	
	output$C.n<-n.c
	output$E.n<-n.e
	
	if(return.d == TRUE){
		output$d<-ds
		output$SE.d<-SE.ds
	}
	
	if (PB == TRUE){
		output$p<-tests
		output$Significant<-output$p < 0.05
	}
	
	return(output)
	
}

