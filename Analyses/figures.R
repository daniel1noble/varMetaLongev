#----------------------------------------------------------#
# Figure script
# Authors: D Noble and A. Senior
# Date: 28.09.16
#----------------------------------------------------------#

# Clear work space
rm(list = ls())

# Load the list of model objects
datObjects <- readRDS("./data/datObjects")

###################################
# Figure 1 - Mean-variance + Funnel plots
#-------------------------------------------------#
pdf(height = 14, width = 14, file = "./output/figures/figure1.pdf")
    par(mfrow=c(2, 2), mar=c(4,8,2.5,2.5))
	    #A)
	    plot(datObjects$data.long$lnMean, datObjects$data.long$lnSD, col=datObjects$data.long$Trt, ylab = "lnSD", xlab = "lnMean", las = 1, cex.lab = 1.5, cex.axis = 1.2, cex = 1.8)
	    #legend
	    points(x = c(2.5,2.5), y = c(5.8, 6.2), col = c("red", "black"), cex = 1.8)
	    text(c("Control", "Diet Restriction"), x = c(2.8,2.8), y = c(5.78, 6.15), adj = c(0,0), cex = 1.2)
	    mtext(side = 2, "A)", adj = 3, las = 1, padj = -10.5, font = 2, cex = 2)

	    par(bty = "n", mar = c(5,12,3,1))
		cex = 1.8
	    #lnVr
		plot(~lnVR, type = "n", data = datObjects$data, ylim = c(0,10), bty = "n", xlab = "lnVR", cex.lab = 1.5, cex.axis = 1.2, xlim = c(-0.2, 0.6))
		abline(v = 0, lty = 2)

		# Pre./post manipulation
		points(c(1,2) ~ PrePost[,1], pch = 16, cex = cex)
		arrows(x0=PrePost[,1] , y0= c(1,2), x1= PrePost[,2] , y1 = c(1,2), length = 0, angle = 90)
		arrows(x0=PrePost[,1] , y0= c(1,2), x1= PrePost[,3] , y1 = c(1,2), length = 0, angle = 90)
		mtext(side  = 2, row.names(PrePost), at = c(1,2), las = 1)
		
		#Adult diet
		points(c(4,5) ~ AdultDiet[,1], pch = 16, cex = cex)
		arrows(x0=AdultDiet[,1] , y0= c(4,5), x1= AdultDiet[,2] , y1 = c(4,5), length = 0, angle = 90)
		arrows(x0=AdultDiet[,1] , y0= c(4,5), x1= AdultDiet[,3] , y1 = c(4,5), length = 0, angle = 90)
		mtext(side  = 2, row.names(AdultDiet), at = c(4,5), las = 1)
		
		#Phylum
		points(c(7,8) ~ Phylum[,1], pch = 16, cex = cex)
		arrows(x0=Phylum[,1] , y0= c(7,8), x1= Phylum[,2] , y1 = c(7,8), length = 0, angle = 90)
		arrows(x0=Phylum[,1] , y0= c(7,8), x1= Phylum[,3] , y1 = c(7,8), length = 0, angle = 90)
		mtext(side  = 2, row.names(Phylum), at = c(7,8), las = 1)
		
		# n
		text(nVals[1:7], x = 0.5, y = c(0,1,2,4,5,7,8), cex = 1.5)
		text("N", x = 0.5, y = 9, font = 2, cex = 1.5)

		# Overall
		load("./output/models/lnVR.1.Rdata")
		mode <- as.numeric(mean(model$Sol[,"(Intercept)"]))
		HPD <- HPDinterval(model$Sol[,"(Intercept)"])

		points(y = 0, x = mode, pch = 18, col = "black", cex = cex + 1.5)
		arrows(x0=mode , y0= 0, x1= HPD[1] , y1 = 0, length = 0, angle = 90)
		arrows(x0=mode , y0= 0, x1= HPD[2] , y1 = 0, length = 0, angle = 90)
		
		#labels
		mtext(side  = 2, c("Overall", "Early Diet", "Adult Diet", "Phylum"), at = c(0, 2.8, 5.8, 8.8), las = 1, font = 2)
		mtext(side = 2, "B)", adj = 4, las = 1, padj = -10.5, font = 2, cex = 2)
	
	#ln CVR
		plot(~lnCVR, type = "n", data = datObjects$data, ylim = c(0,26), bty = "n", xlab = "lnCVR", cex.lab = 1.5, cex.axis = 1.2, xlim = c(-0.2, 0.7))
			abline(v = 0, lty = 2)

			# Pre./post manipulation
			points(c(2,3) ~ PrePostCVR[,1], pch = 16, cex = cex)
			arrows(x0=PrePostCVR[,1] , y0= c(2,3), x1= PrePostCVR[,2] , y1 = c(2,3), length = 0, angle = 90)
			arrows(x0=PrePostCVR[,1] , y0= c(2,3), x1= PrePostCVR[,3] , y1 = c(2,3), length = 0, angle = 90)
			mtext(side  = 2, row.names(PrePostCVR), at = c(2,3), las = 1)
			
			#Adult diet
			points(c(6,7) ~ AdultDietCVR[,1], pch = 16, cex = cex)
			arrows(x0=AdultDietCVR[,1] , y0= c(6,7), x1= AdultDietCVR[,2] , y1 = c(6,7), length = 0, angle = 90)
			arrows(x0=AdultDietCVR[,1] , y0= c(6,7), x1= AdultDietCVR[,3] , y1 = c(6,7), length = 0, angle = 90)
			mtext(side  = 2, row.names(AdultDietCVR), at = c(6,7), las = 1)
			
			#Phylum
			points(c(10,11) ~ PhylumCVR[,1], pch = 16, cex = cex)
			arrows(x0=PhylumCVR[,1] , y0= c(10,11), x1= PhylumCVR[,2] , y1 = c(10,11), length = 0, angle = 90)
			arrows(x0=PhylumCVR[,1] , y0= c(10,11), x1= PhylumCVR[,3] , y1 = c(10,11), length = 0, angle = 90)
			mtext(side  = 2, row.names(PhylumCVR), at = c(10,11), las = 1)
			
			#ManipCVR
			points(c(14,15) ~ ManipCVR[,1], pch = 16, cex = cex)
			arrows(x0=ManipCVR[,1] , y0= c(14,15), x1= ManipCVR[,2] , y1 = c(14,15), length = 0, angle = 90)
			arrows(x0=ManipCVR[,1] , y0= c(14,15), x1= ManipCVR[,3] , y1 = c(14,15), length = 0, angle = 90)
			mtext(side  = 2, row.names(ManipCVR), at = c(14,15), las = 1)

			#Sex
			points(c(18,19,20) ~ SexCVR[,1], pch = 16, cex = cex)
			arrows(x0=SexCVR[,1] , y0= c(18,19,20), x1= SexCVR[,2] , y1 = c(18,19,20), length = 0, angle = 90)
			arrows(x0=SexCVR[,1] , y0= c(18,19,20), x1= SexCVR[,3] , y1 = c(18,19,20), length = 0, angle = 90)
			mtext(side  = 2, row.names(SexCVR), at = c(18,19,20), las = 1)

			#Catchup
			points(c(23,24,25) ~ CatchCVR[,1], pch = 16, cex = cex)
			arrows(x0=CatchCVR[,1] , y0= c(23,24,25), x1= CatchCVR[,2] , y1 = c(23,24,25), length = 0, angle = 90)
			arrows(x0=CatchCVR[,1] , y0= c(23,24,25), x1= CatchCVR[,3] , y1 = c(23,24,25), length = 0, angle = 90)
			mtext(side  = 2, row.names(CatchCVR), at = c(23,24,25), las = 1)

			# Overall
			load("./output/models/lnCVR.1.Rdata")
			mode <- as.numeric(mean(model$Sol[,"(Intercept)"]))
			HPD <- HPDinterval(model$Sol[,"(Intercept)"])

			points(y = 0, x = mode, pch = 18, col = "black", cex = cex + 1.5)
			arrows(x0=mode , y0= 0, x1= HPD[1] , y1 = 0, length = 0, angle = 90)
			arrows(x0=mode , y0= 0, x1= HPD[2] , y1 = 0, length = 0, angle = 90)
			
			#labels
			mtext(side  = 2, c("Overall", "Early Diet", "Adult Diet", "Phylum", "Manipulation Type", "Sex", "Catchup Growth"), at = c(0, 4, 8, 12, 16, 21, 26)+0.2, las = 1, font = 2)
			mtext(side = 2, "C)", adj = 5, las = 1, padj = -10.5, font = 2, cex = 2)

			#N
			text(nVals, x = 0.55, y = c(0,2,3,6,7,10,11,14,15,18,19,20,23,24,25), cex = 1.5)
			text("N", x = 0.55, y = 26, font = 2, cex = 1.5)

	#ln RR
		plot(~lnRR, type = "n", data = datObjects$data, ylim = c(0,26), bty = "n", xlab = "lnRR", cex.lab = 1.5, cex.axis = 1.2, xlim = c(-0.4, 0.3))
			abline(v = 0, lty = 2)

			# Pre./post manipulation
			points(c(2,3) ~ PrePostsolPostlnRR[,1], pch = 16, cex = cex)
			arrows(x0=PrePostsolPostlnRR[,1] , y0= c(2,3), x1= PrePostsolPostlnRR[,2] , y1 = c(2,3), length = 0, angle = 90)
			arrows(x0=PrePostsolPostlnRR[,1] , y0= c(2,3), x1= PrePostsolPostlnRR[,3] , y1 = c(2,3), length = 0, angle = 90)
			mtext(side  = 2, row.names(PrePostsolPostlnRR), at = c(2,3), las = 1)
			
			#Adult diet
			points(c(6,7) ~ AdultDietsolPostlnRR[,1], pch = 16, cex = cex)
			arrows(x0=AdultDietsolPostlnRR[,1] , y0= c(6,7), x1= AdultDietsolPostlnRR[,2] , y1 = c(6,7), length = 0, angle = 90)
			arrows(x0=AdultDietsolPostlnRR[,1] , y0= c(6,7), x1= AdultDietsolPostlnRR[,3] , y1 = c(6,7), length = 0, angle = 90)
			mtext(side  = 2, row.names(AdultDietsolPostlnRR), at = c(6,7), las = 1)
			
			#Phylum
			points(c(10,11) ~ PhylumsolPostlnRR[,1], pch = 16, cex = cex)
			arrows(x0=PhylumsolPostlnRR[,1] , y0= c(10,11), x1= PhylumsolPostlnRR[,2] , y1 = c(10,11), length = 0, angle = 90)
			arrows(x0=PhylumsolPostlnRR[,1] , y0= c(10,11), x1= PhylumsolPostlnRR[,3] , y1 = c(10,11), length = 0, angle = 90)
			mtext(side  = 2, row.names(PhylumsolPostlnRR), at = c(10,11), las = 1)
			
			#ManipRR
			points(c(14,15) ~ ManipsolPostlnRR[,1], pch = 16, cex = cex)
			arrows(x0=ManipsolPostlnRR[,1] , y0= c(14,15), x1= ManipsolPostlnRR[,2] , y1 = c(14,15), length = 0, angle = 90)
			arrows(x0=ManipsolPostlnRR[,1] , y0= c(14,15), x1= ManipsolPostlnRR[,3] , y1 = c(14,15), length = 0, angle = 90)
			mtext(side  = 2, row.names(ManipsolPostlnRR), at = c(14,15), las = 1)

			#Sex
			points(c(18,19,20) ~ SexsolPostlnRR[,1], pch = 16, cex = cex)
			arrows(x0=SexsolPostlnRR[,1] , y0= c(18,19,20), x1= SexsolPostlnRR[,2] , y1 = c(18,19,20), length = 0, angle = 90)
			arrows(x0=SexsolPostlnRR[,1] , y0= c(18,19,20), x1= SexsolPostlnRR[,3] , y1 = c(18,19,20), length = 0, angle = 90)
			mtext(side  = 2, row.names(SexsolPostlnRR), at = c(18,19,20), las = 1)

			#Catchup
			points(c(23,24,25) ~ CatchsolPostlnRR[,1], pch = 16, cex = cex)
			arrows(x0=CatchsolPostlnRR[,1] , y0= c(23,24,25), x1= CatchsolPostlnRR[,2] , y1 = c(23,24,25), length = 0, angle = 90)
			arrows(x0=CatchsolPostlnRR[,1] , y0= c(23,24,25), x1= CatchsolPostlnRR[,3] , y1 = c(23,24,25), length = 0, angle = 90)
			mtext(side  = 2, row.names(CatchsolPostlnRR), at = c(23,24,25), las = 1)

			# Overall
			load("./output/models/lnRR.1.Rdata")
			mode <- as.numeric(mean(model$Sol[,"(Intercept)"]))
			HPD <- HPDinterval(model$Sol[,"(Intercept)"])

			points(y = 0, x = mode, pch = 18, col = "black", cex = cex + 1.5)
			arrows(x0=mode , y0= 0, x1= HPD[1] , y1 = 0, length = 0, angle = 90)
			arrows(x0=mode , y0= 0, x1= HPD[2] , y1 = 0, length = 0, angle = 90)
			
			#labels
			mtext(side  = 2, c("Overall", "Early Diet", "Adult Diet", "Phylum", "Manipulation Type", "Sex", "Catchup Growth"), at = c(0, 4, 8, 12, 16, 21, 26) + 0.2, las = 1, font = 2)
			mtext(side = 2, "D)", adj =5, las = 1, padj = -10.5, font = 2, cex = 2)

			#N
			text(nVals, x = 0.25, y = c(0,2,3,6,7,10,11,14,15,18,19,20,23,24,25), cex = 1.5)
			text("N", x = 0.25, y = 26, font = 2, cex = 1.5)

dev.off()	        

###################################
# Figure 2 - Forest plots
#--------------------------------------------------#

pdf(height = 6, width = 14.5, file = "./output/figures/figure2.pdf")
	par(mfrow = c(1,3), bty = "n", mar = c(5,8,3,1), cex.lab = 2)
		cex = 1.8
	    # C) Funnel plot lnVR

	    MLMAlnVR<-rma.mv(yi = lnVR, mod = ~AdultDiet + ExptLifeStage + ManipType + Sex + CatchUp + Phylum, V = datObjects$VlnVR, random=~1|StudyNo/Map, data=datObjects$data)
	    funnel(MLMAlnVR, yaxis = "seinv", digits = 1, las = 1, xlab = "lnVR", ylab = "Precision (1/s.e.)", cex.axis = 1.2, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"))
	     mtext("A)", adj = -0.25, font = 2)

	     # B) Funnel plot lnCVR

	    MLMAlnCVR<-rma.mv(yi = lnCVR, mods= ~AdultDiet + ExptLifeStage + ManipType + Sex + CatchUp + Phylum, V = datObjects$VlnCVR, random=~1|StudyNo/Map, data=datObjects$data)
	    funnel(MLMAlnCVR, yaxis = "seinv", digits = 1, las = 1, xlab = "lnCVR", ylab = "", cex.axis = 1.2, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"))
	     mtext("B)", adj = -0.30, font = 2)	

	# A) Funnel plot lnRR

	    MLMAlnRR<-rma.mv(yi = lnRR, mods = ~AdultDiet + ExptLifeStage + ManipType + Sex + CatchUp + Phylum, V = datObjects$VlnRR, random=~1|StudyNo/Map, data=datObjects$data)
	    funnel(MLMAlnRR, yaxis = "seinv", digits = 1, las = 1, xlab = "lnRR", ylab = "", cex.axis = 1.2, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"))
	     mtext("C)", adj = -0.30, font = 2)	    

dev.off()
			
			


		
			


