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
pdf(height = 6.5, width = 6.5, file = "./output/figures/figure1.pdf")
    par(mfrow=c(2, 2), mar=c(4,4,1,1))
	    #A)
	    plot(datObjects$data.long$lnMean, datObjects$data.long$lnSD, col=datObjects$data.long$Trt, ylab = "lnSD", xlab = "lnMean", las = 1)
	    #legend
	    points(x = c(2.5,2.5), y = c(5.8, 6.2), col = c("red", "black"))
	    text(c("Control", "Diet Restriction"), x = c(2.8,2.8), y = c(5.7, 6.1), adj = c(0,0))
	    mtext("A)", adj = -0.20, font = 2)

	   # B) Funnel plot lnRR

	    MLMAlnRR<-rma.mv(yi = lnRR, V = datObjects$VlnRR, random=~1|StudyNo/Map, data=datObjects$data)
	    funnel(MLMAlnRR, yaxis = "seinv", digits = 1, las = 1, xlab = "lnRR", ylab = "Precision (1/s.e.)", cex.axis = 0.8)
	     mtext("B)", adj = -0.30, font = 2)	    

	     # C) Funnel plot lnCVR

	    MLMAlnCVR<-rma.mv(yi = lnCVR, V = datObjects$VlnCVR, random=~1|StudyNo/Map, data=datObjects$data)
	    funnel(MLMAlnCVR, yaxis = "seinv", digits = 1, las = 1, xlab = "lnCVR", ylab = "Precision (1/s.e.)", cex.axis = 0.8)
	     mtext("C)", adj = -0.20, font = 2)	

	      # C) Funnel plot lnCVR

	    MLMAlnVR<-rma.mv(yi = lnVR, V = datObjects$VlnVR, random=~1|StudyNo/Map, data=datObjects$data)
	    funnel(MLMAlnVR, yaxis = "seinv", digits = 1, las = 1, xlab = "lnVR", ylab = "Precision (1/s.e.)", cex.axis = 0.8)
	     mtext("D)", adj = -0.30, font = 2)
dev.off()	        

###################################
# Figure 2 - Forest plots
#--------------------------------------------------#

