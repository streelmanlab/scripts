####################
# Fst File to Plot #
####################

# Load modules if installed,
# if not installed, use 'install.packages("package_name")'
#library(swirl)
#library("stringr")
#library(ggplot2)
library(Cairo)


# *** Specify the path to Fst File ***
# Remember that the path must have forward slashes, not backwards
path_to_indel <- "C:/Users/miles_000/Downloads/Research/Tooth/results/bhlhe40_indels.fst"

# Load the data
indel_fst_data <- read.table(path_to_indel, sep="\t", header=FALSE)
# Rename the columns
names(indel_fst_data) <- c("LG", "POS", "FST")
indel_fst_data <- indel_fst_data[as.numeric(as.character(indel_fst_data$POS)) > 7127934,]
indel_fst_data <- indel_fst_data[as.numeric(as.character(indel_fst_data$POS)) < 7650068,]

# Convert the negative fst scores to 0
indel_fst_data[, c(3)][indel_fst_data[, c(3)] < 0] <- 0

# Custom Coloring Scheme
palette(c("grey68", "grey50", "grey0"))
indel_fst_data$COL <- indel_fst_data$FST
indel_fst_data$COL[indel_fst_data$FST > .99] <- 3
indel_fst_data$COL[indel_fst_data$FST < .99] <- 2
indel_fst_data$COL[indel_fst_data$FST < .50] <- 1

# Plot the data
png_name <- "BHLHE40_Bicuspid_vs_Tricuspid_Fst_Indel.png"
Cairo(png_name, type="png", , width = 12, height = 4, units = 'in', res=300)
par(fig=c(0,1,0,0.8))
plot(indel_fst_data$POS, indel_fst_data$FST, col=indel_fst_data$COL, xlab="Position on LG5", ylab=expression("Fst"), xaxt="n", ylim=c(0,1), main="", pch=20)
xaxisat <- c(7127935, 7200000, 7300000, 7400000, 7500000, 7600000, 7650068)
axis(1, at=xaxisat, labels=xaxisat)
#Plot the genes
par(fig=c(0,1,.3,1), new=TRUE)
plot(seq(7127936, 7650068, 1), rep(1,7650069-7127936) , col="white", main="", xlab="", ylab="", yaxt="n", xaxt="n", bty="n")
title(main = expression("BHLHE40 Bicuspid vs Tricuspid Fst - Indels"), line = -0.5)


coord <- list(c(7127935, 7139474, "101474366"), c(7167098,  7170405, "101474077"), c(7200706, 7235571, "MITFA"), c(7247226, 7289133, "FRMD4B"), c(7293110, 7297473, "ARL6IP5"), c(7297541, 7306936, "TRNT1"), c(7307116, 7310492, "AVPR2"), c(7405790, 7419704, "LRRN1"), c(7432720, 7438156, "101470671"), c(7445848, 7526086, "ITPR1"), c(7528132, 7531739, "bhlhe40"), c( 7543994, 7547291, "lncRNA"), c (7549522, 7551078, "lncRNA"), c(7555195, 7558074, "lncRNA"), c(7567775, 7570522, "lncRNA"), c(7574807, 7576428, "lncRNA"), c(7586405, 7589692, "lncRNA"), c(7591989, 7620486, "101480065"), c(7629619, 7635457, "NTN4"), c(7635595, 7638445, "CHCHD4"), c(7645665, 7650070, "SEC61A1"))

xstart = 0
ystart = 0.8
ystop = 1
for (pair in coord) {
  genestart <- pair[1]
  genestop  <- pair[2]
  name      <- pair[3]
  x <- c(genestart, genestart, genestop, genestop)
  y <- c(ystart, ystop, ystop, ystart)
  polygon(x, y, col="green", border=NA)
  text(as.numeric(genestart) + (as.numeric(genestop)-as.numeric(genestart))/2,as.numeric(0.95), labels = name, srt=45, cex=0.8)
}
dev.off()