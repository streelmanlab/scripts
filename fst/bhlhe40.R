####################
# Fst File to Plot #
####################

# Load modules if installed,
# if not installed, use 'install.packages("package_name")'
#library(swirl)
#library("stringr")
#library(ggplot2)


# *** Specify the path to Fst File ***
# Remember that the path must have forward slashes, not backwards
#path_to_fst <- "C:/Users/miles_000/Downloads/Research/Tooth/results/bhlhe40_region.fst"
path_to_fst <- "/mnt/c/Users/miles_000/Downloads/Research/Tooth/results/bhlhe40_region.fst"

# Load the data
fst_data <- read.table(path_to_fst, sep="\t", header=FALSE)
# Rename the columns
names(fst_data) <- c("LG", "POS", "FST")
# Convert the negative fst scores to 0
fst_data[, c(3)][fst_data[, c(3)] < 0] <- 0

# Custom Coloring Scheme
palette(c("grey68", "grey50", "grey0"))
fst_data$COL <- fst_data$FST
fst_data$COL[fst_data$FST > .99] <- 3
fst_data$COL[fst_data$FST < .99] <- 2
fst_data$COL[fst_data$FST < .50] <- 1

# Plot the data
png_name <- "BHLHE40_Bicuspid_vs_Tricuspid_Fst.png"
png(png_name, width = 1350, height = 401)
par(fig=c(0,1,0,0.8))
plot(fst_data$POS, fst_data$FST, col=fst_data$COL, xlab="Position on LG5", ylab=expression("Fst"), main="", pch=20)
#Plot the genes
par(fig=c(0,1,.5,1), new=TRUE)
plot(seq(7127936, 7650068, 1), rep(1,7650069-7127936) , col="white", main="", xlab="", ylab="", yaxt="n", xaxt="n", bty="n")
title(main = expression("BHLHE40 Bicuspid vs Tricuspid Fst"))
coord <- list(c(7127935, 7139474, "LOC101474366"), c(7167098,  7170405, "LOC101474077"), c(7200706, 7235571, "mitfa"), c(7247226, 7289133, "frmd4b"), c(7293110, 7297473, "ARL6IP5"), c(7297541, 7306936, "trnt1"), c(7307116, 7310492, "AVPR2"), c(7405790, 7419704, "lrrn1"), c(7432720, 7438156, "LOC101470671"), c(7445848, 7526086, "ITPR1"), c(7528132, 7531739, "bhlhe40"), c( 7543994, 7547291, "LOC112434879"), c (7549522, 7551078, "LOC112434880"), c(7555195, 7558074, "LOC112434795"), c(7567775, 7570522, "LOC112434846"), c(7574807, 7576428, "LOC112434841"), c(7586405, 7589692, "LOC112434866"), c(7591989, 7620486, "LOC101480065"), c(7629619, 7635457, "NTN4"), c(7635595, 7638445, "CHCHD4"), c(7645665, 7650070, "SEC61A1"))

xstart = 0
ystart = 0
ystop = 0.8
for (pair in coord) {
  genestart <- pair[1]
  genestop  <- pair[2]
  name      <- pair[3]
  x <- c(genestart, genestart, genestop, genestop)
  y <- c(ystart, ystop, ystop, ystart)
  polygon(x, y, col="green", border=NA)
  text(as.numeric(genestart) + (as.numeric(genestop)-as.numeric(genestart))/2,as.numeric(0.93), labels = name, srt=45, cex=0.8)
}
#polygon(c(7127936, 7127936, 7650068, 7650068), c(0,1,1,0), col="green")
dev.copy(png,"/mnt/c/Users/miles_000/Downloads/Research/Tooth/scripts/fst/BHLHE40_Bicuspid_vs_Tricuspid_Fst.png",width=8, height=10, units="in",res=750)
#dev.copy(jpeg,jpeg_name,width=8, height=10, units="in",res=500) 
dev.off()
#par(fig=c(0,1,.5,1), new=TRUE)
#barplot(genes$present, col="green", axes=FALSE, axisnames=FALSE, border=NA, space =NULL)
#barplot(genes$present, col="green", names=gene_names, bty="n", xaxt="n", yaxt="n", ylab="", xlab="")

#xaxisat <- c(3864, 7122, 10841, 13887, 17501, 21472, 27958, 30354, 32449, 35681, 38921, 42326, 45524, 49300, 52739, 56205, 59777, 62721, 65312, 68282, 71749, 75949, 95000)
#xaxisat <- c(0, 3864, 7122, 10841, 13887, 17501, 21472, 27958, 30354, 32449, 35681, 38921, 42326, 45524, 49300, 52739, 56205, 59777, 62721, 65312, 68282, 71749, 75949)
#xaxislabels <- as.character(1:20)
#xaxislabels <- c(xaxislabels, c("22", "23", "unplaced"))
#axis(1, at=xaxisat,labels=xaxislabels)


# Genes
#pos <- seq(7127936, 7650068, 1)
#present <- rep(0, 7650069 - 7127936)
#genes = data.frame(pos, present)
# Upstream genes
#genes$present[genes$pos > 7127935 & genes$pos < 7139474 ] <- 1 
#genes$present[genes$pos > 7167098 & genes$pos < 7170405 ] <- 1 
#genes$present[genes$pos > 7200706 & genes$pos < 7235571 ] <- 1 
#genes$present[genes$pos > 7247226 & genes$pos < 7289133 ] <- 1 
#genes$present[genes$pos > 7293110 & genes$pos < 7297473 ] <- 1 
#genes$present[genes$pos > 7297541 & genes$pos < 7306936 ] <- 1 
#genes$present[genes$pos > 7307116 & genes$pos < 7310492 ] <- 1 
#genes$present[genes$pos > 7405790 & genes$pos < 7419704 ] <- 1 
#genes$present[genes$pos > 7432720 & genes$pos < 7438156 ] <- 1 
#genes$present[genes$pos > 7445848 & genes$pos < 7526086 ] <- 1
# Downstream Genes
#genes$present[genes$pos > 7543994 & genes$pos < 7547291 ] <- 1
#genes$present[genes$pos > 7549522 & genes$pos < 7551078 ] <- 1
#genes$present[genes$pos > 7555195 & genes$pos < 7558074 ] <- 1
#genes$present[genes$pos > 7567775 & genes$pos < 7570522 ] <- 1
#genes$present[genes$pos > 7574807 & genes$pos < 7576428 ] <- 1
#genes$present[genes$pos > 7586405 & genes$pos < 7589692 ] <- 1
#genes$present[genes$pos > 7591989 & genes$pos < 7620486 ] <- 1
#genes$present[genes$pos > 7629619 & genes$pos < 7635457 ] <- 1
#genes$present[genes$pos > 7635595 & genes$pos < 7638445 ] <- 1
#genes$present[genes$pos > 7645665 & genes$pos < 7650070 ] <- 1
#gene_names <- c("LOC101474366", "LOC101474077", "mitfa", "frmd4b", "LOC101471523", "trnt1", "LOC101471235", "lrrn1", "LOC101470671", "LOC101470371", "LOC112434879", "LOC112434880", "LOC112434795", "LOC112434846", "LOC112434841", "LOC112434866", "LOC101480065", "LOC101479783", "LOC101467661", "LOC101467346")
