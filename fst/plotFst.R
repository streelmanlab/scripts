####################
# Fst File to Plot #
####################

# Load modules if installed,
# if not installed, use 'install.packages("package_name")'
library(swirl)
library("stringr")
library(ggplot2)


# *** Specify the path to Fst File ***
# Remember that the path must have forward slashes, not backwards
path_to_fst <- "C:/Users/miles_000/Downloads/Research/Tooth/results/fst.vcf.windowed.weir.fst"

# Load the data
fst_data <- read.table(path_to_fst, sep="\t", header=TRUE)
# Convert the negative fst scores to 0
fst_data[, c(5)][fst_data[, c(5)] < 0] <- 0
# Add a column that just has row number
fst_data$ID <- seq.int(nrow(fst_data))
# Add a column that for LG numbers
fst_data$LG <- sub("^NW.*", "36802.1", fst_data$CHROM)
fst_data$LG <- sub("^NC_027944.1", "36802.1", fst_data$LG)
fst_data$LG <- sub("^NC_0", "", fst_data$LG)
#fst_data$ID2 <- as.numeric(fst_data$LG)
#fst_data$ID2 <- fst_data$ID2 - 36779.1
#fst_data$LG <- as.numeric(fst_data$LG)
# Add a column for Zfst
fst_data$Zfst <- ( fst_data$WEIGHTED_FST - mean(fst_data$WEIGHTED_FST) ) / sd(fst_data$WEIGHTED_FST)
fst_data <- fst_data[,c("ID","Zfst", "LG", "CHROM", "BIN_START","BIN_END", "WEIGHTED_FST")]




# Plot the data
jpeg("UMD2a_Bicuspid_vs_Tricuspid_Zfst_edit.jpg", width = 1350, height = 401)
plot(fst_data$ID, fst_data$Zfst, col=as.numeric(as.character(fst_data$LG)), xlab="Linkage Group", ylab=expression("Z"[fst]), xaxt="n", main="", pch=20 )
title(main = expression("UMD2a Bicuspid vs Tricuspid Z"[fst]))

#xaxisat <- c(3864, 7122, 10841, 13887, 17501, 21472, 27958, 30354, 32449, 35681, 38921, 42326, 45524, 49300, 52739, 56205, 59777, 62721, 65312, 68282, 71749, 75949, 95000)
xaxisat <- c(0, 3864, 7122, 10841, 13887, 17501, 21472, 27958, 30354, 32449, 35681, 38921, 42326, 45524, 49300, 52739, 56205, 59777, 62721, 65312, 68282, 71749, 75949)
xaxislabels <- as.character(1:20)
xaxislabels <- c(xaxislabels, c("22", "23", "unplaced"))
axis(1, at=xaxisat,labels=xaxislabels)
dev.off()

