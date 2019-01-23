####################
# Fst File to Plot #
####################

# Load modules if installed,
# if not installed, use 'install.packages("package_name")'
library("Cairo")
library("stringr")


# *** Specify the path to Fst File ***
# Remember that the path must have forward slashes, not backwards
path_to_snps <- "C:/Users/miles_000/Downloads/Research/Tooth/results/snps.fst"
path_to_indels <- "C:/Users/miles_000/Downloads/Research/Tooth/results/indels.fst"

# Load the data
snps_data <- read.table(path_to_snps, sep="\t", header=TRUE)
indels_data <- read.table(path_to_indels, sep="\t", header=TRUE)
# Convert the negative fst scores to 0
snps_data[, c(5)][snps_data[, c(5)] < 0] <- 0
indels_data[, c(5)][indels_data[, c(5)] < 0] <- 0
# Make the indels be plotted below the SNPs by making all their Zfst score inverted
# Add a column that just has row number
snps_data$ID <- seq.int(nrow(snps_data))
indels_data$ID <- rep(0, 93664)

indels_data <- readRDS(file = "indels_data.rds")
# Only to be run first time, save the data, then load it
# ~20 minute runtime
# for (i in 1:nrow(indels_data)) {
#   indel_start = indels_data[i,2]
#   for (j in i:nrow(snps_data)) {
#     snps_start = snps_data[j,2]
#     snps_id = snps_data[j,7]
#     if (indel_start == snps_start) {
#       indels_data[i,7] <- snps_id
#       break
#     }
#   }
# }
#indels_data <- indels_data[indels_data$ID != 0,]
#saveRDS(indels_data, file = "indels_data.rds")
# Add a column for Zfst
snps_data$Zfst <- (( snps_data$WEIGHTED_FST - mean(snps_data$WEIGHTED_FST) ) / sd(snps_data$WEIGHTED_FST)) + 1
indels_data$Zfst <- -(( indels_data$WEIGHTED_FST - mean(indels_data$WEIGHTED_FST) ) / sd(indels_data$WEIGHTED_FST)+1)

total <- rbind(snps_data, indels_data)

total$LG <- sub("^NW.*", "36802.1", total$CHROM)
total$LG <- sub("^NC_027944.1", "36802.1", total$LG)
total$LG <- sub("^NC_0", "", total$LG)

# Plot the data
image_name <- "UMD2a_Bicuspid_vs_Tricuspid_Zfst_Separated.png"
#svg("UMD2a_Bicuspid_vs_Tricuspid_Zfst_Separated.svg", width = 1350, height = 401)
Cairo(image_name, type="png", , width = 12, height = 4, units = 'in', res=300)
palette("default")
plot(total$ID, total$Zfst, col=as.numeric(as.character(total$LG)), xlab="Linkage Group", ylab=expression("Z"[fst]), xaxt="n", yaxt="n", bty="n", main="", pch=20, cex=0.70)
abline(h = 0)
title(main = expression("UMD2a Bicuspid vs Tricuspid Z"[fst]))

text(2, 15, labels="SNPs", cex=0.8)
text(2, -15, labels="Indels", cex=0.8)

yaxisat <- c(-20, 0, 20)
yaxislabels <- c("high", "low", "high")
axis(2, at=yaxisat, labels=yaxislabels)

xaxisat <- c(0, 3864, 7122, 10841, 13887, 17501, 21472, 27958, 30354, 32449, 35681, 38921, 42326, 45524, 49300, 52739, 56205, 59777, 62721, 65312, 68282, 71749, 75949)
xaxislabels <- as.character(1:20)
xaxislabels <- c(xaxislabels, c("22", "23", "unplaced"))
axis(1, at=xaxisat,labels=xaxislabels, cex.axis=0.75)
#dev.copy(png,"/mnt/c/Users/miles_000/Downloads/Research/Tooth/scripts/fst/BHLHE40_Bicuspid_vs_Tricuspid_Fst.png",width=8, height=10, units="in",res=750)
dev.off()
