setwd("C:/Users/Jonathan/Documents/UCLA/Pelligrini/projects/geomx/")

library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(readxl)
library(openxlsx)
library(knitr)
library(ggplot2)

# DCC filenames 
DCCFiles <- c(dir("current_data/geomx_files/DCC_batch_4", pattern = ".dcc$",full.names = TRUE, recursive = TRUE),
              dir("current_data/geomx_files/DCC_batch_5", pattern = ".dcc$",full.names = TRUE, recursive = TRUE),
              dir("current_data/geomx_files/DCC_batch_3", pattern = ".dcc$",full.names = TRUE, recursive = TRUE),
              dir("current_data/geomx_files/DCC_batch_1", pattern = ".dcc$",full.names = TRUE, recursive = TRUE),
              dir("current_data/geomx_files/DCC_batch_2", pattern = ".dcc$",full.names = TRUE, recursive = TRUE))

# Configuration filename
PKCFiles <- "current_data/geomx_files/Hs_R_NGS_WTA_v1.0.pkc"

# Annotation filename
SampleAnnotationFile <- "current_data/geomx_files/anno/geomx_annotations_all_batches.xlsx"

# Load data 
batched_data <- readNanoStringGeoMxSet(dccFiles = DCCFiles,
                                       pkcFiles = PKCFiles,
                                       phenoDataFile = SampleAnnotationFile,
                                       phenoDataSheet = "Sheet1",
                                       phenoDataDccColName = "Sample_ID",
                                       protocolDataColNames = c("ROI","Segment label"))

# Offset and flag 
batched_data <- shiftCountsOne(batched_data, useDALogic = TRUE)

QC_params <- list(minSegmentReads = 1000, 
       percentTrimmed = 80,    
       percentStitched = 80,  
       percentAligned = 75,    
       percentSaturation = 50, # Minimum sequencing saturation 
       minNegativeCount = 1,   # Minimum negative control counts 
       maxNTCCount = 9000,     # Maximum counts observed in NTC well 
       minNuclei = 20,        
       minArea = 1000)        

batched_data <- setSegmentQCFlags(batched_data, 
                                  qcCutoffs = QC_params)  

# Summarize QC results
QCResults <- protocolData(batched_data)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults$QCStatus <- apply(QCResults, 1, function(x) ifelse(sum(x) == 0, "PASS", "WARNING"))
QC_Summary["TOTAL FLAGS", ] <- c(sum(QCResults[, "QCStatus"] == "PASS"),
                                 sum(QCResults[, "QCStatus"] == "WARNING"))

# Calculate negative geometric means
negativeGeoMeans <- esBy(negativeControlSubset(batched_data), 
       GROUP = "Module", 
       FUN = function(x) { 
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
       }) 
protocolData(batched_data)[["NegGeoMean"]] <- negativeGeoMeans

pkcs <- annotation(batched_data)
modules <- gsub(".pkc", "", pkcs)
negCols <- paste0("NegGeoMean_", modules)
pData(batched_data)[, negCols] <- sData(batched_data)[["NegGeoMean"]]

# Detatch negative geometric mean column ahead of aggregateCounts call
pData(batched_data) <- pData(batched_data)[, !colnames(pData(batched_data)) %in% negCols]
batched_data <- batched_data[, QCResults$QCStatus == "PASS"]


# Set BioProbe QC flags based on probe ratio and Grubb’s test outlier detection
batched_data <- setBioProbeQCFlags(batched_data, 
                                   qcCutoffs = list(minProbeRatio = 0.1,
                                                    percentFailGrubbs = 20), 
                                   removeLocalOutliers = TRUE)

# Define QC metrics for each probe based on outcomes of QC tests
ProbeQCResults <- fData(batched_data)[["QCFlags"]]
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))

# Subset data to include only probes that passed QC checks
ProbeQCPassed <- subset(batched_data,
                        fData(batched_data)[["QCFlags"]][, "LowProbeRatio"] == FALSE &
                          fData(batched_data)[["QCFlags"]][, "GlobalGrubbsOutlier"] == FALSE)


batched_data <- ProbeQCPassed 

# Aggregate counts for the final data set
target_batched_data <- aggregateCounts(batched_data)
pData(target_batched_data)$roi <- protocolData(target_batched_data)[["ROI"]]

# filter out data sets we have ignored and DCC control files 
index <- !grepl("TBL 7", pData(target_batched_data)$`Scan name`) & !grepl("TBL 60", pData(target_batched_data)$`Scan name`) & 
  !grepl("No Template Control", pData(target_batched_data)$`Segment tags`)


# Separate out metadata and counts 
batch_meta<-cbind(pData(target_batched_data)[index,],sData(target_batched_data)[index,])
batch_count<-exprs(target_batched_data)[,index]

# Separate out negative probes 
negativeProbefData <- subset(fData(batched_data), CodeClass == "Negative")
batch_ngpb <- exprs(batched_data[fData(batched_data)$TargetName %in% negativeProbefData$TargetName,])
batch_ngpb <- batch_ngpb[,index]

# Comprehensive row names: region, sample, ROI
batch_names <- pData(target_batched_data)[colnames(batch_count)[grepl("dcc",colnames(batch_count))],c("Segment tags","Scan name","roi")]
batch_names$roi <- sprintf("%03d",batch_names$roi)
batch_names$`Scan name` <- gsub("PTB","PTB ",batch_names$`Scan name`)
batch_names <- apply(batch_names,1,function(x) paste(x, collapse =" | "))

# update names 
colnames(batch_count) <- batch_names[colnames(batch_count)]

# Data structure for nuclei counts 
nuclei <- rbind(batch_meta[c("AOINucleiCount","AOISurfaceArea","ROICoordinateX","ROICoordinateY","Aligned")])
rownames(nuclei) <- c(batch_names[rownames(batch_meta["AOINucleiCount"])])
nuclei["AlignedPercent"] <- unname(unlist(batch_meta["Aligned (%)"]))

# Calculate LOQ 
pkcs <- annotation(batched_data)
modules <- gsub(".pkc", "", pkcs)
module <- modules[1]

# LOQ SD threshold and minimum value
cutoff <- 2
minLOQ <- 2

LOQ <- data.frame(row.names = colnames(batch_count))

NegGeoMean <- apply(batch_ngpb,2,ngeoMean)
NegGeoSD <- apply(batch_ngpb,2,ngeoSD)

LOQ <- data.frame(row.names = colnames(batch_count))
LOQ[, module] <-pmax(minLOQ, NegGeoMean * NegGeoSD ^ cutoff)

# Threshold counts against LOQ 
LOQ_Mat <- c()
Mat_i <- apply(batch_count,1,FUN = function(x) { x > LOQ})
LOQ_Mat <- t(rbind(LOQ_Mat, Mat_i))
colnames(LOQ_Mat) <- colnames(batch_count)

# Filter out spots by number of genes detected 
GenesDetected <- colSums(LOQ_Mat, na.rm = TRUE)
GeneDetectionRate <- GenesDetected / nrow(batch_count)
DetectionThreshold <- cut(GeneDetectionRate,
                          breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
                          labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))

batch_count <- batch_count[,GeneDetectionRate >= .05]
batch_ngpb <- batch_ngpb[,GeneDetectionRate >= .05]
nuclei <- nuclei[GeneDetectionRate >= .05,]

# Filter out genes by spot detection 
LOQ_Mat <- LOQ_Mat[, colnames(batch_count)]
DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
DetectionRate <- DetectedSegments / ncol(batch_count)

batch_count <- batch_count[DetectionRate >= 0.05,]

# Q3 normalization
qs <- apply(batch_count, 2, function(x) stats::quantile(x, 0.75))
q3 <- sweep(batch_count, 2L, qs / ngeoMean(qs), FUN = "/")

# Save to CSV
write.table(nuclei,"current_data/nuclei_2024b.csv",sep=",",col.names=NA)
write.table(batch_ngpb,"current_data/geomx_negprobes_2024b.csv",sep=",",col.names=NA)
write.table(batch_count,"current_data/geomx_integrated_raw_data_2024b.csv",sep=",",col.names=NA)
write.table(q3,"current_data/geomx_integrated_q3_data_2024b.csv",sep=",",col.names=NA)
