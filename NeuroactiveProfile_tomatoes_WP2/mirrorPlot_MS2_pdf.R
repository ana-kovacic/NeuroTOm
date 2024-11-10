library(MetaboCoreUtils)
library(MsCoreUtils)
library(Spectra)
library(readxl)

filterHigherPrecursorPeaks <- function(x, precursorMz, tolerance = 0.05, ...){
  pmz <- precursorMz + tolerance
  x[x[, "mz"] < pmz, , drop = FALSE]
}

cmps <- readxl::read_xlsx("data/AnnotatedComp_tomatoProfile_filtered.xlsx")
cmps$mass <- calculateMass(cmps$formula)
cmps$POS <- as.numeric(mass2mz(cmps$mass, "[M+H]+"))
cmps$NEG <- as.numeric(mass2mz(cmps$mass, "[M-H]-"))

excluded <- data.frame(
  rbind(
    c("4-Aminobenzoic acid", "NEG"),
    c("5'-Deoxy-5'-(methylthio)adenosine", "NEG"),
    c("Aspartame", "NEG"),
    c("Guanine", "NEG"),
    c("L-Glutamic acid", "NEG"),
    c("L-glutamine", "NEG"),
    c("L-Tyrosine", "NEG"),
    c("Tryptophan", "NEG"),
    c("5-caffeoylshikimic acid", "POS"),
    c("Citric acid", "POS"),
    c("Naringin", "NEG"),
    c("Resveratrol", "NEG"),
    c("Resveratrol", "POS"),
    c("Rutin", "POS"),
    c("Tryptamine", "NEG")
     )
) 
colnames(excluded) <- c("name", "POL")

excluded_std <- data.frame(
  rbind(
   c("Ascorbic Acid", "POS"),
    c("Caffeic acid", "POS"),
    c("Hesperetin", "POS"),
    c("L-Dopa", "POS"),
    c("Quercetin", "NEG")
  )
) 
colnames(excluded_std) <- c("name", "POL")

load("rdata/POS_combined_data_dda.RData")
standards <- combined_spectra[!grepl(paste(
  c(paste0("_", c("D1", "D2", "D3", "D4", "D5", "D6", "DO1", "DO2", "DO3", "DO4", "DO5", "DO6", 
                "DOP", "DP", "PF", "PFP", "PS", "PSP"), "_"), "Blank", "QC", "PPB"),collapse = "|"),
  combined_spectra$dataOrigin)]
standards <- addProcessing(standards, filterHigherPrecursorPeaks, 
                                  spectraVariables = c("precursorMz"))
idx <- which(unlist(lapply(mz(standards), length)) == 0)
if(length(idx) > 0){
  standards <- standards[-idx]
}
standards <- standards[order(standards$entropy)]
std_POS <- standards

combined_spectra <- combined_spectra[grepl(paste(paste0("_", c("D1", "D2", "D3", "D4", "D5", "D6", "DO1", "DO2", "DO3", "DO4", "DO5", "DO6", 
                                                               "DOP", "DP", "PF", "PFP", "PS", "PSP"), "_"),collapse = "|"),
                                           combined_spectra$dataOrigin)]
combined_spectra <- addProcessing(combined_spectra, filterHigherPrecursorPeaks, 
                                  spectraVariables = c("precursorMz"))
idx <- which(unlist(lapply(mz(combined_spectra), length)) == 0)
if(length(idx) > 0){
  combined_spectra <- combined_spectra[-idx]
}
combined_spectra <- combined_spectra[order(combined_spectra$entropy)]
ms2_POS <- combined_spectra

load("rdata/NEG_combined_data_dda.RData")
#unique(basename(combined_spectra$dataOrigin))
standards <- combined_spectra[!grepl(paste(
  c(paste0("_", c("D1", "D2", "D3", "D4", "D5", "D6", "DO1", "DO2", "DO3", "DO4", "DO5", "DO6", 
                  "DOP", "DP", "PF", "PFP", "PS", "PSP"), "_"), "Blank", "QC", "PPB"),collapse = "|"),
  combined_spectra$dataOrigin)]
standards <- addProcessing(standards, filterHigherPrecursorPeaks, 
                           spectraVariables = c("precursorMz"))
idx <- which(unlist(lapply(mz(standards), length)) == 0)
if(length(idx) > 0){
  standards <- standards[-idx]
}
standards <- standards[order(standards$entropy)]
std_NEG <- standards

combined_spectra <- combined_spectra[grepl(paste(paste0("_", c("D1", "D2", "D3", "D4", "D5", "D6", "DO1", "DO2", "DO3", "DO4", "DO5", "DO6", 
                                                               "DOP", "DP", "PF", "PFP", "PS", "PSP"), "_"),collapse = "|"),
                                           combined_spectra$dataOrigin)]
combined_spectra <- addProcessing(combined_spectra, filterHigherPrecursorPeaks, 
                                  spectraVariables = c("precursorMz"))
idx <- which(unlist(lapply(mz(combined_spectra), length)) == 0)
if(length(idx) > 0){
  combined_spectra <- combined_spectra[-idx]
}
combined_spectra <- combined_spectra[order(combined_spectra$entropy)]
ms2_NEG <- combined_spectra

rm(standards, combined_spectra, idx)

pdf("output/MSMS.pdf", paper = "a4", height = 4*3, width = 2*4)
par(mfrow = c(3, 2)) #put three if pdf
for(i in seq(nrow(cmps))){
  for(p in c("POS", "NEG")){
    if(!paste(cmps$name[i], p) %in% paste(excluded$name, excluded$POL)) {
      ms2 <- get(paste("ms2", p, sep = "_"))
      if(is.na(cmps[i, paste0(p, "_file")])){ #if it is not NA, do next lines
        ms2sub <- filterPrecursorMzRange(ms2, as.numeric(cmps[i,p]) + 0.05*c(-1, 1)) 
        ms2sub <- filterRt(ms2sub, as.numeric(cmps$RT[i])*60 + 10*c(-1, 1)) 
      } else {
        idx <- which(grepl(unlist(cmps[i, paste0(p, "_file")]),basename(ms2$dataOrigin))& 
                       ms2$scanIndex == unlist(cmps[i, paste0(p, "_scan")]))
        if(length(idx) == 1){
          ms2sub <- ms2[idx]
        } else {
          print(paste("Warning!! Check compound", cmps$name[i], p, "!!!"))
        }
      }
      if(cmps$level[i] == "L1"){
        stds <- get(paste("std", p, sep = "_"))
        if(is.na(cmps[i, paste0("std_", p, "_file")])) {
          stdsub <- filterPrecursorMzRange(stds, as.numeric(cmps[i,p]) + 0.05*c(-1, 1)) 
          stdsub <- filterRt(stdsub, as.numeric(cmps$RT[i])*60 + 10*c(-1, 1)) 
          if(paste(cmps$name[i],p) %in% paste(excluded_std$name, excluded_std$POL)){
            stdsub <- filterPrecursorMzRange(stdsub, (as.numeric(cmps[i,p]) + 1) + 0.05*c(-1, 1)) #precursor + 1, since noone will be it will be empty
          }
        } else {
          idx <- which(grepl(unlist(cmps[i, paste0("std_", p, "_file")]),basename(stds$dataOrigin))& 
                         stds$scanIndex == unlist(cmps[i, paste0("std_", p, "_scan")]))
          if(length(idx) == 1){
            stdsub <- stds[idx]
          } else {
            print(paste("Warning!! Check compound", cmps$name[i], p, "!!!"))
          }
        }
        
        yl <- c(-1.1, 1.1)
        if(length(stdsub) > 0){
          sps2 <- data.frame(
            mz = mz(stdsub)[[1]],
            intensity = intensity(stdsub)[[1]]
          )
          sps2$intrel <- sps2$intensity / max(sps2$intensity)
          mycol <- 2
        } else {
          yl <- c(0, 1.1)
          mycol <- 1
        }
      } else if(cmps$level[i] == "L2"){
        tmp <- unlist(cmps[i, paste0("std_", p, "_file")])
        if(!is.na(tmp)){
          sps2 <- read_xlsx("data/L2_MSMS.xlsx", sheet = tmp)
          yl <- c(-1.1, 1.1)
          sps2$intrel <- sps2$intensity / max(sps2$intensity)
          sps2$mz <- as.numeric(sps2$mz)
          mycol <- 4
        } else {
          yl <- c(0, 1.1)
          mycol <- 1
        }
      } else {
        yl <- c(0, 1.1)
        mycol <- 1
      }
      if(length(ms2sub) > 0){
        sps <- data.frame(
          mz = mz(ms2sub)[[1]],
          intensity = intensity(ms2sub)[[1]]
        )
        sps$intrel <- sps$intensity / max(sps$intensity)
        if(p == "POS"){
          ad <- "[M+H]+"
        } else if(p == "NEG"){
          ad <- "[M-H]-"
        }
        plot(sps$mz, sps$intrel, type = "h", ylim = yl,
             xlab = "m/z", ylab = "relative intensity", bty = "l", 
             xlim = c(60, precursorMz(ms2sub)[1] + 2),
             main = paste(cmps$name[i], ad, "\n",
                          sprintf("%.4f", precursorMz(ms2sub)[1]), "@", 
                          sprintf("%.2f", rtime(ms2sub)[1]/60)))
        idx <- which(sps$intrel > 0.1)
        text(sps$mz[idx], sps$intrel[idx], sprintf("%.4f", sps$mz[idx]), pos = 3, 
             cex = 0.7)
        if(exists("sps2")){
          abline(h = 0)
          lines(sps2$mz, sps2$intrel*(-1), col = 1, type = "h")
          idx <- which(sps2$intrel > 0.1)
          text(sps2$mz[idx], sps2$intrel[idx]*(-1), 
               sprintf("%.4f", sps2$mz[idx]), pos = 1, cex = 0.7, col = 1)
          jo <- join(x = sps[,1], y = sps2$mz, tolerance = 0.05)
          jo_idx <- jo$x[which(!is.na(jo$x) & !is.na(jo$y))]
          if(length(jo_idx) > 0){
            points(sps$mz[jo_idx], sps$intrel[jo_idx], col = mycol, type = "h")
            if(length(jo_idx) > 5){
              jo_idx <- jo_idx[which(sps$intrel[jo_idx] > 0.2)]
            }
            if(length(jo_idx) > 0){
              points(sps$mz[jo_idx], sps$intrel[jo_idx], col = mycol, pch = 16)
              text(sps$mz[jo_idx], sps$intrel[jo_idx], 
                   sprintf("%.4f", sps$mz[jo_idx]), pos = 3, cex = 0.7, col = mycol)
            }
          }
          jo_idx <- jo$y[which(!is.na(jo$x) & !is.na(jo$y))]
          if(length(jo_idx) > 0){
            points(sps2$mz[jo_idx], sps2$intrel[jo_idx]*(-1), col = mycol, type = "h")
            if(length(jo_idx) > 5){
              jo_idx <- jo_idx[which(sps$intrel[jo_idx] > 0.2)]
            }
            if(length(jo_idx) > 0){
              points(sps2$mz[jo_idx], sps2$intrel[jo_idx]*(-1), col = mycol, pch = 16)
              text(sps2$mz[jo_idx], sps2$intrel[jo_idx]*(-1), 
                   sprintf("%.4f", sps2$mz[jo_idx]), pos = 1, cex = 0.7, col = mycol)
            }
          }
        }
        
      } else {
        plot(0,0, type = "n", xlab = "", ylab = "", bty = "n", xaxt = "n", yaxt = "n")
      }
    } else {
      plot(0,0, type = "n", xlab = "", ylab = "", bty = "n", xaxt = "n", yaxt = "n")
    }
    if(exists("sps2")){
      rm(sps2)
    }
  }
}
dev.off() #closing pdf
