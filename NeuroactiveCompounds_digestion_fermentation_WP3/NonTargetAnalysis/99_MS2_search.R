library(Spectra)
scale_int <- function(x, ...) {
  maxint <- max(x[, "intensity"], na.rm = TRUE)
  x[, "intensity"] <- 100 * x[, "intensity"] / maxint
  x
}

p <- "POS"
load(paste0("rdata/", p, "_MS2spectra_tomato_corrected.RData"))


# precursor m/z - RT -----

mymz <- 914.5972
myrt <- 9.93

(ms2sub <- filterPrecursorMzRange(ms2, mymz + 0.05 * c(-1, 1))) #HIGH PPM of precursor because of DDA acquisition in centorid
(ms2sub <- filterRt(ms2sub, myrt*60 + 10 * c(-1, 1)))
#(ms2sub <- ms2sub[grep("45_PS", basename(ms2sub$dataOrigin))])


# plot -----

ms2sub <- ms2sub[order(ms2sub$entropy, ms2sub$entropy_norm)]
ms2sub <- ms2sub[1:40]
for(i in rev(seq(length(ms2sub)))){
  sps <- data.frame(
    mz = mz(ms2sub)[[i]],
    i = intensity(ms2sub)[[i]]
  )
  if(nrow(sps) > 0){
    plot(sps$mz, sps$i, type = "h", bty = "l", xlab = "m/z", ylab = "intensity", 
         xlim = c(70, precursorMz(ms2sub)[i]+20),
         main = paste(
           "i", i, "-",
           round(precursorMz(ms2sub)[i], 4), "@", round(rtime(ms2sub)[i]/60, 2),
           "\n", gsub(".*RES_", "RES_", 
                      basename(ms2sub$dataOrigin[i])), "-", ms2sub$scanIndex[i]#,
           #"\n entropy = ", round(ms2sub$entropy[i], 2), 
           #"- norm =", round(ms2sub$entropy_norm[i], 2)
         ))
    idx <- which(sps$i / max(sps$i) > 0.1)
    text(sps$mz[idx], sps$i[idx], round(sps$mz[idx], 4), cex = 0.7)
    print(paste(
      "i", i, "-",
      "\n", basename(ms2sub$dataOrigin[i]), "-", ms2sub$scanIndex[i]))
  }
}



# contains m/z ----
ms2n <- addProcessing(ms2, scale_int)
ms2n <- filterIntensity(ms2n, intensity = c(10, Inf))

mymz <- 191.0506
myrt <- 5.93
(ms2sub <- filterRt(ms2n, myrt*60 + 10 * c(-1, 1)))
(ms2sub <- ms2sub[containsMz(ms2sub, mz = mymz, tolerance = 0.01)])
