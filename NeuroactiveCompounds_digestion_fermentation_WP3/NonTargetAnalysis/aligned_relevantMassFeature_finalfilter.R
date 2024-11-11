library(Spectra)
library(MetaboCoreUtils)
library(xcms)
library(MsCoreUtils)
library(Rdisop)
library("dplyr")
library(openxlsx)
library(DT)

library(MetaboAnnotation)
#for aligning with exclusion list

data_tomato <- read.csv("RelevantMF_Blank_Tomato_finalFilter.csv")
data_tomato$RT <- data_tomato$RT*60
#table with exclusion list
data_inulin <- read.csv("RelevantMF_Inulin_Tomato_finalFilter.csv")
data_inulin$RT <- data_inulin$RT*60

#compound_db$excatmass <- calculateMass(compound_db$formula)
#compound_db$positive <- as.numeric(mass2mz(compound_db$excatmass, "[M+H]+"))
#compound_db$negativo <- as.numeric(mass2mz(compound_db$excatmass, "[M-H]-"))


#cros mz and rt
pks_rt_match <- matchValues(
  data_inulin, data_tomato, param = MzRtParam(ppm = 30, toleranceRt = 12),
  mzColname= c("mz","mz"),
   rtColname=c("RT", "RT")
)
pks_rt_match

#SAVE THE WHOLE TABLE
tb_rt_match <- data.frame(matchedData(pks_rt_match))
tb_rt_match$RT <- tb_rt_match$RT/60
write.csv(tb_rt_match, "MF_inulin_BLANK_TOMATO_ALIGNEMENT_POS_WHOLETABLE.csv")

#filter the table with only the matched ones
tb_rt_match <- data.frame(matchedData(pks_rt_match))
tb_rt_match_only <- tb_rt_match[!is.na(tb_rt_match$target_name),]
tb_rt_match_only$RT <- tb_rt_match_only$RT/60

#filter the table with only the NON matched ones
tb_rt_match <- data.frame(matchedData(pks_rt_match))
tb_rt_NOTmatch <- tb_rt_match[is.na(tb_rt_match$target_name),]
tb_rt_NOTmatch$RT <- tb_rt_NOTmatch$RT/60
#save table as xlxs file
write.xlsx(tb_rt_match_only, "MF_inulin_2CONSEQTime_fragmentMS2_POS.xlsx")
write.xlsx(tb_rt_NOTmatch, "MF_inulin_BLANK_TOMATO_ALIGNEMENT_POS_WHOLETABLE_NOTMATCHED.xlsx")


