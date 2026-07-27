### make circos plot of QTL
### try final goal of one track for each analysis, one plot for each trait
library(gap)
library(circlize)
library(qtl)
circos.initializeWithIdeogram()
text(0, 0, "default", cex = 1)
circos.info()
circos.clear()
genome_bed <- read.table("/Users/asherhudson/Downloads/b73_v5.bed")

### load RIL QTL LOD scores

ril_results <-readRDS(file = "data/RIL_rqtl.RDS")
str(ril_results[[1]][[1]])
result <- ril_results[[1]][[1]]
perms <- mqmprocesspermutation(ril_results[[2]][[1]])
thrs <- summary(perms)
lod_qtl <- data.frame(chr=result$chr,start=result$`pos (cM)`,end=result$`pos (cM)`+1,lod=result$`LOD SLB_WMD_BLUE`)
lod_qtl <- lod_qtl %>%
  mutate(value2 = ifelse(lod>thrs[1],lod,NA),
         value3 = ifelse(lod<thrs[1],lod,NA))
lod_qtl$chr <- paste("chr", lod_qtl$chr,sep="")
IcM_bed <- lod_qtl %>%
  select(chr,end) %>%
  group_by(chr) %>%
  filter(end==max(end)) %>%
  mutate(start=0) %>%
  relocate(start, .after = "chr")
IcM_bed <- data.frame(IcM_bed)
circos.initializeWithIdeogram(IcM_bed)
circos.genomicTrackPlotRegion(lod_qtl,
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value,numeric.column=3, col="black")
                                circos.genomicLines(region, value,numeric.column=2, col="red")
                              })

### Make function to make one data with LOD scores for one trait from all analyses
ril_results <- readRDS(file = "data/RIL_rqtl.RDS")
b73_results <- readRDS(file = "data/B73_hybrids_rqtl.RDS")
b73_mph_results <- readRDS(file = "data/B73_hybrids_mph_rqtl.RDS")
mo17_results <- readRDS(file = "data/Mo17_hybrids_rqtl.RDS")
mo17_mph_results <- readRDS(file = "data/Mo17_hybrids_mph_rqtl.RDS")

effect_df <- read.csv("analyses/qtl_effects_whole_genome.csv")

slb_results <- list(ril_results[[1]][[1]],
                    b73_results[[1]][[1]],
                    b73_mph_results[[1]][[2]],
                    mo17_results[[1]][[1]],
                    mo17_mph_results[[1]][[1]])
slb_perms <- list(ril_results[[2]][[1]],
                  b73_results[[2]][[1]],
                  b73_mph_results[[2]][[2]],
                  mo17_results[[2]][[1]],
                  mo17_mph_results[[2]][[1]])
slb_effect_df <- effect_df %>% 
  select(contains("chr") | contains("pos..cM.") | 
           (contains("SLB") & contains("effect")))  %>%
  select(-c(SLB_WMD_BLUE_ril_effect_centered,SLB_WMD_BLUE_b73_mph_effect))
slb_combined_qtl <- combined_qtl %>%
  filter(grepl("SLB",trait)) %>%
  filter(!is.na(qtl)) %>%
  select(chr1,pos1) %>%
  distinct() %>%
  rename(chr = "chr1", pos..cM. = "pos1") %>%
  mutate(chr = paste("chr", chr,sep=""))

make_effect_circos(slb_results,slb_perms,slb_effect_df,"SLB_WMD_BLUE", slb_combined_qtl)
circos.clear()

gls_results <- list(ril_results[[1]][[2]],
                    b73_results[[1]][[2]],
                    b73_mph_results[[1]][[4]],
                    mo17_results[[1]][[2]],
                    mo17_mph_results[[1]][[2]])
gls_perms <- list(ril_results[[2]][[2]],
                  b73_results[[2]][[2]],
                  b73_mph_results[[2]][[4]],
                  mo17_results[[2]][[2]],
                  mo17_mph_results[[2]][[2]])


gls_effect_df <- effect_df %>% 
  select(contains("chr") | contains("pos..cM.") | 
           (contains("GLS") & contains("effect")))  %>%
  select(-c(GLS_WMD_BLUE_b73_mph_effect))
gls_combined_qtl <- combined_qtl %>%
  filter(grepl("GLS",trait)) %>%
  filter(!is.na(qtl)) %>%
  select(chr1,pos1) %>%
  distinct() %>%
  rename(chr = "chr1", pos..cM. = "pos1") %>%
  mutate(chr = paste("chr", chr,sep=""))
gls_combined_qtl <- gls_combined_qtl %>%
  filter(chr!="chr5")

make_effect_circos(gls_results,gls_perms,gls_effect_df,"GLS_WMD_BLUE", gls_combined_qtl)

circos.clear()

dta_results <- list(ril_results[[1]][[3]],
                    b73_results[[1]][[3]],
                    b73_mph_results[[1]][[6]],
                    mo17_results[[1]][[3]],
                    mo17_mph_results[[1]][[3]])
dta_perms <- list(ril_results[[2]][[3]],
                  b73_results[[2]][[3]],
                  b73_mph_results[[2]][[6]],
                  mo17_results[[2]][[3]],
                  mo17_mph_results[[2]][[3]])

dta_effect_df <- effect_df %>% 
  select(contains("chr") | contains("pos..cM.") | 
           (contains("DTA") & contains("effect")))  %>%
  select(-c(DTA_BLUE_b73_mph_effect))
dta_combined_qtl <- combined_qtl %>%
  filter(grepl("DTA",trait)) %>%
  filter(!is.na(qtl)) %>%
  select(chr1,pos1) %>%
  distinct() %>%
  rename(chr = "chr1", pos..cM. = "pos1") %>%
  mutate(chr = paste("chr", chr,sep=""))


make_effect_circos(dta_results,dta_perms,dta_effect_df,"DTA_BLUE", dta_combined_qtl)

circos.clear()

ph_results <- list(ril_results[[1]][[4]],
                   b73_results[[1]][[4]],
                   b73_mph_results[[1]][[8]],
                   mo17_results[[1]][[4]],
                   mo17_mph_results[[1]][[4]])
ph_perms <- list(ril_results[[2]][[4]],
                 b73_results[[2]][[4]],
                 b73_mph_results[[2]][[8]],
                 mo17_results[[2]][[4]],
                 mo17_mph_results[[2]][[4]])


ph_effect_df <- effect_df %>% 
  select(contains("chr") | contains("pos..cM.") | 
           (contains("PH_BLUE") & contains("effect")))  %>%
  select(-c(PH_BLUE_b73_mph_effect))
ph_combined_qtl <- combined_qtl %>%
  filter(grepl("PH_BLUE",trait)) %>%
  filter(!is.na(qtl)) %>%
  select(chr1,pos1) %>%
  distinct() %>%
  rename(chr = "chr1", pos..cM. = "pos1") %>%
  mutate(chr = paste("chr", chr,sep=""))


make_effect_circos(ph_results,ph_perms,ph_effect_df,"PH_BLUE", ph_combined_qtl)

circos.clear()

eh_results <- list(ril_results[[1]][[5]],
                   b73_results[[1]][[5]],
                   b73_mph_results[[1]][[10]],
                   mo17_results[[1]][[5]],
                   mo17_mph_results[[1]][[5]])
eh_perms <- list(ril_results[[2]][[5]],
                 b73_results[[2]][[5]],
                 b73_mph_results[[2]][[10]],
                 mo17_results[[2]][[5]],
                 mo17_mph_results[[2]][[5]])


eh_effect_df <- effect_df %>% 
  select(contains("chr") | contains("pos..cM.") | 
           (contains("EH_BLUE") & contains("effect")))  %>%
  select(-c(EH_BLUE_b73_mph_effect))
eh_combined_qtl <- combined_qtl %>%
  filter(grepl("EH_BLUE",trait)) %>%
  filter(!is.na(qtl)) %>%
  select(chr1,pos1) %>%
  distinct() %>%
  rename(chr = "chr1", pos..cM. = "pos1") %>%
  mutate(chr = paste("chr", chr,sep=""))


make_effect_circos(eh_results,eh_perms,eh_effect_df,"EH_BLUE", eh_combined_qtl)
circos.clear()

make_effect_circos <- function(results,perms,effect_df,trait, combined_qtl){
  results <- lapply(results, function(x) data.frame(x))
  lod_qtl <- data.frame(chr=results[[1]][[1]],start=results[[1]][[2]],end=results[[1]][[2]]+1,
                        ril_lod=results[[1]][[3]],
                        b73_lod=results[[2]][[3]],
                        b73_mph_lod=results[[3]][[3]],
                        mo17_lod=results[[4]][[3]],
                        mo17_mph_lod=results[[5]][[3]])
  lod_qtl <- lod_qtl %>%
    mutate(ril_lod_sig=ril_lod,
           b73_lod_sig=b73_lod,
           b73_mph_lod_sig=b73_mph_lod,
           mo17_lod_sig=mo17_lod,
           mo17_mph_lod_sig=mo17_mph_lod)
  
  perms <- lapply(perms,function(x) mqmprocesspermutation(x))
  perms <- lapply(perms, function(x) summary(x))
  
  lod_qtl <- lod_qtl %>%
    mutate(ril_lod_sig = ifelse(ril_lod_sig>(perms[[1]][[1]]),ril_lod_sig,NA),
           ril_lod = ifelse(ril_lod<(perms[[1]][[1]]),ril_lod,NA),
           b73_lod_sig = ifelse(b73_lod_sig>perms[[2]][[1]],b73_lod_sig,NA),
           b73_lod = ifelse(b73_lod<perms[[2]][[1]],b73_lod,NA),
           b73_mph_lod_sig = ifelse(b73_mph_lod_sig>perms[[3]][[1]],b73_mph_lod_sig,NA),
           b73_mph_lod = ifelse(b73_mph_lod<perms[[3]][[1]],b73_mph_lod,NA),
           mo17_lod_sig = ifelse(mo17_lod_sig>perms[[4]][[1]],mo17_lod_sig,NA),
           mo17_lod = ifelse(mo17_lod<perms[[4]][[1]],mo17_lod,NA),
           mo17_mph_lod_sig = ifelse(mo17_mph_lod_sig>perms[[5]][[1]],mo17_mph_lod_sig,NA),
           mo17_mph_lod = ifelse(mo17_mph_lod<perms[[5]][[1]],mo17_mph_lod,NA))
  lod_qtl <- data.frame(lod_qtl)
  lod_qtl$chr <- paste("chr",lod_qtl$chr,sep="")
  
  
  # convert effects to NA depending on significance
  
  tmp=grep("effect",colnames(effect_df))
  
  effect_df2 <- cbind(
    effect_df,
    setNames(effect_df[,tmp],paste0(colnames(effect_df)[tmp],"no2"))
  )
  
  effect_df2[which(is.na(lod_qtl$ril_lod)),grep("ril_effect$",colnames(effect_df2))] <- NA
  effect_df2[which(is.na(lod_qtl$b73_lod)),grep("b73_effect$",colnames(effect_df2))] <- NA
  effect_df2[which(is.na(lod_qtl$b73_mph_lod)),grep("b73_mph_effect$",colnames(effect_df2))] <- NA
  effect_df2[which(is.na(lod_qtl$mo17_lod)),grep("mo17_effect$",colnames(effect_df2))] <- NA
  effect_df2[which(is.na(lod_qtl$mo17_mph_lod)),grep("mo17_mph_effect$",colnames(effect_df2))] <- NA
  
  effect_df2[which(is.na(lod_qtl$ril_lod_sig)),grep("ril_effectno2",colnames(effect_df2))] <- NA
  effect_df2[which(is.na(lod_qtl$b73_lod_sig)),grep("b73_effectno2",colnames(effect_df2))] <- NA
  effect_df2[which(is.na(lod_qtl$b73_mph_lod_sig)),grep("b73_mph_effectno2",colnames(effect_df2))] <- NA
  effect_df2[which(is.na(lod_qtl$mo17_lod_sig)),grep("mo17_effectno2",colnames(effect_df2))] <- NA
  effect_df2[which(is.na(lod_qtl$mo17_mph_lod_sig)),grep("mo17_mph_effectno2",colnames(effect_df2))] <- NA
  
  ymax <- effect_df2 %>% select(contains(trait) & contains("effect")) %>% unlist() %>% abs() %>% max(na.rm=TRUE) 
  ymax <- ymax*1.1
  
  effect_df2 <- effect_df2 %>%
    mutate(end=pos..cM. + 1)
  effect_df2 <- effect_df2 %>%
    relocate(end, .after = pos..cM.)
  effect_df2$chr <- paste("chr", effect_df2$chr, sep="")
  
  combined_qtl2 <- combined_qtl
  combined_qtl2$peak <- TRUE
  effect_df2 <- merge(effect_df2,combined_qtl2,all=TRUE)
  effect_df2$peak[which(effect_df2$peak==TRUE)] <- rep(ymax,length(which(effect_df2$peak==TRUE)))
  
  effect_df2 <- effect_df2 %>%
    mutate(end=pos..cM. + 1)
  
  circos.par("start.degree" = 70)
  circos.par("gap.degree" = c(rep(1,9),40))
  par(cex=1)
  circos.initializeWithIdeogram(IcM_bed,
                                plotType = c("axis", "labels"))
  #circos.text(x=1,y=0,adj = c(degree(317.5), degree(0)), labels="Chr")
  
  circos.genomicTrackPlotRegion(effect_df2, ylim=c(-ymax,ymax), track.height=0.1,
                                panel.fun = function(region, value, ...) {
                                  circos.genomicLines(region, value,numeric.column=1, col="black")
                                  circos.genomicLines(region, value,numeric.column=6, col="red",type="o", cex=0.20, pch=16, pt.col="red")
                                  circos.genomicLines(region, value,numeric.column=11, col="red",type="h", lty=3, ylim = c(-ymax,ymax))
                                  circos.segments(x0=0,y0=0,x1=CELL_META$xlim,y1=0,lty=2)
                                  
                                })
  
  circos.text(x=1,y=0,adj = c(degree(317.5), degree(0)), labels="RIL", font=2,cex=.5)
  
  circos.genomicTrackPlotRegion(effect_df2,ylim=c(-ymax,ymax), track.height=0.1,
                                panel.fun = function(region, value, ...) {
                                  circos.genomicLines(region, value,numeric.column=2, col="black")
                                  circos.genomicLines(region, value,numeric.column=7, col="red",type="o", cex=0.20, pch=16, pt.col="red")
                                  circos.genomicLines(region, value,numeric.column=11, col="red",type="h", lty=3, ylim = c(-ymax,ymax))
                                  
                                  circos.segments(x0=0,y0=0,x1=CELL_META$xlim,y1=0,lty=2)
                                  
                                })
  
  circos.text(x=1,y=0,adj = c(degree(317.5), degree(0)), labels="B73 BC", font=2,cex=.5)
  
  circos.genomicTrackPlotRegion(effect_df2,ylim=c(-ymax,ymax), track.height=0.1, 
                                panel.fun = function(region, value, ...) {
                                  circos.genomicLines(region, value,numeric.column=4, col="black")
                                  circos.genomicLines(region, value,numeric.column=9, col="red",type="o", cex=0.20, pch=16, pt.col="red")
                                  circos.genomicLines(region, value,numeric.column=11, col="red",type="h", lty=3, ylim = c(-ymax,ymax))
                                  
                                  circos.segments(x0=0,y0=0,x1=CELL_META$xlim,y1=0,lty=2)
                                })
  
  circos.text(x=1,y=0,adj = c(degree(317.5), degree(0)), labels="B73 MPH", font=2,cex=.5)
  
  circos.genomicTrackPlotRegion(effect_df2,ylim=c(-ymax,ymax), track.height=0.1, 
                                panel.fun = function(region, value, ...) {
                                  circos.genomicLines(region, value,numeric.column=3, col="black")
                                  circos.genomicLines(region, value,numeric.column=8, col="red",type="o", cex=0.20, pch=16, pt.col="red")
                                  circos.genomicLines(region, value,numeric.column=11, col="red",type="h", lty=3, ylim = c(-ymax,ymax))
                                  
                                   circos.segments(x0=0,y0=0,x1=CELL_META$xlim,y1=0,lty=2)
                                })
  
  circos.text(x=1,y=0,adj = c(degree(317.5), degree(0)), labels="Mo17 BC", font=2,cex=.5)
  
  circos.genomicTrackPlotRegion(effect_df2,ylim=c(-ymax,ymax), track.height=0.1, 
                                panel.fun = function(region, value, ...) {
                                  circos.genomicLines(region, value,numeric.column=5, col="black")
                                  circos.genomicLines(region, value,numeric.column=10, col="red",type="o", cex=0.20, pch=16, pt.col="red")
                                  circos.genomicLines(region, value,numeric.column=11, col="red",type="h", lty=3, ylim = c(-ymax,ymax))
                                  
                                  circos.segments(x0=0,y0=0,x1=CELL_META$xlim,y1=0,lty=2)
                                })
  
  circos.text(x=1,y=0,adj = c(degree(317.5), degree(0)), labels="Mo17 MPH", font=2,cex=.5)
  
}


results <- slb_results
perms <- slb_perms
results <- lapply(results, function(x) data.frame(x))
lod_qtl <- data.frame(chr=results[[1]][[1]],start=results[[1]][[2]],end=results[[1]][[2]]+1,
                      ril_lod=results[[1]][[3]],
                      b73_lod=results[[2]][[3]],
                      b73_mph_lod=results[[3]][[3]],
                      mo17_lod=results[[4]][[3]],
                      mo17_mph_lod=results[[5]][[3]])
lod_qtl <- lod_qtl %>%
  mutate(ril_lod_sig=ril_lod,
         b73_lod_sig=b73_lod,
         b73_mph_lod_sig=b73_mph_lod,
         mo17_lod_sig=mo17_lod,
         mo17_mph_lod_sig=mo17_mph_lod)

perms <- lapply(perms,function(x) mqmprocesspermutation(x))
perms <- lapply(perms, function(x) summary(x))

lod_qtl <- lod_qtl %>%
  mutate(ril_lod_sig = ifelse(ril_lod_sig>(perms[[1]][[1]]),ril_lod_sig,NA),
         ril_lod = ifelse(ril_lod<(perms[[1]][[1]]),ril_lod,NA),
         b73_lod_sig = ifelse(b73_lod_sig>perms[[2]][[1]],b73_lod_sig,NA),
         b73_lod = ifelse(b73_lod<perms[[2]][[1]],b73_lod,NA),
         b73_mph_lod_sig = ifelse(b73_mph_lod_sig>perms[[3]][[1]],b73_mph_lod_sig,NA),
         b73_mph_lod = ifelse(b73_mph_lod<perms[[3]][[1]],b73_mph_lod,NA),
         mo17_lod_sig = ifelse(mo17_lod_sig>perms[[4]][[1]],mo17_lod_sig,NA),
         mo17_lod = ifelse(mo17_lod<perms[[4]][[1]],mo17_lod,NA),
         mo17_mph_lod_sig = ifelse(mo17_mph_lod_sig>perms[[5]][[1]],mo17_mph_lod_sig,NA),
         mo17_mph_lod = ifelse(mo17_mph_lod<perms[[5]][[1]],mo17_mph_lod,NA))
lod_qtl <- data.frame(lod_qtl)
lod_qtl$chr <- paste("chr",lod_qtl$chr,sep="")


circos.par("start.degree" = 70)
circos.par("gap.degree" = c(rep(1,9),40))
par(cex=1)
circos.initializeWithIdeogram(IcM_bed,
                              plotType = c("axis", "labels"))
#circos.text(x=1,y=0,adj = c(degree(317.5), degree(0)), labels="Chr")

circos.genomicTrackPlotRegion(lod_qtl, ylim=c(0,9), track.height=0.1,
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value,numeric.column=1, col="black")
                                circos.genomicLines(region, value,numeric.column=6, col="red",type="o", cex=0.20, pch=16, pt.col="red")
                                
                              })

circos.text(x=1,y=-1,adj = c(degree(317.5), degree(0)), labels="RIL", font=2,cex=.5)

circos.genomicTrackPlotRegion(lod_qtl,ylim=c(0,9), track.height=0.1,
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value,numeric.column=2, col="black")
                                circos.genomicLines(region, value,numeric.column=7, col="red",type="o", cex=0.20, pch=16, pt.col="red")
                              })

circos.text(x=1,y=-2,adj = c(degree(317.5), degree(0)), labels="B73 BC", font=2,cex=.5)

circos.genomicTrackPlotRegion(lod_qtl,ylim=c(0,9), track.height=0.1, 
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value,numeric.column=3, col="black")
                                circos.genomicLines(region, value,numeric.column=8, col="red",type="o", cex=0.20, pch=16, pt.col="red")
                              })

circos.text(x=1,y=-2,adj = c(degree(317.5), degree(0)), labels="B73 MPH", font=2,cex=.5)

circos.genomicTrackPlotRegion(lod_qtl,ylim=c(0,9), track.height=0.1, 
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value,numeric.column=4, col="black")
                                circos.genomicLines(region, value,numeric.column=9, col="red",type="o", cex=0.20, pch=16, pt.col="red")
                              })

circos.text(x=1,y=-2,adj = c(degree(317.5), degree(0)), labels="Mo17 BC", font=2,cex=.5)

circos.genomicTrackPlotRegion(lod_qtl,ylim=c(0,9), track.height=0.1, 
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value,numeric.column=5, col="black")
                                circos.genomicLines(region, value,numeric.column=10, col="red",type="o", cex=0.20, pch=16, pt.col="red")
                              })

circos.text(x=1,y=-2,adj = c(degree(317.5), degree(0)), labels="Mo17 MPH", font=2,cex=.5)



### Make a gap and add labels for either letters or analysis name
### Maybe put gap between chr 1 and 10 at top
### Change bp to IcM

