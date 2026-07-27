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

### working code

# bed = generateRandomBed(nc = 2)
# head(bed, n = 2)
# 
# col_fun = colorRamp2(c(-1, 0, 1), c("green", "black", "red"))
# 
# 
# circos.initializeWithIdeogram(genome_bed)
# 
# 
# circos.genomicTrack(bed, numeric.column = 4, 
#                     panel.fun = function(region, value, ...) {
#                       # numeric.column is automatically passed to `circos.genomicPoints()`
#                       circos.genomicPoints(region, value, ...)
#                     })
# circos.genomicHeatmap(bed, col = col_fun, side = "inside", border = "white")
# 
# 
# 
# circos.initializeWithIdeogram(genome_bed)
# bed <- data.frame(chr="chr1",start=seq(1,1000000,by=10000),end=seq(5000,1000000,by=10000),value1=rnorm(100))
# col_fun = colorRamp2(c(-1, 0, 1), c("green", "black", "red"))
# circos.genomicHeatmap(bed, col = col_fun, side = "inside", border = "white")
# circos.clear()


circos.initializeWithIdeogram(genome_bed)
bed <- data.frame(chr="chr1",start=seq(1,300000000,by=1000000),end=seq(5000,300000000,by=1000000),value1=rnorm(300))
bed <- rbind(bed,
             data.frame(chr="chr2",start=seq(1,200000000,by=1000000),end=seq(5000,200000000,by=1000000),value1=rnorm(200)))
# bed <- bed %>%
#   mutate(clr=ifelse(value1 > 2, "red", "black"))
bed <- bed %>%
  mutate(value2 = ifelse(value1>.5,value1,NA),
         value3 = ifelse(value1<.5,value1,NA))

# circos.genomicLines(bed)
#                     
#                     
# circos.genomicTrackPlotRegion(bed, ylim = c(-3, 3),
#                     panel.fun = function(region, value, ...) {
#                       circos.genomicLines(region, value,pt.col=value$clr,pch=16)
#                     })

circos.genomicTrackPlotRegion(bed, ylim = c(-3, 3.1),
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value,numeric.column=3, col="black")
                                circos.genomicLines(region, value,numeric.column=2, col="red")
                              })
circos.clear()

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
combined_qtl <- read.csv("analyses/rqtl_combined_qtl_fitqtl.csv")


### Make function to make one data with LOD scores for one trait from all analyses
ril_results <- readRDS(file = "data/RIL_rqtl.RDS")
b73_results <- readRDS(file = "data/B73_hybrids_rqtl.RDS")
b73_mph_results <- readRDS(file = "data/B73_hybrids_mph_rqtl.RDS")
mo17_results <- readRDS(file = "data/Mo17_hybrids_rqtl.RDS")
mo17_mph_results <- readRDS(file = "data/Mo17_hybrids_mph_rqtl.RDS")

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
slb_combined_qtl <- combined_qtl %>%
  filter(grepl("SLB",trait)) %>%
  filter(!is.na(qtl)) %>%
  select(chr1,pos1) %>%
  distinct() %>%
  rename(chr = "chr1", start = "pos1") %>%
  mutate(chr = paste("chr", chr,sep=""))


make_lod_circos(slb_results,slb_perms, slb_combined_qtl)
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


gls_combined_qtl <- combined_qtl %>%
  filter(grepl("GLS",trait)) %>%
  filter(!is.na(qtl)) %>%
  select(chr1,pos1) %>%
  distinct() %>%
  rename(chr = "chr1", start = "pos1") %>%
  mutate(chr = paste("chr", chr,sep="")) %>%
  filter(chr != "chr5")


make_lod_circos(gls_results,gls_perms, gls_combined_qtl)
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


dta_combined_qtl <- combined_qtl %>%
  filter(grepl("DTA",trait)) %>%
  filter(!is.na(qtl)) %>%
  select(chr1,pos1) %>%
  distinct() %>%
  rename(chr = "chr1", start = "pos1") %>%
  mutate(chr = paste("chr", chr,sep="")) 


make_lod_circos(dta_results,dta_perms, dta_combined_qtl)
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

ph_combined_qtl <- combined_qtl %>%
  filter(grepl("PH_BLUE",trait)) %>%
  filter(!is.na(qtl)) %>%
  select(chr1,pos1) %>%
  distinct() %>%
  rename(chr = "chr1", start = "pos1") %>%
  mutate(chr = paste("chr", chr,sep="")) 

make_lod_circos(ph_results,ph_perms,ph_combined_qtl)
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
eh_combined_qtl <- combined_qtl %>%
  filter(grepl("EH_BLUE",trait)) %>%
  filter(!is.na(qtl)) %>%
  select(chr1,pos1) %>%
  distinct() %>%
  rename(chr = "chr1", start = "pos1") %>%
  mutate(chr = paste("chr", chr,sep="")) 

make_lod_circos(eh_results,eh_perms,eh_combined_qtl)
circos.clear()

make_lod_circos <- function(results,perms,combined_qtl){
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
  ymax <- lod_qtl %>% select(contains("sig")) %>% unlist() %>% max(na.rm=TRUE) %>% ceiling()
  
  # add column for QTL peaks
  combined_qtl2 <- combined_qtl
  combined_qtl2$peak <- TRUE
  lod_qtl <- merge(lod_qtl,combined_qtl2,all=TRUE, by = c("chr","start"))
  lod_qtl$peak[which(lod_qtl$peak==TRUE)] <- rep(ymax,length(which(lod_qtl$peak==TRUE)))
  
  lod_qtl <- lod_qtl %>%
    mutate(end=start + 1)
  
  circos.par("start.degree" = 70)
  circos.par("gap.degree" = c(rep(1,9),40))
  par(cex=1)
  circos.initializeWithIdeogram(IcM_bed,
                                plotType = c("axis", "labels"))
  #circos.text(x=1,y=0,adj = c(degree(317.5), degree(0)), labels="Chr")
  
  circos.genomicTrackPlotRegion(lod_qtl, ylim=c(0,ymax), track.height=0.1,
                                panel.fun = function(region, value, ...) {
                                  circos.genomicLines(region, value,numeric.column=1, col="black")
                                  circos.genomicLines(region, value,numeric.column=6, col="red",type="o", cex=0.20, pch=16, pt.col="red")
                                  circos.genomicLines(region, value,numeric.column=11, col="red",type="h", lty=3, ylim = c(0,ymax))
                                  
                                  
                                })
  
  circos.text(x=1,y=-1,adj = c(degree(317.5), degree(0)), labels="RIL", font=2,cex=.5)
  
  circos.genomicTrackPlotRegion(lod_qtl,ylim=c(0,ymax), track.height=0.1,
                                panel.fun = function(region, value, ...) {
                                  circos.genomicLines(region, value,numeric.column=2, col="black")
                                  circos.genomicLines(region, value,numeric.column=7, col="red",type="o", cex=0.20, pch=16, pt.col="red")
                                  circos.genomicLines(region, value,numeric.column=11, col="red",type="h", lty=3, ylim = c(0,ymax))
                                })
  
  circos.text(x=1,y=-2,adj = c(degree(317.5), degree(0)), labels="B73 BC", font=2,cex=.5)
  
  circos.genomicTrackPlotRegion(lod_qtl,ylim=c(0,ymax), track.height=0.1, 
                                panel.fun = function(region, value, ...) {
                                  circos.genomicLines(region, value,numeric.column=3, col="black")
                                  circos.genomicLines(region, value,numeric.column=8, col="red",type="o", cex=0.20, pch=16, pt.col="red")
                                  circos.genomicLines(region, value,numeric.column=11, col="red",type="h", lty=3, ylim = c(0,ymax))
                                })
  
  circos.text(x=1,y=-2,adj = c(degree(317.5), degree(0)), labels="B73 MPH", font=2,cex=.5)
  
  circos.genomicTrackPlotRegion(lod_qtl,ylim=c(0,ymax), track.height=0.1, 
                                panel.fun = function(region, value, ...) {
                                  circos.genomicLines(region, value,numeric.column=4, col="black")
                                  circos.genomicLines(region, value,numeric.column=9, col="red",type="o", cex=0.20, pch=16, pt.col="red")
                                  circos.genomicLines(region, value,numeric.column=11, col="red",type="h", lty=3, ylim = c(0,ymax))
                                })
  
  circos.text(x=1,y=-2,adj = c(degree(317.5), degree(0)), labels="Mo17 BC", font=2,cex=.5)
  
  circos.genomicTrackPlotRegion(lod_qtl,ylim=c(0,ymax), track.height=0.1, 
                                panel.fun = function(region, value, ...) {
                                  circos.genomicLines(region, value,numeric.column=5, col="black")
                                  circos.genomicLines(region, value,numeric.column=10, col="red",type="o", cex=0.20, pch=16, pt.col="red")
                                  circos.genomicLines(region, value,numeric.column=11, col="red",type="h", lty=3, ylim = c(0,ymax))
                                })
  
  circos.text(x=1,y=-2,adj = c(degree(317.5), degree(0)), labels="Mo17 MPH", font=2,cex=.5)
  
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

