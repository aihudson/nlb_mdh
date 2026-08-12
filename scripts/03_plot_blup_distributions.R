# libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(cowplot)
library(gridGraphics)

library(conflicted)
conflicts_prefer(base::as.data.frame, dplyr::select, dplyr::filter, dplyr::mutate)

# read in data
blups <- read.csv("analyses/IBM_NLB_BLUPs.csv")

# make datasets for RILs and backcrosses
ril_blups <- blups %>%
  filter(grepl("M0", Line)) %>%
  filter(!grepl("x", Line))
ril_blups$Population <- "RIL"
mo17_blups <- blups %>%
  filter(grepl("Mo17", Line)) %>%
  filter(grepl("x", Line))
mo17_blups$Population <- "Mo17 BC"
b73_blups <- blups %>%
  filter(grepl("B73", Line)) %>%
  filter(grepl("x", Line))
b73_blups$Population <- "B73 BC"

mo17_blup <- blups %>%
  filter(Line == "Mo17") 
mo17_blup$Population <- "Mo17"
b73_blup <- blups %>%
  filter(Line == "B73") 
b73_blup$Population <- "B73"
bxm_blup <- blups %>%
  filter(Line == "B73xMo17") 
bxm_blup$Population <- "B73 X Mo17"

blups <- rbind(ril_blups,mo17_blups,b73_blups)
original_lines <- rbind(mo17_blup, b73_blup, bxm_blup)

blups_long <- gather(blups, condition, blup, NLB_WMD_BLUP, factor_key=TRUE)


make_figure <- function(df, original_lines, trait){
  trait_name <- case_when(
    trait == "NLB_WMD_BLUP" ~ "NLB resistance BLUP (100 - WMD)"
  )
  bw <- 2 * IQR(blups[[trait]], na.rm=TRUE) / length(which(!is.na(blups[[trait]])))^(1/3)
  b73_line <- original_lines %>%
  filter(Line == "B73") %>%
  select(!!trait)
  b73_line <- as.numeric(b73_line)
  mo17_line <- original_lines %>%
  filter(Line == "Mo17") %>%
  select(!!trait)
  mo17_line <- as.numeric(mo17_line)
  bxm_line <- original_lines %>%
   filter(Line == "B73xMo17") %>%
   select(!!trait)
  bxm_line <- as.numeric(bxm_line)
  group.colors <- c(`B73 BC` = "#F8766D", `RIL` = "#00BA38", `Mo17 BC` ="#619CFF")
  
  plot <- ggplot(blups, aes(.data[[trait]], colour = Population)) +
    geom_freqpoly(binwidth=bw, linewidth = 2) +
    labs(x = trait_name,
         y = "Count") +
    theme(text = element_text(size=20)) +
    scale_color_manual(values=group.colors) 
  # geom_vline(xintercept = b73_line, color = "black") +
  #  annotate("text", label = "B73", x=b73_line+3, y = 51) +
  # geom_vline(xintercept = mo17_line, color = "black") +
  #   annotate("text", label = "Mo17", x=mo17_line-3, y = 51) +
  # geom_vline(xintercept = bxm_line, color = "black") +
  #   annotate("text", label = "B73xMo17", x=bxm_line+4.5, y = 51)
  return(plot)
}



(plot_list <- lapply(list("NLB_WMD_BLUP"), function(x){ make_figure(blups, original_lines, x) }))
plot <- plot_list[[1]]

pdf("figures/blup_distributions.pdf",height=10)
plot
dev.off()

# plot heterosis


### MPH

make_figure <- function(df, original_lines, trait){
  #trait_name <- case_when(
  #trait == "SLB_WMD_BLUE_MPH_PCT" ~ "Southern leaf blight resistance",
  #trait == "GLS_WMD_BLUE_MPH_PCT" ~ "Gray leaf spot resistance",
  #trait == "DTA_BLUE_MPH_PCT" ~ "Days to anthesis",
  # trait == "PH_BLUE_MPH_PCT" ~ "Plant height",
  # trait == "EH_BLUE_MPH_PCT" ~ "Ear height"
  #)
  trait_name <- case_when(
    trait == "SLB_WMD_BLUE_MPH_PCT" ~ "SLB",
    trait == "NLB_WMD_BLUP_MPH_PCT" ~ "NLB resistance MPH (%)",
    trait == "GLS_WMD_BLUE_MPH_PCT" ~ "GLS",
    trait == "DTA_BLUE_MPH_PCT" ~ "DTA",
    trait == "PH_BLUE_MPH_PCT" ~ "PH",
    trait == "EH_BLUE_MPH_PCT" ~ "EH"
  )
  bw <- 2 * IQR(df[[trait]], na.rm=TRUE) / length(which(!is.na(df[[trait]])))^(1/3)
  group.colors <- c(`B73 BC` = "#F8766D", `RIL` = "#00BA38", `Mo17 BC` ="#619CFF")
  plot <- ggplot(df, aes(.data[[trait]], colour = Population)) +
    geom_freqpoly(binwidth=bw, linewidth = 2, show.legend = FALSE) +
    labs(x = trait_name,
         y = "Count",
         title = "MPH") +
    theme(text = element_text(size=20),
          title = element_text(size=14),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_x_continuous(labels = scales::percent, limits = c(-0.3,0.3)) +
    scale_color_manual(values=group.colors)
  return(plot)
}



plot_list <- lapply(list("NLB_WMD_BLUP_MPH_PCT"), function(x){ make_figure(blups, original_lines, x) })


### BPH

make_figure_bph <- function(df, original_lines, trait){
  #trait_name <- case_when(
  #trait == "SLB_WMD_BLUE_MPH_PCT" ~ "Southern leaf blight resistance",
  #trait == "GLS_WMD_BLUE_MPH_PCT" ~ "Gray leaf spot resistance",
  #trait == "DTA_BLUE_MPH_PCT" ~ "Days to anthesis",
  # trait == "PH_BLUE_MPH_PCT" ~ "Plant height",
  # trait == "EH_BLUE_MPH_PCT" ~ "Ear height"
  #)
  trait_name <- case_when(
    trait == "SLB_WMD_BLUE_BPH_PCT" ~ "SLB",
    trait == "NLB_WMD_BLUP_BPH_PCT" ~ "NLB resistance BPH (%)",
    trait == "GLS_WMD_BLUE_BPH_PCT" ~ "GLS",
    trait == "DTA_BLUE_BPH_PCT" ~ "DTA",
    trait == "PH_BLUE_BPH_PCT" ~ "PH",
    trait == "EH_BLUE_BPH_PCT" ~ "EH"
  )
  bw <- 2 * IQR(df[[trait]], na.rm=TRUE) / length(which(!is.na(df[[trait]])))^(1/3)
  group.colors <- c(`B73 BC` = "#F8766D", `RIL` = "#00BA38", `Mo17 BC` ="#619CFF")
  plot <- ggplot(df, aes(.data[[trait]], colour = Population)) +
    geom_freqpoly(binwidth=bw, linewidth = 2, show.legend = FALSE) +
    labs(x = trait_name,
         y = "Count",
         title = "BPH") +
    theme(text = element_text(size=20),
          title = element_text(size=14),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_x_continuous(labels = scales::percent, limits = c(-0.3,0.3)) +
    scale_color_manual(values=group.colors)
  return(plot)
}

bph_plot_list <- lapply(list("NLB_WMD_BLUP_BPH_PCT"), function(x){ make_figure_bph(blups, original_lines, x) })

grab_legend <- function(df, original_lines, trait){
  #trait_name <- case_when(
  #trait == "SLB_WMD_BLUE_MPH_PCT" ~ "Southern leaf blight resistance",
  #trait == "GLS_WMD_BLUE_MPH_PCT" ~ "Gray leaf spot resistance",
  #trait == "DTA_BLUE_MPH_PCT" ~ "Days to anthesis",
  # trait == "PH_BLUE_MPH_PCT" ~ "Plant height",
  # trait == "EH_BLUE_MPH_PCT" ~ "Ear height"
  #)
  trait_name <- case_when(
    trait == "SLB_WMD_BLUE_BPH_PCT" ~ "SLB",
    trait == "NLB_WMD_BLUP_BPH_PCT" ~ "NLB resistance BPH (%)",
    trait == "GLS_WMD_BLUE_BPH_PCT" ~ "GLS",
    trait == "DTA_BLUE_BPH_PCT" ~ "DTA",
    trait == "PH_BLUE_BPH_PCT" ~ "PH",
    trait == "EH_BLUE_BPH_PCT" ~ "EH"
  )
  bw <- 2 * IQR(df[[trait]], na.rm=TRUE) / length(which(!is.na(df[[trait]])))^(1/3)
  group.colors <- c(`B73 BC` = "#F8766D", `RIL` = "#00BA38", `Mo17 BC` ="#619CFF")
  plot <- ggplot(df, aes(.data[[trait]], colour = Population)) +
    geom_freqpoly(binwidth=bw, linewidth = 2) +
    scale_color_manual(values=group.colors)

  legend <- get_legend(plot)
  return(legend)
}

legend <- lapply(list("NLB_WMD_BLUP_BPH_PCT"), function(x){ grab_legend(blups, original_lines, x) })


### Combining MPH and BPH in one figure


(figure <- plot_grid(plot_list[[1]], bph_plot_list[[1]], legend[[1]], ncol = 3, rel_widths = c(2, 2, 1))  + 
  theme(plot.margin = margin(0, 1, 0, 1, "cm"))
)


 figure_combined <- ggdraw(figure) + 
  draw_label("Count", x = 0, y = 0.5, vjust = 1.5, angle = 90, size = 14)

ggsave("figures/mph_bph_distributions.pdf", figure_combined, height=3)

make_figure <- function(df, original_lines, trait){
  #trait_name <- case_when(
  #trait == "SLB_WMD_BLUE_MPH_PCT" ~ "Southern leaf blight resistance",
  #trait == "GLS_WMD_BLUE_MPH_PCT" ~ "Gray leaf spot resistance",
  #trait == "DTA_BLUE_MPH_PCT" ~ "Days to anthesis",
  # trait == "PH_BLUE_MPH_PCT" ~ "Plant height",
  # trait == "EH_BLUE_MPH_PCT" ~ "Ear height"
  #)
  trait_name <- case_when(
    trait == "SLB_WMD_BLUE_MPH_PCT" ~ "SLB",
    trait == "GLS_WMD_BLUE_MPH_PCT" ~ "GLS",
    trait == "DTA_BLUE_MPH_PCT" ~ "DTA",
    trait == "PH_BLUE_MPH_PCT" ~ "PH",
    trait == "EH_BLUE_MPH_PCT" ~ "EH"
  )
  bw <- 2 * IQR(df[[trait]], na.rm=TRUE) / length(which(!is.na(df[[trait]])))^(1/3)
  bxm_line <- original_lines %>%
    filter(Line == "B73 X Mo17") %>%
    select(!!trait)
  bxm_line <- as.numeric(bxm_line)
  group.colors <- c(`B73 BC` = "#F8766D", `RIL` = "#00BA38", `Mo17 BC` ="#619CFF")
  plot <- ggplot(df, aes(.data[[trait]], colour = Population)) +
    geom_freqpoly(binwidth=bw, linewidth = 2, show.legend = FALSE) +
    geom_vline(xintercept = 1, color = "blue") +
    labs(x = trait_name,
         y = "Count") +
    theme(text = element_text(size=20),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_x_continuous(labels = scales::percent, limits = c(-.3,.4)) +
    scale_y_continuous(limits =  c(0,50)) +
    scale_color_manual(values=group.colors)
  return(plot)
}



plot_list <- lapply(list("SLB_WMD_BLUE_MPH_PCT", "GLS_WMD_BLUE_MPH_PCT", "DTA_BLUE_MPH_PCT", "PH_BLUE_MPH_PCT", "EH_BLUE_MPH_PCT"), function(x){ make_figure(blues, original_lines, x) })

figure <- plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]], ncol = 1)


make_figure_bph <- function(df, original_lines, trait){
  #trait_name <- case_when(
  #trait == "SLB_WMD_BLUE_MPH_PCT" ~ "Southern leaf blight resistance",
  #trait == "GLS_WMD_BLUE_MPH_PCT" ~ "Gray leaf spot resistance",
  #trait == "DTA_BLUE_MPH_PCT" ~ "Days to anthesis",
  # trait == "PH_BLUE_MPH_PCT" ~ "Plant height",
  # trait == "EH_BLUE_MPH_PCT" ~ "Ear height"
  #)
  trait_name <- case_when(
    trait == "SLB_WMD_BLUE_BPH_PCT" ~ "SLB",
    trait == "GLS_WMD_BLUE_BPH_PCT" ~ "GLS",
    trait == "DTA_BLUE_BPH_PCT" ~ "DTA",
    trait == "PH_BLUE_BPH_PCT" ~ "PH",
    trait == "EH_BLUE_BPH_PCT" ~ "EH"
  )
  bw <- 2 * IQR(blues[[trait]], na.rm=TRUE) / length(which(!is.na(blues[[trait]])))^(1/3)
  bxm_line <- original_lines %>%
    filter(Line == "B73 X Mo17") %>%
    select(!!trait)
  bxm_line <- as.numeric(bxm_line)
  group.colors <- c(`B73 BC` = "#F8766D", `RIL` = "#00BA38", `Mo17 BC` ="#619CFF")
  plot <- ggplot(blues, aes(.data[[trait]], colour = Population)) +
    geom_freqpoly(binwidth=bw, linewidth = 2, show.legend = FALSE) +
    geom_vline(xintercept = 1, color = "blue") +
    labs(x = trait_name,
         y = "Count") +
    theme(text = element_text(size=20),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_x_continuous(labels = scales::percent, limits = c(-.3,.4)) +
    scale_y_continuous(limits =  c(0,50)) +
    scale_color_manual(values=group.colors)
  return(plot)
}

bph_plot_list <- lapply(list("SLB_WMD_BLUE_BPH_PCT", "GLS_WMD_BLUE_BPH_PCT", "DTA_BLUE_BPH_PCT", "PH_BLUE_BPH_PCT", "EH_BLUE_BPH_PCT"), function(x){ make_figure_bph(blues, original_lines, x) })

figure_bph <- plot_grid(bph_plot_list[[1]], bph_plot_list[[2]], bph_plot_list[[3]], bph_plot_list[[4]], bph_plot_list[[5]], ncol = 1)

figure_r1 <- plot_grid(plot_list[[1]], bph_plot_list[[1]], ncol=2)
figure_r2 <- plot_grid(plot_list[[2]], bph_plot_list[[2]], ncol=2)
figure_r3 <- plot_grid(plot_list[[3]], bph_plot_list[[3]], ncol=2)
figure_r4 <- plot_grid(plot_list[[4]], bph_plot_list[[4]], ncol=2)
figure_r5 <- plot_grid(plot_list[[5]], bph_plot_list[[5]], ncol=2)
figure_combined <- plot_grid(figure_r1,figure_r2, figure_r3, figure_r4, figure_r5, labels = c("      SLB", "      GLS", "DTA", " PH", " EH"), ncol = 1,vjust=1.2,scale = .85, hjust=-1.6, label_x=c(-0.185,-0.185,-0.05,-0.045,-0.045))
pdf("figures/mph_bph_distributions_paper.pdf", height = 10)
figure_combined
dev.off()

figure <- plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]], labels = c('A', 'B', 'C', 'D','E'), label_size = 12, ncol = 1)


figure_combined <- plot_grid(figure,figure_bph,ncol=2)
pdf("figures/mph_bph_distributions_paper.pdf", height = 10)
figure_combined
dev.off()
