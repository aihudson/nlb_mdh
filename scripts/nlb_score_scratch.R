data_1 <- read.csv("~/Downloads/2024-07-24-01-22-39_24NP_ForTablet_table.csv")
data_2 <- read.csv("~/Downloads/2024-08-01-07-11-09_24NP_plots_table.csv")
data_merge <- merge(data_1,data_2, by = "Plot.ID")
data_merge  <- data_merge %>% 
  mutate(NLB.080124.ah.score = as.numeric(gsub("[^0-9.]", "",NLB.080124.ah))) %>%
  filter(Population != "H100")
data_merge %>%
  group_by(Population) %>%
  summarise(cor = cor(DiseasePCTRating1,NLB.080124.ah.score,use = "complete.obs"))
ggplot(data = data_merge, mapping = aes(x=DiseasePCTRating1,y=NLB.080124.ah.score)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(vars(Population))
ggplot(data = data_merge, mapping = aes(x=DiseasePCTRating1,y=NLB.080124.ah.score)) +
  geom_point() +
  geom_smooth() +
  facet_grid(vars(Population))
model <- lm(NLB.080124.ah.score~DiseasePCTRating1, data = data_merge)

