
library(tidyverse)
library(openxlsx)
library(patchwork)

# load data ---------------------------------------------------------------

load('./xiamen.RData')

# panel A -----------------------------------------------------------------

datafile <-DataCtValue |>
  group_by(OnsetDate) |>
  summarise(n = length(unique(CaseID)),
            .groups = 'drop')

FigA <- ggplot(data = datafile)+
  geom_col(mapping = aes(x = as.numeric(OnsetDate - min(OnsetDate)), y = n),
           position = "stack",
           fill = 'grey', color = 'black',
           width = 1)+
  scale_x_continuous(breaks = seq(0, 40, 10),
                     limits = c(-5, 40),
                     expand = c(0, 0))+
  scale_y_continuous(breaks = seq(0, 15, 5),
                     limits = c(0, 15),
                     expand = c(0, 0))+
  theme_classic()+
  theme(panel.grid.major = element_line(color = 'grey'),
        axis.text.x = element_blank())+
  labs(y = 'Incidence',
       x = NULL,
       title = 'A')

# panel D -----------------------------------------------------------------

datafile <- DataCtValue |>
  pivot_longer(cols = c(CtN, CtORF1ab),
               names_to = 'Ct',
               values_to = 'value')

FigD <- ggplot(data = datafile)+
  stat_smooth(mapping = aes(x = as.numeric(SampleDate - OnsetDate),
                            y = value,
                            color = Ct,
                            fill = Ct),
              alpha = 0.3,
              method = lm,
              formula = y ~ splines::bs(x, df = 4))+
  geom_hline(yintercept = 30)+
  coord_cartesian(ylim = c(15, 40))+
  scale_x_continuous(breaks = seq(-5, 20, 5),
                     limits = c(-5, 20),
                     expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  scale_color_manual(values = c("#017098", "#987501"),
                     labels = c('N gene', 'ORF gene'))+
  scale_fill_manual(values = c("#017098", "#987501"))+
  theme_classic()+
  theme(panel.grid.major.x = element_line(color = 'grey'),
        axis.text.x = element_blank(),
        panel.grid.major = element_line(color = 'grey'),
        legend.position = c(1, 0.5),
        legend.justification = c(1, 1))+
  labs(y = 'Cycle Threshold Value',
       x = NULL,
       color = NULL,
       title = 'D')+
  guides(fill = 'none',
         color = guide_legend(override.aes = list(fill = 'white',
                                                  title = NULL)))

# panel B -----------------------------------------------------------------

Ct_complete <- function(i){
  # i = 3
  datafile <- DataCtValue |> filter(CaseID == unique(DataCtValue$CaseID)[i])
  # remove false negative of CtN
  LastPositiveDate <- datafile |>
    filter(CtN < 35)
  LastPositiveDate <- max(LastPositiveDate$SampleDate, na.rm = T)
  datafile <- datafile[!(datafile$SampleDate < LastPositiveDate & datafile$CtN >= 35),]
  # remove false negative of CtORF
  LastPositiveDate <- datafile |>
    filter(CtORF1ab < 35)
  LastPositiveDate <- max(LastPositiveDate$SampleDate, na.rm = T)
  datafile <- datafile[!(datafile$SampleDate < LastPositiveDate & datafile$CtORF1ab >= 35),]
  datafile$CaseID <- i
  return(datafile)
}
datafile <- lapply(1:length(unique(DataCtValue$CaseID)), Ct_complete)
datafile <- do.call('rbind', datafile)

library(paletteer)
colors <- c('#F39B7FFF', '#8491B4FF')

FigB <- ggplot(data = datafile)+
  geom_line(mapping = aes(x = as.numeric(SampleDate - min(datafile$OnsetDate, na.rm = T)),
                          y = CaseID,
                          group = CaseID,
                          color = CtN),
            size = 1)+
  scale_y_reverse()+
  scale_color_viridis_c(direction = -1,
                        option = 'D',
                        limits = c(10, 45))+
  # scale_color_gradientn(colors = colors,
  #                       limits = c(10, 45),
  #                       )+
  scale_x_continuous(breaks = seq(0, 40, 10),
                     limits = c(-5, 40),
                     expand = c(0, 0))+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.x = element_line(color = 'grey'),
        legend.position = c(1, 1),
        legend.justification = c(1, 1))+
  labs(y = 'Infections',
       x = NULL,
       color = 'N gene',
       title = 'B')

# panel C -----------------------------------------------------------------

FigC <- ggplot(data = datafile)+
  geom_line(mapping = aes(x = as.numeric(SampleDate - min(OnsetDate, na.rm = T)),
                          y = CaseID,
                          group = CaseID,
                          color = CtORF1ab),
            size = 1)+
  scale_y_reverse()+
  scale_color_viridis_c(direction = -1,
                        option = 'D',
                        limits = c(10, 45))+
  scale_x_continuous(breaks = seq(0, 40, 10),
                     limits = c(-5, 40),
                     expand = c(0, 0))+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.x = element_line(color = 'grey'),
        legend.position = c(1, 1),
        legend.justification = c(1, 1))+
  labs(y = 'Infections',
       x = 'Days since start of outbreak',
       color = 'ORF gene',
       title = 'C')

# panel E -----------------------------------------------------------------

datafile <- DataCtValue |>
  pivot_longer(cols = c(CtN, CtORF1ab),
               names_to = 'Ct',
               values_to = 'value')

datafile$VaccineGroup <- datafile$VaccineDose

colors <- c("#E64B35FF", "#3C5488FF", "#00A087FF")

FigE <- datafile |>
  filter(Ct == "CtN") |>
  ggplot()+
  stat_smooth(mapping = aes(x = as.numeric(SampleDate - OnsetDate),
                            y = value,
                            color = VaccineGroup,
                            fill = VaccineGroup),
              alpha = 0.3,
              method = lm,
              formula = y ~ splines::bs(x, df = 4))+
  geom_hline(yintercept = 30)+
  coord_cartesian(ylim = c(15, 40))+
  scale_x_continuous(breaks = seq(-5, 20, 5),
                     limits = c(-5, 20),
                     expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  scale_color_manual(values = colors,
                     labels = c('Unfully Vaccinated', 'Fully Vaccinated', 'Booster dose'))+
  scale_fill_manual(values = colors)+
  theme_classic()+
  theme(panel.grid.major.x = element_line(color = 'grey'),
        axis.text.x = element_blank(),
        panel.grid.major = element_line(color = 'grey'),
        legend.position = c(1, 0.5),
        legend.justification = c(1, 1))+
  labs(y = 'Cycle Threshold Value',
       x = NULL,
       color = NULL,
       title = 'E')+
  guides(fill = 'none',
         color = guide_legend(override.aes = list(fill = 'white',
                                                  title = NULL)))


FigF <- datafile |>
  filter(Ct == "CtORF1ab") |>
  ggplot()+
  stat_smooth(mapping = aes(x = as.numeric(SampleDate - OnsetDate),
                            y = value,
                            color = VaccineGroup,
                            fill = VaccineGroup),
              alpha = 0.3,
              method = lm,
              formula = y ~ splines::bs(x, df = 4))+
  geom_hline(yintercept = 30)+
  coord_cartesian(ylim = c(15, 40))+
  scale_x_continuous(breaks = seq(-5, 20, 5),
                     limits = c(-5, 20),
                     expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  scale_color_manual(values = colors,
                     labels = c('Unfully Vaccinated', 'Fully Vaccinated', 'Booster dose'))+
  scale_fill_manual(values = colors)+
  theme_classic()+
  theme(panel.grid.major.x = element_line(color = 'grey'),
        panel.grid.major = element_line(color = 'grey'),
        legend.position = c(1, 0.5),
        legend.justification = c(1, 1))+
  labs(y = 'Cycle Threshold Value',
       x = 'Days from Onset Date',
       color = NULL,
       title = 'F')+
  guides(fill = 'none',
         color = guide_legend(override.aes = list(fill = 'white',
                                                  title = NULL)))

FigA + FigB + FigC + FigD + FigE + FigF+
  plot_layout(ncol = 2, byrow = F)

ggsave('./test.pdf', width = 10, height = 8, device = cairo_pdf)

# bs predict --------------------------------------------------------------

testvalue <- seq(-3, 12, 0.05)

spN <- datafile |>
  filter(Ct == "CtN") |>
  mutate(x = as.numeric(SampleDate - OnsetDate))
spN <- lm(value ~ splines::bs(x, df = 4), data = spN)
summary(spN)

outcome <- predict(spN, data.frame(x = testvalue), interval = "predict") |>
  as.data.frame() |>
  mutate(x = testvalue,
         fit = abs(fit -30),
         lwr = abs(lwr -30),
         upr = abs(upr - 30))

spO <- datafile |>
  filter(Ct == "CtORF1ab") |>
  mutate(x = as.numeric(SampleDate - OnsetDate))
spO <- lm(value ~ splines::bs(x, df = 4), data = spO)
summary(spO)

outcome <- predict(spO, data.frame(x = testvalue), interval = "predict") |>
  as.data.frame() |>
  mutate(x = testvalue,
         fit = abs(fit -30),
         lwr = abs(lwr -30),
         upr = abs(upr - 30)) |>
  arrange(fit)


for (i in levels(datafile$VaccineGroup)) {
  spN <- datafile |>
    filter(Ct == "CtN" & VaccineGroup == i) |>
    mutate(x = as.numeric(SampleDate - OnsetDate))
  spN <- lm(value ~ splines::bs(x, df = 4), data = spN)
  print(i)
  print(summary(spN))
  outcome <- predict(spN, data.frame(x = testvalue), interval = "predict") |>
    as.data.frame() |>
    mutate(x = testvalue,
           fit = abs(fit -30),
           lwr = abs(lwr -30),
           upr = abs(upr - 30)) |>
    arrange(fit)
  print(head(outcome))
}

for (i in levels(datafile$VaccineGroup)) {
  spO <- datafile |>
    filter(Ct == "CtORF1ab" & VaccineGroup == i) |>
    mutate(x = as.numeric(SampleDate - OnsetDate))
  spO <- lm(value ~ splines::bs(x, df = 4), data = spO)
  print(i)
  print(summary(spO))
  outcome <- predict(spO, data.frame(x = testvalue), interval = "predict") |>
    as.data.frame() |>
    mutate(x = testvalue,
           fit = abs(fit -30),
           lwr = abs(lwr -30),
           upr = abs(upr - 30)) |>
    arrange(fit)
  print(head(outcome))
}

