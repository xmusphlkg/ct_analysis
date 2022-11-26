
# packages ----------------------------------------------------------------

library(mgcv)
library(tidyverse)

# load data ---------------------------------------------------------------

set.seed(20220831)
load('./omicron.RData')

datafile$VaccineGroup <- factor(datafile$VaccineDose,
                                levels = as.character(0:3),
                                labels = c('B', 'B', 'B', 'A'))
datafile$Ct <- factor(datafile$Ct,
                      levels = c('CtN', 'CtORF1ab'),
                      labels = c('N gene', 'ORF gene'))
colors <- ggsci::pal_npg()(2)

# figA --------------------------------------------------------------------

FigA <- ggplot(data = datafile)+
  stat_smooth(mapping = aes(x = as.numeric(SampleDate - OnsetDate),
                            y = value,
                            color = VaccineGroup,
                            fill = VaccineGroup),
              alpha = 0.3,
              method = lm,
              formula = y ~ splines::bs(x, df = 4))+
  geom_hline(yintercept = 30)+
  coord_cartesian(ylim = c(15, 40))+
  scale_x_continuous(breaks = seq(0, 30, 10),
                     limits = c(-5, 30),
                     expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  scale_color_manual(values = colors,
                     labels = c('Booster', 'Unbooster'))+
  scale_fill_manual(values = colors)+
  theme_classic()+
  theme(panel.grid.major.x = element_line(color = 'grey'),
        panel.grid.major = element_line(color = 'grey'),
        legend.position = c(1, 0.5),
        legend.justification = c(1, 1))+
  facet_wrap(.~Ct,
             scales = 'free')+
  labs(y = 'Cycle Threshold Value',
       x = 'Days from Onset Date',
       color = NULL,
       title = 'A')+
  guides(fill = 'none',
         color = guide_legend(override.aes = list(fill = 'white',
                                                  title = NULL)))

FigB <- ggplot(data = datafile)+
  stat_smooth(mapping = aes(x = as.numeric(SampleDate - OnsetDate),
                            y = value,
                            color = CaseType,
                            fill = CaseType),
              alpha = 0.3,
              method = lm,
              formula = y ~ splines::bs(x, df = 4))+
  geom_hline(yintercept = 30)+
  coord_cartesian(ylim = c(15, 40))+
  scale_x_continuous(breaks = seq(0, 30, 10),
                     limits = c(-5, 30),
                     expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  scale_color_manual(values = colors)+
  scale_fill_manual(values = colors)+
  theme_classic()+
  theme(panel.grid.major.x = element_line(color = 'grey'),
        panel.grid.major = element_line(color = 'grey'),
        legend.position = c(1, 0.5),
        legend.justification = c(1, 1))+
  facet_wrap(.~Ct,
             scales = 'free')+
  labs(y = 'Cycle Threshold Value',
       x = 'Days from Onset Date',
       color = NULL,
       title = 'B')+
  guides(fill = 'none',
         color = guide_legend(override.aes = list(fill = 'white',
                                                  title = NULL)))

FigA + FigB +
  plot_layout(ncol = 1)

ggsave('./test2.pdf', width = 8, height = 6, device = cairo_pdf)
