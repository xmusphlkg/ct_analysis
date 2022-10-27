
# packages ----------------------------------------------------------------

library(mgcv)
library(tidyverse)

# load data ---------------------------------------------------------------

set.seed(20220831)

load('../script/xiamen.RData')
source('./function.R')

# analysis Ct value -------------------------------------------------------

DataCtValue <- DataCtValue |>
  mutate(DateSeq = as.numeric(SampleDate - OnsetDate),
         VaccineDose = case_when(
           VaccineDose <= 1 ~ 'Unfully Vaccine',
           VaccineDose == 2 ~ 'Fully Vaccine',
           VaccineDose == 3 ~ 'Booster Dose',
           TRUE ~ as.character(VaccineDose)
         ))

DataCtValue |>
  ggplot(mapping = aes(x = DateSeq, y = CtORF1ab, color = CaseType))+
  geom_line(mapping = aes(group = CaseID),
            alpha = 0.1)+
  stat_smooth()

DataCtValue |>
  ggplot(mapping = aes(x = DateSeq, y = CtN, color = CaseType))+
  geom_line(mapping = aes(group = CaseID),
            alpha = 0.1)+
  geom_hline(yintercept = 35)+
  stat_smooth()+
  scale_x_continuous(breaks = seq(-2, 22, 2),
                     limits = c(-2, 22),
                     expand = c(0, 0))+
  theme_bw()

ggsave(filename = 'cttest.png', width = 5, height = 5)

DataCtValue |>
  ggplot(mapping = aes(x = DateSeq, y = CtORF1ab, color = VaccineDose))+
  geom_line(mapping = aes(group = CaseID),
            alpha = 0.1)+
  stat_smooth()

DataCtValue |>
  ggplot(mapping = aes(x = DateSeq, y = CtN, color = VaccineDose))+
  geom_line(mapping = aes(group = CaseID),
            alpha = 0.1)+
  stat_smooth(method = lm, formula = y ~ splines::bs(x, df = 6))

DataCtValue |>
  ggplot(mapping = aes(x = DateSeq, y = CtN, color = VaccineDose))+
  geom_line(mapping = aes(group = CaseID),
            alpha = 0.1)+
  stat_smooth(method = gam, formula = y ~ s(x))
