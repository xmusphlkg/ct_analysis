
# library -----------------------------------------------------------------

library(openxlsx)
library(tidyverse)

# data --------------------------------------------------------------------

datafile <- read.xlsx('../data/20220820.xlsx') |>
  mutate_at(vars(contains('date')), convertToDate) |>
  mutate(vaccine = as.character(vaccine)) |>
  filter(lineage == 'BA.2.76')
datafile_transmission <- datafile |>
  select(id, vaccine, ide_type,infector) |>
  left_join(datafile[,c('name', 'id', 'vaccine', 'ide_type')], by = c('infector' = 'name')) |>
  select(-infector)

# Ct Value ----------------------------------------------------------------

DataCtValue <- read.xlsx('../data/xiamen_omicron.xlsx')
DataCtValue <- DataCtValue |>
  mutate(SampleDate = convertToDate(SampleDate),
         CtORF1ab = as.numeric(str_replace_all(CtORF1ab, "-", "40")),
         CtN = as.numeric(str_replace_all(CtN, "-", "40"))) |>
  select(CaseID, Name, Age, Gender, SampleDate, CtORF1ab, CtN) |>
  left_join(datafile[,c('name', 'vaccine', 'type', 'date')], by = c(Name = 'name')) |>
  rename(c(VaccineDose = 'vaccine',
           CaseType = 'type',
           OnsetDate = 'date')) |>
  mutate(CaseType = case_when(
    CaseType == '轻型' ~ 'Mild',
    CaseType == '普通型' ~ 'Moderate',
    TRUE ~ as.character(CaseType)
  )) |>
  filter(!is.na(SampleDate) & !is.na(CaseType)) |>
  select(-Name)

save.image('./xiamen.RData')
