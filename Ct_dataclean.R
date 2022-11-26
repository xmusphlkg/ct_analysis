
# library -----------------------------------------------------------------

library(openxlsx)
library(tidyverse)

# data --------------------------------------------------------------------

datafile <- read.xlsx('../data/20220820.xlsx') |>
  mutate_at(vars(contains('date')), convertToDate) |>
  mutate(vaccine = as.character(vaccine)) |>
  filter(lineage == 'BA.2.76' & id != 3)
datafile_transmission <- datafile |>
  select(id, vaccine, ide_type,infector) |>
  left_join(datafile[,c('name', 'id', 'vaccine', 'ide_type')], by = c('infector' = 'name')) |>
  select(-infector)

# Ct Value ----------------------------------------------------------------

DataCtValue <- read.xlsx('../data/xiamen_omicron.xlsx')
DataVaccine <- read.xlsx('../data/casevaccine.xlsx')

DataVaccine <- DataVaccine |>
  mutate_at(vars(contains('Date')), convertToDateTime) |>
  left_join(datafile[,c('name', 'age', 'type', 'gender', 'date')], by = c(Name = 'name')) |>
  filter(!is.na(gender)) |>
  mutate(type = factor(type,
                       levels = c("轻型", "普通型"),
                       labels = c("Mild", "Moderate")),
         gender = factor(gender,
                         levels = 1:2,
                         labels = c('Male', 'Female')),
         VaccineType = str_remove_all(paste(FirstVaccineProduce,
                                            SecondVaccineProduce,
                                            ThirdVaccineProduce,
                                            sep = '_'),
                                      'NA'),
         VaccineType = str_remove_all(VaccineType, "_"),
         VaccineType = sapply(VaccineType, FUN = function(x){
           paste(sort(unique(strsplit(x, "")[[1]])), collapse = '')
         }),
         vaccine = if_else(VaccineType == 'C',
                           vaccine + 1,
                           vaccine),
         vaccine = factor(vaccine,
                          levels = 0:3,
                          labels = c('U', 'U', '2 dose', 'Booster dose')))

DataCtValue <- DataCtValue |>
  mutate(SampleDate = convertToDate(SampleDate),
         CtORF1ab = as.numeric(str_replace_all(CtORF1ab, "-", "40")),
         CtN = as.numeric(str_replace_all(CtN, "-", "40"))) |>
  select(CaseID, Name, Age, Gender, SampleDate, CtORF1ab, CtN) |>
  left_join(datafile[,c('name', 'type', 'date')], by = c(Name = 'name')) |>
  left_join(DataVaccine[,c('Name', 'vaccine', 'VaccineType')], by = c("Name")) |>
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
