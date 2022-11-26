


Data <- DataValue |>
  mutate(LastVaccine = max(first))
