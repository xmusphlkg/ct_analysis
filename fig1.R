
# packages ----------------------------------------------------------------

library(tidyverse)
library(openxlsx)
library(patchwork)
library(ggraph)
library(igraph, include.only = 'graph_from_data_frame')
library(ggsci)
library(MASS, include.only = 'fitdistr')
library(epitrix)
library(EpiEstim)

# load data ---------------------------------------------------------------

set.seed(20220831)

load('../script/xiamen.RData')
source('./function.R')

datafile_info <- datafile |>
  select(id, type, vaccine, control) |>
  rename('type' = 'type',
         'vaccine' = 'vaccine') |>
  mutate(label = id,
         vaccine = factor(vaccine,
                          levels = as.character(3:0),
                          labels = c('A', 'B', 'C', 'C')))
fill_color <- pal_nejm()(4)

# network plot ------------------------------------------------------------

fig <- graph_from_data_frame(d = datafile_transmission[!is.na(datafile_transmission$id.y),c('id.y', 'id.x')],
                             vertices = datafile_info,
                             directed = T)
fig_a <- ggraph(fig,layout = "kk")+
  geom_edge_link(arrow = arrow(length = unit(1.5, 'mm')),
                 end_cap = circle(4, 'mm'),
                 check_overlap = F,
                 width = 0.7,
                 show.legend = F)+
  geom_node_circle(aes(fill = control,
                       linetype = vaccine,
                       r = 0.35),
                   size = 1,
                   show.legend = T)+
  geom_node_text(aes(label = label),
                 colour = 'white',
                 show.legend = F) +
  coord_equal(clip = 'off')+
  scale_fill_manual(values = fill_color[2:1])+
  scale_linetype_manual(name = 'Vaccination status',
                        labels = c('Booster Dose', 'Fully Vaccinated', 'Unfully Vaccinated'),
                        values = c('solid', 'longdash', 'dotted'))+
  theme_graph()+
  theme(
    legend.position = 'none',
    # plot.title = element_text(size = 16, hjust = 0, vjust = 0, face = 'bold'),
    plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"),
    # panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 12, hjust = .5, vjust = 0, face = 'bold')
  )+
  labs(title = 'A')

library(Cairo)
ggsave('../outcome/fig1_a.pdf', fig_a, width = 6, height = 6, device = cairo_pdf)

# epicurve ----------------------------------------------------------------

datafile$lineage[is.na(datafile$lineage)] <- 'Unknow'

order <- c("BA.2.76", "BA.2.3.7", "BA.5.2",  "Unknow")

fig_b <- ggplot(data = filter(datafile, lineage == order[1]))+
  geom_col(mapping = aes(x = date, y = 1, fill = factor(lineage)),
           color = 'white',
           width = 1,
           position = position_stack(reverse = TRUE))+
  scale_fill_manual(values = fill_color)+
  scale_y_continuous(expand = expansion(add = c(0, 3)),
                     breaks = seq(0, 16, 4))+
  scale_x_date(date_labels = "%m/%d",
               date_breaks = "2 day",
               expand = expansion(add = c(2, 1)))+
  coord_equal(ratio = 1)+
  theme_classic()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text = element_text(color = 'black'),
        # legend.position = c(0.2, 0.8),
        legend.position = 'none',
        legend.justification = c(0, 1),
        legend.background = element_rect(color = 'black', size = 1),
        # legend.margin = margin(),
        plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"))+
  labs(x = "Onset date",
       y = "Number of infections",
       title = 'B',
       fill = "Lineage\nof VOC")

# ggsave('../outcome/fig1_b.pdf', width = 5, height = 4, device = cairo_pdf)

# serial interval ---------------------------------------------------------

infector_276 <- datafile |>
  filter(lineage == 'BA.2.76') |>
  select(infector) |>
  unique() |>
  left_join(datafile[,c('name', 'date', 'lineage')], by = c(infector = 'name')) |>
  filter(lineage == 'BA.2.76')

datafile_serial <- datafile |>
  filter(infector %in% infector_276$infector) |>
  select(name, infector, date) |>
  left_join(infector_276[,c("infector", "date")], by = 'infector') |>
  mutate(date_seq = as.numeric(date.x - date.y))

si_dist <- fit_best(datafile_serial$date_seq)

# ggplot(data = data.frame(x = c(0, 7)), aes(x))+
#   geom_histogram(data = datafile_serial,
#                  mapping = aes(x = date_seq,
#                                y = ..density..),
#                  binwidth = 1,
#                  color = 'grey',
#                  alpha = 0.7) +
#   stat_function(fun = dgamma, n = 100,
#                 args = list(shape = si_dist$shape,
#                             scale = 1/si_dist$rate),
#                 color = 'red')+
#   scale_y_continuous(limits = c(0, 0.5),
#                      expand = c(0, 0))+
#   scale_x_continuous(breaks = 0:7)+
#   theme_classic()+
#   theme(panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank(),
#         axis.text = element_text(color = 'black'),
#         plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"))+
#   labs(x = "Serial Interval",
#        y = "Relative frequency",
#        title = 'C')

# incubation --------------------------------------------------------------

datafile_incubat <- datafile |>
  filter(lineage == 'BA.2.76') |>
  filter(infector %in% infector_276$infector) |>
  select(name, infector, exposedate1, exposedate2, date) |>
  rename(dateonset = 'date') |>
  mutate(exposedate1 = if_else(is.na(exposedate1),
                               exposedate2,
                               exposedate1),
         date_seq_1 = dateonset - exposedate1,
         date_seq_2 = dateonset - exposedate2) |>
  filter(!is.na(exposedate1))
datafile_incubat$expose_date <- mapply(expose_date, datafile_incubat$exposedate1, datafile_incubat$exposedate2)

ib_dist <- fit_gamma_incubation_dist(datafile_incubat,
                                     dateonset,
                                     exposedate1,
                                     exposedate2)

fig_c <- ggplot(data = data.frame(x = c(0, 7)), aes(x))+
  geom_histogram(data = datafile_serial,
                 mapping = aes(x = date_seq,
                               y = ..density..,
                               fill = 'Serial Interval'),
                 binwidth = 1,
                 color = 'white',
                 alpha = 0.9) +
  stat_function(fun = dgamma, n = 100,
                args = list(shape = si_dist$shape,
                            scale = 1/si_dist$rate),
                mapping = aes(color = 'Fitted Serial Interval'))+
  stat_function(fun = dgamma, n = 100,
                args = list(shape = ib_dist$distribution$parameters$shape,
                            scale = ib_dist$distribution$parameters$scale),
                mapping = aes(color = 'Fitted Incubation Period'))+
  scale_y_continuous(limits = c(0, 0.5),
                     expand = c(0, 0))+
  scale_x_continuous(breaks = 0:7)+
  scale_color_manual(values = fill_color)+
  scale_fill_manual(values = 'grey')+
  theme_classic()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = c(1, 1),
        legend.justification = c(1, 1),
        axis.text = element_text(color = 'black'),
        plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"))+
  labs(x = "Days",
       y = "Relative frequency",
       title = 'C',
       colour = NULL,
       fill = NULL)

# Time-varing reproductive number -----------------------------------------

datafile_rt <- datafile |>
  filter(lineage == 'BA.2.76') |>
  group_by(date) |>
  count() |>
  as.data.frame() |>
  complete(date = seq.Date(min(date), max(date)+1, by = 'day'),
           fill = list(n = 0))

start_dates <- seq(2, nrow(datafile_rt) - 4)
end_dates <- start_dates + 4

config_lit <- make_config(list(
  mean_si = si_dist$mean,
  std_si = si_dist$sd,
  t_start = start_dates,
  t_end = end_dates
))
names(datafile_rt) <- c('dates', 'I')
date_st <- min(datafile_rt$dates)

epiestim_res_lit <- estimate_R(incid = datafile_rt,
                               method = "parametric_si",
                               config = config_lit)
outcome <- epiestim_res_lit$R
outcome$date <- (outcome$t_start + outcome$t_end) / 2 + date_st
outcome$t_start <- date_st + outcome$t_start
outcome$t_end <- date_st + outcome$t_end

fig_d <- ggplot(data = outcome,
       mapping = aes(x = date))+
  geom_ribbon(
    mapping = aes(
      ymin = `Quantile.0.025(R)`,
      ymax = `Quantile.0.975(R)`
    ),
    fill = fill_color[1],
    alpha = 0.3,
    show.legend = F
  ) +
  geom_line(
    mapping = aes(y = `Mean(R)`),
    color = fill_color[1],
    show.legend = F
  ) +
  geom_line(
    mapping = aes(y = 1),
    color = 'black',
    linetype = "dashed",
    show.legend = F
  )+
  scale_y_continuous(expand = expansion(add = c(0, 0)),
                     breaks = seq(0, 4, 1),
                     limits = c(0, 4)
                     )+
  scale_x_date(date_labels = "%m/%d",
               date_breaks = "2 day",
               expand = expansion(add = c(0, 0)))+
  theme_classic()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text = element_text(color = 'black'),
        legend.position = 'none',
        legend.justification = c(0, 1),
        legend.background = element_rect(color = 'black', size = 1),
        # legend.margin = margin(),
        plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"))+
  labs(x = "Onset date",
       y = "Time-varing reproductive number",
       title = 'D')

fig_a + fig_b + fig_c + fig_d&
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        text = element_text(size = 12),
        title = element_text(size = 16))

ggsave('../outcome/fig.pdf', width = 10, height = 10, device = cairo_pdf)
