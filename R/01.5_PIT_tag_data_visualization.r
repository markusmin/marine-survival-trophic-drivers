# 01.5 PIT tag data visualizations

# Description: This script takes the output of 01_PIT_tag_data.Rmd and creates visualizations
# for the supplemental figures for Min et al. 2026.

library(tidyverse)
library(here)

# load the SAR data - the output of 01_PIT_tag_data.Rmd
SAR_data <- read.csv(here::here("model_inputs", "chinook_det_hist.csv"))


# inspect distribution of outmigration timing

SAR_data %>% 
  mutate(BON_juv_det_date = substr(BON_juv_det_time, 1, 10)) %>% 
  mutate(BON_juv_det_date = ymd(BON_juv_det_date)) %>% 
  mutate(outmigration_date = yday(BON_juv_det_date)) -> SAR_data

# plot distribution of outmigration timing
outmigration_timing_plot <- ggplot(SAR_data, aes(x = outmigration_date)) +
  geom_histogram() +
  facet_wrap(~group + run_name)

ggsave(here::here("figures", "outmigration_timing_plot_by_group.png"), outmigration_timing_plot,
       height = 12, width = 12)


### SUPPLEMENTAL RESULTS FIGURE ###

# comparison of run timing for Snake River Fall and Upper Columbia Summer/Fall
SAR_data %>% 
  mutate(stock = paste0(group, " ", run_name)) -> SAR_data
subset(SAR_data, stock %in% c("Upper Columbia Summer", "Upper Columbia Fall", "Snake River Fall")) -> SAR_data_UCSF_SRF
SAR_data_UCSF_SRF %>% 
  mutate(paper_group = ifelse(stock %in% c("Upper Columbia Summer", "Upper Columbia Fall"), "Upper Columbia Summer/Fall", stock)) -> SAR_data_UCSF_SRF

SRF_UCSF_outmigration_timing_plot <- ggplot(SAR_data_UCSF_SRF, aes(x = outmigration_date)) +
  geom_histogram(bins = 60) +
  geom_vline(xintercept = 167, lty = 2, color = "red") +
  facet_wrap(~paper_group, ncol = 1) +
  xlab("Outmigration Date (Julian Day)") +
  ylab("Count") +
  theme(panel.grid.major = element_line(color = "gray90"),
        panel.background = element_rect(fill = "white", color = NA),
        panel.border = element_rect(color = NA, fill=NA, linewidth=0.4),
        legend.key.height = unit(1.25, "cm"),
        legend.key.width = unit(1.25, "cm"),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 20),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        strip.background = element_rect(fill = "white"),
        axis.title.x = element_text(size = 20, margin = margin(t = 10)),
        axis.title.y = element_text(size = 20, margin = margin(r = 10)),
        # these plot margins are to leave space for the population name on the big figure
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2),"cm"))

ggsave(here::here("figures", "paper_figures", "outmigration_timing_plot_SRF_UCSF.png"), SRF_UCSF_outmigration_timing_plot,
       height = 12, width = 12)




# inspect distribution of years with adult returns
subset(SAR_data, adult_det == 1) -> SAR_adult_det

adult_returns_plot <- ggplot(SAR_adult_det, aes(x = run_year)) +
  geom_histogram(breaks = seq(1998, 2025, 1)) +
  facet_wrap(~group + run_name)

ggsave(here::here("figures", "adult_returns_plot.png"), adult_returns_plot,
       height = 12, width = 12)

SAR_adult_det %>% 
  group_by(group, run_name, run_year) %>% 
  summarise(N = sum(adult_det)) -> adult_det_by_year

# inspect distribution of juvenile run years

adult_returns_plot <- ggplot(SAR_data, aes(x = run_year)) +
  geom_histogram(breaks = seq(1998, 2025, 1)) +
  facet_wrap(~group + run_name)

ggsave(here::here("figures", "annual_juveniles_BON_det_plot.png"), adult_returns_plot,
       height = 12, width = 12)


#### Create visualizations for our two target groups: Upper Columbia and Snake River Fall Chinook ####

### show outmigration timing, with June 15 clearly demarcated

SRF_outmigration_plot <- ggplot(subset(SAR_data, run_name == "Fall" & group == "Snake River"), aes(x = outmigration_date)) +
  geom_histogram() +
  geom_vline(xintercept = 167, lty = 2, color = "red") +
  annotate(geom = "text", x = 168, y = 150000, label = "June 15", hjust = 0, color = "red") +
  xlab("Outmigration Date (Julian Day)") +
  ylab("Count") +
  ggtitle("Snake River Fall Chinook")

ggsave(here::here("figures", "PIT_tag_data", "SRF_outmigration_timing_plot.png"), SRF_outmigration_plot,
       height = 8, width = 8)

UCF_outmigration_plot <- ggplot(subset(SAR_data, run_name == "Fall" & group == "Upper Columbia"), aes(x = outmigration_date)) +
  geom_histogram() +
  geom_vline(xintercept = 167, lty = 2, color = "red") +
  annotate(geom = "text", x = 168, y = 5000, label = "June 15", hjust = 0, color = "red") +
  xlab("Outmigration Date (Julian Day)") +
  ylab("Count") +
  ggtitle("Upper Columbia Fall Chinook")

ggsave(here::here("figures", "PIT_tag_data", "UCF_outmigration_timing_plot.png"), UCF_outmigration_plot,
       height = 8, width = 8)

UCS_outmigration_plot <- ggplot(subset(SAR_data, run_name == "Summer" & group == "Upper Columbia"), aes(x = outmigration_date)) +
  geom_histogram() +
  geom_vline(xintercept = 167, lty = 2, color = "red") +
  annotate(geom = "text", x = 168, y = 23000, label = "June 15", hjust = 0, color = "red") +
  xlab("Outmigration Date (Julian Day)") +
  ylab("Count") +
  ggtitle("Upper Columbia Summer Chinook")

ggsave(here::here("figures", "PIT_tag_data", "UCS_outmigration_timing_plot.png"), UCS_outmigration_plot,
       height = 8, width = 8)

### show data availability by year - juvenile and adult counts
detection_colors <- c("Adult" = "black",
                      "Juvenile" = "black")

detection_fills <- c("Adult" = "black",
                      "Juvenile" = "white")

SAR_data %>% 
  group_by(group, run_name, run_year) %>% 
  summarise(N_adult = sum(adult_det),
            N_juv = sum(juv_det)) %>% 
  pivot_longer(cols = c(N_adult, N_juv)) %>% 
  mutate(name = ifelse(name == "N_adult", "Adult", "Juvenile")) %>% 
  dplyr::rename(Detections = name) -> SAR_det_by_year

UCS_detections_plot <- ggplot(subset(SAR_det_by_year, run_name == "Summer" & group == "Upper Columbia"), 
                              aes(x = run_year, y = value, color = Detections, fill = Detections)) +
  scale_fill_manual(values = detection_fills) +
  scale_color_manual(values = detection_colors) +
  geom_bar(stat = "identity", position = "dodge2") +
  scale_x_continuous(expand = c(0,0), limits = c(1997, 2026)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 18000)) +
  ylab("Detections") +
  xlab("Run Year") +
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill="white", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.position = c(0.8, 0.8),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))

ggsave(here::here("figures", "PIT_tag_data", "UCS_detections_plot.png"), UCS_detections_plot,
       height = 8, width = 12)

SRF_detections_plot <- ggplot(subset(SAR_det_by_year, run_name == "Fall" & group == "Snake River"), 
                              aes(x = run_year, y = value, color = Detections, fill = Detections)) +
  scale_fill_manual(values = detection_fills) +
  scale_color_manual(values = detection_colors) +
  geom_bar(stat = "identity", position = "dodge2") +
  scale_x_continuous(expand = c(0,0), limits = c(1997, 2026)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 85000)) +
  ylab("Detections") +
  xlab("Run Year") +
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill="white", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.position = c(0.8, 0.8),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))

ggsave(here::here("figures", "PIT_tag_data", "SRF_detections_plot.png"), SRF_detections_plot,
       height = 8, width = 12)




#### Make table 1 for paper ####

# load the SAR data
SAR_data <- read.csv(here::here("model_inputs", "chinook_det_hist.csv"))

# inspect distribution of outmigration timing

SAR_data %>% 
  mutate(BON_juv_det_date = substr(BON_juv_det_time, 1, 10)) %>% 
  mutate(BON_juv_det_date = ymd(BON_juv_det_date)) %>% 
  mutate(outmigration_date = yday(BON_juv_det_date)) -> SAR_data

#### Snake River Fall Chinook

# extract Snake River Fall Chinook from SAR data
SRF <- subset(SAR_data, run_name == "Fall" & group == "Snake River")

#### subset by outmigration date
hist(SRF$outmigration_date)
summary(SRF$outmigration_date)

# restrict the analysis to fish that could have theoretically been caught in the June JSOES survey.
# we will call this April 15 - June 15 (for now!)
# day 105 is April 15 in a non-leap year and day 167 is June 15 in a leap year
# SRF_subset <- filter(SRF, outmigration_date >= 105 & outmigration_date <= 167)

# let's just have a final cutoff - based on the fact that fish can hang out
# off the coast for a while, it's less justified to have a front-end cutoff date
SRF_subset <- filter(SRF, outmigration_date <= 167)
nrow(SRF_subset)/nrow(SRF)
# this leaves us with 54% of our initial dataset


#### Join with the transport data

transport_data_export <- read.csv(here::here("model_inputs", "chinook_transport.csv"))

# join the juvenile transport data with the SAR data
SRF_subset %>% 
  left_join(transport_data_export, join_by(tag_code == tag_id)) -> SRF_subset

SRF_subset$transport[is.na(SRF_subset$transport)] <- 0

# inspect number of transported fish
table(SRF_subset$transport)


#### Drop the few natural origin fish

# fish with rear type U or H are hatchery fish
SRF_subset %>% 
  mutate(rear_numeric = ifelse(rear_type_code %in% c("H", "U"), 1, 0)) -> SRF_subset

# ok, so they're all hatchery. That's perhaps not surprising, given what we see in the JSOES survey is basically all hatchery
# Let's only keep hatchery fish
SRF_subset <- subset(SRF_subset, rear_type_code != "W")


#### Upper Columbia Summer/Fall Chinook

# extract Snake River Fall Chinook from SAR data
UCSF <- subset(SAR_data, run_name %in% c("Summer", "Fall") & group == "Upper Columbia")

#### subset by outmigration date
hist(UCSF$outmigration_date)
summary(UCSF$outmigration_date)

# restrict the analysis to fish that could have theoretically been caught in the June JSOES survey.
# we will call this April 15 - June 15 (for now!)
# day 105 is April 15 in a non-leap year and day 167 is June 15 in a leap year
# UCSF_subset <- filter(UCSF, outmigration_date >= 105 & outmigration_date <= 167)

# let's just have a final cutoff - based on the fact that fish can hang out
# off the coast for a while, it's less justified to have a front-end cutoff date
UCSF_subset <- filter(UCSF, outmigration_date <= 167)
nrow(UCSF_subset)/nrow(UCSF)
# this leaves us with 72% of our initial dataset

#### Drop the few natural origin fish (28 of ~97,500)

table(UCSF_subset$rear_type_code)

# fish with rear type U or H are hatchery fish
UCSF_subset %>% 
  mutate(rear_numeric = ifelse(rear_type_code %in% c("H", "U"), 1, 0)) -> UCSF_subset

# ok, so they're all hatchery. That's perhaps not surprising, given what we see in the JSOES survey is basically all hatchery
# Let's only keep hatchery fish
UCSF_subset <- subset(UCSF_subset, rear_type_code != "W")


#### Reformat for table 1 for the paper

SRF_subset %>% 
  bind_rows(UCSF_subset) -> SRF_UCSF_subset

as.data.frame(table(SRF_subset$run_year)) %>% 
  dplyr::rename(run_year = Var1, juv_det = Freq) %>% 
  mutate(stock_group = "Snake River Fall") -> SRF_subset_juv_table

as.data.frame(table(subset(SRF_subset, adult_det == 1)$run_year)) %>% 
  dplyr::rename(run_year = Var1, adult_det = Freq) %>% 
  mutate(stock_group = "Snake River Fall") -> SRF_subset_adult_table

SRF_subset_juv_table %>% 
  left_join(SRF_subset_adult_table, by = c("run_year", "stock_group")) %>% 
  mutate(adult_det = ifelse(is.na(adult_det), 0, adult_det)) -> SRF_subset_table

as.data.frame(table(UCSF_subset$run_year)) %>% 
  dplyr::rename(run_year = Var1, juv_det = Freq) %>% 
  mutate(stock_group = "Upper Columbia Summer/Fall") -> UCSF_subset_juv_table

as.data.frame(table(subset(UCSF_subset, adult_det == 1)$run_year)) %>% 
  dplyr::rename(run_year = Var1, adult_det = Freq) %>% 
  mutate(stock_group = "Upper Columbia Summer/Fall") -> UCSF_subset_adult_table

UCSF_subset_juv_table %>% 
  left_join(UCSF_subset_adult_table, by = c("run_year", "stock_group")) %>% 
  mutate(adult_det = ifelse(is.na(adult_det), 0, adult_det)) -> UCSF_subset_table

SRF_subset_table %>% 
  mutate(SAR = round(adult_det/juv_det,3)) %>% 
  mutate("Snake River Fall Run" = paste0(adult_det, "/", juv_det, " (", SAR, ")")) -> SRF_subset_table

UCSF_subset_table %>% 
  mutate(SAR = round(adult_det/juv_det,3)) %>% 
  mutate("Upper Columbia Summer/Fall Run" = paste0(adult_det, "/", juv_det, " (", SAR, ")")) -> UCSF_subset_table

# join them all together
UCSF_subset_table %>% 
  dplyr::select(run_year, "Upper Columbia Summer/Fall Run") -> UCSF_for_join

SRF_subset_table %>% 
  dplyr::select(run_year, "Snake River Fall Run") -> SRF_for_join

SRF_for_join %>% 
  left_join(UCSF_for_join, by = "run_year") %>% 
  mutate(run_year = as.numeric(as.character(run_year))) %>% 
  filter(run_year >= 2001 & run_year <= 2021) -> PIT_tag_table

write.csv(PIT_tag_table, here::here("figures", "paper_figures", "SAR_data_table.csv"), row.names = FALSE)

# second version of this table
SRF_subset_table %>% 
  dplyr::rename(SR_returning_adults = adult_det, 
                SR_outmigrants = juv_det,
                SR_SAR = SAR) %>% 
  dplyr::select(run_year, SR_outmigrants, SR_returning_adults, SR_SAR) -> SRF_paper_table

UCSF_subset_table %>% 
  dplyr::rename(UC_returning_adults = adult_det, 
                UC_outmigrants = juv_det,
                UC_SAR = SAR) %>% 
  dplyr::select(run_year, UC_outmigrants, UC_returning_adults, UC_SAR) -> UCSF_paper_table

SRF_paper_table %>% 
  left_join(UCSF_paper_table, by = "run_year") %>% 
  mutate(run_year = as.numeric(as.character(run_year))) %>% 
  filter(run_year >= 2001 & run_year <= 2021) -> PIT_tag_table_2

write.csv(PIT_tag_table_2, here::here("figures", "paper_figures", "SAR_data_table_v2.csv"), row.names = FALSE)

