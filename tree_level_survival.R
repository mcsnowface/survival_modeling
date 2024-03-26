
# CLEAN UP
rm(list=ls())

# packages
library(MASS)
library(tidyverse)
library(ggeffects)
library(lme4)
library(patchwork)

# read the MPB PSP selection with attack year
# the file was originally called "PSPs_coord.csv", I renamed it for myself to find it faster
psps_attack <- read_csv("PSPs_attack_year.csv") %>% 
  # make plot number a character which is needed for the join below
  mutate(company_plot_number = as.character(company_plot_number)) 

# Mountain pine beetle PSPs
mpb_psps <- psps_attack$company_plot_number

# read PYGI tree level master file and immediately filter for only MPB PSPs
# I didn't include the file in the GitHub repository due to its size. I can be found on the Growth and Yield Google drive
# this runs on the master file I created and stored in the PGYI (R edition) folder.
read_csv("G:\\Shared drives\\Growth & Yield Lab\\Data Sets\\PGYI (R edition)\\pgyi_tree_level_master.csv") %>% 
  filter(company_plot_number %in% mpb_psps) -> dbh_mpb_psps

# join attack year file
dbh_mpb_psps %>% 
  left_join(psps_attack, by = c("company_plot_number")) -> dbh_mpb_psps_attack

# check whether number of PSP matches with 236
length(unique(dbh_mpb_psps_attack$company_plot_number))


# data preparation and cleaning
dbh_mpb_psps_attack %>%
  # remove any dbh measurements with NAs
  # at present I don't have a comparable way of imputing dbh measurements they way the SAS code does!
  filter(!is.na(dbh)) %>%
  # filter for alive and dead trees and exclude border trees
  filter(condition_code1 %in% c(0:12),
         tree_type != "B",!is.na(dbh)) %>%
  # create a new survival column and create grouping variables for species
  mutate(
    survival = if_else(condition_code1 %in% c(1, 2), 0, 1),
    species_groups = case_when(
      species %in% c("Aw", "Pb", "Bw", "Ax") ~ "Deciduous",
      species %in% c("Pl", "Pw", "Pf", "Pj", "Px") ~ "Pine",
      species %in% c("Sw", "Se", "Sb", "Sx") ~ "Spruce",
      species %in% c("Fb", "Fa", "Fd") ~ "Fir",
      species %in% c("Lt", "Lw", "La", "Ls") ~ "Larch",
      .default = "Others"
    ),
    # balt_groups = case_when(
    #   species %in% c("Fb", "Fa", "Fd", "Sw", "Se", "Sb", "Sx") ~ "spruce_fir",
    #   species %in% c("Pl", "Pw", "Pf", "Pj", "Px") ~ "pine",
    #   species %in% c("Aw", "Pb", "Bw", "Ax") ~ "deciduous"
    # ),
    # create some more variables that may be needed later in the analysis
    measurement_interval = measurement_year - establishment_year,
    tree_id = paste(company_plot_number, tree_number, sep = "-"),
    basal_area = ((pi * ((
      dbh / 2
    ) ^ 2)) / 10000),
    expansion_factor = case_when(
      tree_type %in% c("T", "ET") ~ 10000 / tree_plot_area,
      tree_type %in% c("S1", "ES1") ~ 10000 / sapling_plot_area,
      tree_type %in% c("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10") ~ 10000 / (regen_plot_area * number_regen_plots),
      .default = NA
    ),
    basal_area2 = basal_area * expansion_factor
  ) %>%
  group_by(measurement_year, company_plot_number) %>%
  arrange(desc(basal_area2), .by_group = TRUE) %>%
  mutate(balt = lag(cumsum(basal_area2), default = 0)) %>%
  # group_by(measurement_year, company_plot_number, balt_groups) %>% 
  # mutate(balt2 = lag(cumsum(basal_area2), default = 0)) %>%
  ungroup() %>%
  # select only those columns that are needed (just for convenience)
  dplyr::select(
    company_plot_number,
    establishment_year,
    measurement_year,
    measurement_number,
    measurement_interval,
    attack_yr,
    tree_number,
    tree_id,
    species,
    species_groups,
    dbh,
    basal_area,
    basal_area2,
    balt,
    survival
  ) -> dbh_sub


dbh_sub %>% 
  group_by(company_plot_number,establishment_year) %>% 
  summarise(surv = mean(survival)) %>% pull(surv) %>% boxplot()


# this is a small function I wrote
# this function will select measurement years (x) that whose differences (y) (attack_yr - measurement_year) are larger than zero
# it then can pick the required number of prior measurements (n) for the final selection
psp_selection <-  function(x, y, n) {
  s <- x[y > 0]
  u <- unique(s)
  n <- as.integer(n)
  sort(u, decreasing = TRUE)[n]
}

# for this this graph I included all measurements to create an overview
# if I specific subset is required according to the function above (psp_selection)
# uncomment the code below (between "dbh_sub %>%" and   "mutate(" )
dbh_sub %>%
  # group by plot to make sure everything happens within each plot separately
  #group_by(company_plot_number) %>%
  # filter out company plots that contained negative diffs, i.e. only post attack measurements
  # if this is not being done the filter step with the psp_selection function will throw an error
  #filter(!all(diff < 0)) %>%
  #filter(measurement_year >= psp_selection(measurement_year, diff, 4)) %>%
  # add another set of grouping variables on the new selection to differentiate between pre attack, attack, and post attack data points
  mutate(
    attack_label = case_when(
      measurement_year >= attack_yr ~ "Post attack",
      measurement_year < attack_yr ~ "Pre attack",
      .default = NA
    ),
    attack_label = factor(attack_label, levels = c("Pre attack", "Post attack"))
  ) %>%
  group_by(company_plot_number, attack_label) %>%
  # this one is needed for the adjusted logistic function (see logistic regression analyses below)
  mutate(measurement_interval2 = measurement_year - min(measurement_year) + 1) %>%
  # filter step to exclude cases where establishment_year is measurement_year
  # as well as to exclude measurements recorded during the attack year
  filter(measurement_year != establishment_year) %>%
  ungroup() %>% 
  # create diff helper variable
  mutate(diff_pre = attack_yr - measurement_year,
         diff_post = measurement_year - attack_yr) -> dbh0

cmi <- read_csv("MPB_climate_1991-2020_normals.csv") %>% 
  mutate(company_plot_number = as.character(company_plot_number))


cmi %>% 
  filter(company_plot_number %in% unique(dbh0$company_plot_number)) -> cmi_sub

dbh0 %>% 
  left_join(cmi_sub)-> dbh

dbh %>% filter(attack_label == "Pre attack") -> dbh_pre


# get a sense of the time intervals between measurements
boxplot(dbh_pre$measurement_interval2)
min(dbh_pre$measurement_interval2)


dbh %>% 
  filter(attack_label == "Post attack") -> dbh_post

# get a sense of the time intervals between measurements
boxplot(dbh_post$measurement_interval2)
min(dbh_post$measurement_interval2)


# data containing min / max intervals by PSP
dbh %>% 
  group_by(attack_label, company_plot_number, attack_yr) %>% 
  summarise(y_min = min(measurement_year),
            y_max = max(measurement_year)) %>% 
  ungroup() -> dbh_plot1

# data containing distinct measurement/establishment years by company plot number for geom_point below
dbh %>% 
  group_by(attack_label, company_plot_number, measurement_year, establishment_year, attack_yr) %>% 
  distinct(measurement_year) %>% 
  ungroup()-> dbh_plot2

dbh_plot1 %>% 
  ggplot(aes(y = reorder(company_plot_number, -attack_yr))) +
  geom_errorbarh(aes(xmin = y_min, xmax = y_max, color = factor(attack_yr)), linewidth = 1, alpha = 0.5) +
  geom_vline(aes(xintercept = attack_yr, color = factor(attack_yr)), linewidth = 1) +
  geom_point(data = dbh_plot2, aes(y = company_plot_number, x = measurement_year, fill = factor(attack_yr)),shape = 21,show.legend=FALSE) +
  # when un-commenting the line below the establishment year will also be shown in the plot
  #geom_point(data = dbh_plot2, aes(y = company_plot_number, x = establishment_year),shape = 21, color = "black",fill = "red", size = 2) +
  scale_x_continuous(breaks = seq(min(dbh_plot2$establishment_year, na.rm = T), max(dbh_plot2$measurement_year), by = 1)) +
  facet_wrap(~attack_label, scales = "free_x") +
  theme(axis.text.x=element_text(angle=90, vjust=.5, hjust=1)) + 
  labs(color = "MPB Attack", x = "Measurement Year", y = "Company Plot Number") -> p1

# print results (takes a few minutes because of size)
png(
  "MPB_PSP_selection_overview.png",
  res = 300,
  width = 18,
  height = 20,
  units = "in"
)
print(p1)
dev.off()

# data containing min / max intervals by PSP
dbh_pre %>% 
  group_by(attack_label, company_plot_number, attack_yr) %>% 
  summarise(y_min = min(measurement_year),
            y_max = max(measurement_year)) %>% 
  ungroup() -> dbh_plot1_pre

# data containing distinct measurement/establishment years by company plot number
# for geom_plot
dbh_pre %>% 
  group_by(attack_label, company_plot_number, measurement_year, establishment_year, attack_yr) %>% 
  distinct(measurement_year) %>% 
  ungroup()-> dbh_plot2_pre

dbh_plot1_pre %>% 
  ggplot(aes(y = reorder(company_plot_number, -attack_yr))) +
  geom_errorbarh(aes(xmin = y_min, xmax = y_max, color = factor(attack_yr)), linewidth = 1, alpha = 0.5) +
  geom_vline(aes(xintercept = attack_yr, color = factor(attack_yr)), linewidth = 1) +
  geom_point(data = dbh_plot2_pre, aes(y = company_plot_number, x = measurement_year, fill = factor(attack_yr)),shape = 21,show.legend=FALSE) +
  # when un-commenting the line below the establishment year will also be shown in the plot
  #geom_point(data = dbh_plot2, aes(y = company_plot_number, x = establishment_year),shape = 21, color = "black",fill = "red", size = 2) +
  scale_x_continuous(breaks = seq(min(dbh_plot2$establishment_year, na.rm = T), max(dbh_plot2$measurement_year), by = 1)) +
  theme(axis.text.x=element_text(angle=90, vjust=.5, hjust=1)) + 
  labs(color = "MPB Attack", x = "Measurement Year", y = "Company Plot Number") -> p2

# print results (takes a few minutes because of size)
png(
  "MPB_PSP_selection_pre.png",
  res = 300,
  width = 18,
  height = 20,
  units = "in"
)
print(p2)
dev.off()

# data containing min / max intervals by PSP
dbh_post %>% 
  group_by(attack_label, company_plot_number, attack_yr) %>% 
  summarise(y_min = min(measurement_year),
          y_max = max(measurement_year)) %>% 
  ungroup() -> dbh_plot1_post

# data containing distinct measurement/establishment years by company plot number
# for geom_plot
dbh_post %>% 
  group_by(attack_label, company_plot_number, measurement_year, establishment_year, attack_yr) %>% 
  distinct(measurement_year) %>% 
  ungroup()-> dbh_plot2_post

dbh_plot1_post %>% 
  ggplot(aes(y = reorder(company_plot_number, -attack_yr))) +
  geom_errorbarh(aes(xmin = y_min, xmax = y_max, color = factor(attack_yr)), linewidth = 1, alpha = 0.5) +
  geom_vline(aes(xintercept = attack_yr, color = factor(attack_yr)), linewidth = 1) +
  geom_point(data = dbh_plot2_post, aes(y = company_plot_number, x = measurement_year, fill = factor(attack_yr)),shape = 21,show.legend=FALSE) +
  # when un-commenting the line below the establishment year will also be shown in the plot
  #geom_point(data = dbh_plot2, aes(y = company_plot_number, x = establishment_year),shape = 21, color = "black",fill = "red", size = 2) +
  scale_x_continuous(breaks = seq(min(dbh_plot2$establishment_year, na.rm = T), max(dbh_plot2$measurement_year), by = 1)) +
  theme(axis.text.x=element_text(angle=90, vjust=.5, hjust=1)) + 
  labs(color = "MPB Attack", x = "Measurement Year", y = "Company Plot Number") -> p3


# print results (takes a few minutes because of size)
png(
  "MPB_PSP_selection_post.png",
  res = 300,
  width = 18,
  height = 20,
  units = "in"
)
print(p3)
dev.off()



######### LOGISTIC REGRESSION #########

# Ben Bolker's hack to model exposure adjusted logistic function (straight copy and paste)
# https://rpubs.com/bbolker/logregexp
# Also read here for the number successes and failures required for the model to run well
# https://github.com/lme4/lme4/issues/179#issuecomment-428411405

logexp <- function(exposure = 1) {
  ## hack to help with visualization, post-prediction etc etc
  get_exposure <- function() {
    if (exists("..exposure", envir = .GlobalEnv))
      return(get("..exposure", envir = .GlobalEnv))
    exposure
  }
  linkfun <- function(mu)
    qlogis(mu ^ (1 / get_exposure()))
  ## FIXME: is there some trick we can play here to allow
  ##   evaluation in the context of the 'data' argument?
  linkinv <- function(eta)
    plogis(eta) ^ get_exposure()
  logit_mu_eta <- function(eta) {
    ifelse(abs(eta) > 30, .Machine$double.eps,
           exp(eta) / (1 + exp(eta)) ^ 2)
  }
  mu.eta <- function(eta) {
    get_exposure() * plogis(eta) ^ (get_exposure() - 1) *
      logit_mu_eta(eta)
  }
  valideta <- function(eta)
    TRUE
  link <- paste("logexp(", deparse(substitute(exposure)), ")",
                sep = "")
  structure(
    list(
      linkfun = linkfun,
      linkinv = linkinv,
      mu.eta = mu.eta,
      valideta = valideta,
      name = link
    ),
    class = "link-glm"
  )
}

#### Ben Bolker's hack end ####

#######################################

# subset 

# GLM models no random effects
# select species and 
species <- c("Aw", "Sw", "Sb", "Pl")
attack_label <-
  unique(dbh$attack_label) %>% 
  droplevels()

# create empty data frames to fill in regression output
d_glm <- data.frame()
d_glm2 <- data.frame()

# this for loop will run models one species at a time and within each attack_label separately
for (k in attack_label) {
  dbh %>% filter(attack_label == k) -> dbh2
  for (i in species) {
    dbh2 %>% filter(species == i) -> d
    m <-
      glm(
        survival ~ scale(dbh) + scale(dbh ^ 2) + scale(balt)  + scale(CMI),
        family = binomial(logexp(d$measurement_interval2)),
        data = d
      )
    ..exposure <- mean(d$measurement_interval2)
    dd <-
      data.frame(ggeffects::ggpredict(m, terms = c("dbh [all]")))
    dd2 <-
      data.frame(ggeffects::ggpredict(m, terms = c("balt [all]")))
    dd$species <- d$species[1]
    dd$attack_label <- d$attack_label[1]
    dd2$species <- d$species[1]
    dd2$attack_label <- d$attack_label[1]
    rm(..exposure)
    d_glm <- rbind(d_glm, dd)
    d_glm2 <- rbind(d_glm2, dd2)
  }
}

# plots results
ggplot(d_glm, aes(x, predicted, group = species)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              color = "grey90",
              alpha = 0.05) +
  geom_line(aes(color = species), linewidth = 1) +
  scale_y_continuous(labels = scales::percent, breaks = seq(0,1,0.25), limits = c(0,1)) +
  labs(x = "DBH (cm)", y = "Survival", color = "Species") +
  theme_bw() +
  facet_wrap(~attack_label) +
  ggtitle("GLM") -> p1_glm

ggplot(d_glm2, aes(x, predicted, group = species)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              color = "grey90",
              alpha = 0.05) +
  geom_line(aes(color = species), linewidth = 1) +
  scale_y_continuous(labels = scales::percent, breaks = seq(0,1,0.25), limits = c(0,1)) +
  labs(x = "BALT", y = "Survival", color = "Species") +
  theme_bw() +
  facet_wrap(~attack_label) +
  ggtitle("GLM") -> p2_glm


p1_glm | p2_glm



#### Now for the GLMM #######

species <- c("Aw", "Sw", "Sb", "Pl")
attack_label <- unique(dbh$attack_label) %>% 
  droplevels()

d_glmer <- data.frame()
d_glmer2 <- data.frame()

for (k in attack_label) {
  dbh %>% filter(attack_label == k) -> dbh2
  for (i in species) {
    dbh2 %>% filter(species == i) -> d
    m <-
      glmer(
        survival ~ scale(dbh) +  scale(dbh ^ 2) + scale(balt) + scale(CMI) + (1 | company_plot_number) ,
        family = binomial(logexp(d$measurement_interval2)),
        glmerControl(optimizer = "Nelder_Mead", optCtrl = list(maxfun = 2e4)),
        data = d
      )
    ..exposure <- mean(d$measurement_interval2)
    dd <- data.frame(ggeffects::ggpredict(m, terms = c("dbh [all]")))
    dd2 <- data.frame(ggeffects::ggpredict(m, terms = c("balt [all]")))
    dd$species <- d$species[1]
    dd$attack_label <- d$attack_label[1]
    dd2$species <- d$species[1]
    dd2$attack_label <- d$attack_label[1]
    rm(..exposure)
    d_glmer <- rbind(d_glmer, dd)
    d_glmer2 <- rbind(d_glmer2, dd2)
  }
}

# plot results
ggplot(d_glmer, aes(x, predicted, group = species)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              color = "grey90",
              alpha = 0.05) +
  geom_line(aes(color = species), linewidth = 1) +
  scale_y_continuous(labels = scales::percent, breaks = seq(0,1,0.25), limits = c(0,1)) +
  labs(x = "DBH (cm)", y = "Survival", color = "Species") +
  facet_wrap(~attack_label) +
  theme_bw() +
  ggtitle("GLMM\n(random effect: plot)") -> p1_glmer


ggplot(d_glmer2, aes(x, predicted, group = species)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              color = "grey90",
              alpha = 0.05) +
  geom_line(aes(color = species), linewidth = 1) +
  scale_y_continuous(labels = scales::percent, breaks = seq(0,1,0.25), limits = c(0,1)) +
  labs(x = "BALT", y = "Survival", color = "Species") +
  facet_wrap(~attack_label) +
  theme_bw() +
  ggtitle("GLMM\n(random effect: plot)") -> p2_glmer

p1_glmer | p2_glmer

# patch plots together
p <- (p1_glm | p1_glmer) / (p2_glm | p2_glmer)

# print results
png(
  "survival_post_attack_cmi.png",
  res = 900,
  width = 12,
  height = 10,
  units = "in"
)
print(p)
dev.off()

