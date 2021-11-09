# new case 1 figures

# setup =============

library(lme4)
library(ggtext)
library(emmeans)
library(tidymodels)
source("R/diversity_data_prep.R")
# data prep ============

# data mongering 

uninvaded_sites <- all_scales %>% 
  mutate(year = as.numeric(year),
         uniqueid = paste0(year+1,plotID,scale,subplotID, site)) %>%
  filter(nspp_exotic == 0) %>% 
  dplyr::select(year, plotID, scale, subplotID, site,uniqueid, 
                shannon_total, nspp_total, shannon_native, nspp_native)

uniqueids <- uninvaded_sites$uniqueid

next_year<-all_scales %>% 
  mutate(year = as.numeric(year),
         uniqueid = paste0(year,plotID,scale,subplotID, site))%>%
  filter(uniqueid %in% uniqueids) %>%
  dplyr::select(uniqueid, next_shannon_total=shannon_total, 
                next_nspp_total=nspp_total, 
                next_nspp_exotic = nspp_exotic,
                next_shannon_native=shannon_native, 
                next_nspp_native = nspp_native,
                next_shannon_exotic = shannon_exotic) %>%
  mutate(invaded = ifelse(next_nspp_exotic > 0, 1, 0))

prev_year_div <- left_join(next_year, uninvaded_sites)

# models ====================
p1_mods <- all_scales %>% 
  mutate(invaded = ifelse(invaded=="invaded", 1,0)) %>%
  group_nest(scale)  %>%
  mutate(model = map(data, ~glm(invaded ~ nspp_native, 
                               data=.x, family="binomial")))%>%
  mutate(fit = map(model, ~predict(.x,se.fit=TRUE, type = "response")[1]),
         se = map(model, ~predict(.x,se.fit=TRUE, type = "response")[2]))

p1dat <- unnest(p1_mods, cols=c(fit, se)) 

p1_dat<-p1dat %>% unnest(data) %>%
  mutate(fit = (unlist(p1dat$fit)),
         se = unlist(p1dat$se)) 


# plots ===================
# didn't bother with the other ones, since it looks like everything will be 
# identical to the originals.
p1<-ggplot(p1_dat, 
           aes(x = nspp_native, y=invaded, color = scale)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  geom_line(aes(y=fit, group = scale), color = "grey", lty=2)+
  theme_classic() +
  theme(legend.position = c(1,0),
        legend.justification = c(1,0),
        panel.border = element_rect(fill = NA, size=0.75),
        legend.background = element_rect(fill="transparent"))+
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.5,1)) +
  scale_color_viridis_d(option = "B") +
  xlab("Native Species Richness") +
  ylab("P(Invaded)")

# tables ------------

p1_tabs <- all_scales %>% 
  mutate(invaded = ifelse(invaded=="invaded", 1,0)) %>%
  group_nest(scale)  %>%
  mutate(model = map(data, ~glm(invaded ~ nspp_native, 
                                data=.x, family="binomial")),
         tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(formula = map(model, formula) %>% as.character,
         formula = str_c("Figure 2b: ",formula)) %>%
  dplyr::select(-data,-model)

p2_tabs<-all_scales  %>% 
  mutate(invaded = ifelse(invaded=="invaded", 1,0))%>%
  group_nest(scale)  %>%
  mutate(model = map(data, ~glm(nspp_exotic ~ nspp_native, 
                                data=.x, family = "quasipoisson")),
         tidied = map(model, tidy)) %>%
  unnest(tidied)%>%
  mutate(formula = map(model, formula) %>% as.character,
         formula = str_c("Figure 2a: ", formula)) %>%
  dplyr::select(-data,-model)


p3_tabs<-prev_year_div %>%
  group_nest(scale)  %>%
  mutate(model = map(data, ~glm(invaded ~ nspp_native, 
                                data=.x, family = "binomial")),
         tidied = map(model, tidy)) %>%
  unnest(tidied)%>%
  mutate(formula = map(model, formula) %>% as.character,
         formula = str_c("Figure 2d: ", formula)) %>%
  dplyr::select(-data,-model)


p4_tabs<-prev_year_div %>%
  group_nest(scale)  %>%
  mutate(model = map(data, ~glm(next_nspp_exotic ~ nspp_native, 
                                data=.x, family = "quasipoisson")),
         tidied = map(model, tidy)) %>%
  unnest(tidied)%>%
  mutate(formula = map(model, formula) %>% as.character,
         formula = str_c("Figure 2c: ", formula)) %>%
  dplyr::select(-data,-model)


bigtab <- bind_rows(p1_tabs, p2_tabs, p3_tabs,p4_tabs) %>%
  dplyr::select(formula, scale, term, estimate, std.error, statistic, p.value) %>%
  mutate(term = str_replace_all(term,"nspp_native","Native Species Richness"),
         formula = str_replace_all(formula, "nspp_native", "Native Species Richness"),
         formula = str_replace_all(formula, "invaded", "P(Invaded)"),
         formula = str_replace_all(formula, "nspp_exotic", "Non-Native Species Richness"),
         formula = str_replace_all(formula, "next_",""),
         estimate = signif(estimate,2),
         std.error = signif(std.error,(2)),
         statistic = signif(statistic,(2)),
         p.value = signif(p.value,(2)),
         p.value = ifelse(p.value < 0.001, "< 0.0001", p.value)) %>%
  arrange(formula, scale)

write_csv(bigtab, "tables/split_models.csv")
