


#### libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(mgcv)
library(popbio)



#### read data file
dat_raw <- read_csv('data/lemna-temp-light-2013.csv') %>% 
  mutate(treat = paste(light, temp, sep = '-'))



#### subset to fronds not contaminated prior to death
dat <- filter(dat_raw, !is.na(date_death))



#### light test
dat_light <- read_csv('data/light-test-2013-09-24.csv') %>%
  setNames(c('N', 'Time', 'PAR', 'Event'))

dat_light_summ <- dat_light %>% 
  group_by(Event) %>% 
  summarize(par_mean = mean(PAR))

mean(dat_light_summ$par_mean[c(1, 3, 5)]) # mean PAR bright
mean(dat_light_summ$par_mean[c(2, 4, 6)]) # mean PAR dim



#### extract reproduction data from main data frame
# number of offspring produced each day by each frond
dat_raw_repro <- subset(dat, select = May_09_2013:Jun_23_2013) %>% 
  apply(2, as.integer)




#### determine dates of first and last reproduction
# get column name (i.e. date) associated with first (min) or last (max) repro
GetDateFirstRepro <- function (x) colnames(dat_raw_repro)[min(which(x > 0))]
GetDateLastRepro <- function (x) colnames(dat_raw_repro)[max(which(x > 0))]
dat$date_first_repro <- apply(dat_raw_repro, 1, GetDateFirstRepro)
dat$date_last_repro <- apply(dat_raw_repro, 1, GetDateLastRepro)

# convert dates to R's date format
dat$date_birth <- as.Date(dat$date_birth, format = "%b_%d_%Y")
dat$date_first_repro <- as.Date(dat$date_first_repro, format = "%b_%d_%Y")
dat$date_last_repro <- as.Date(dat$date_last_repro, format = "%b_%d_%Y")




#### calculate frond-level demographic statistics
# 'ffr' is from first reproduction, 'latency' is latency to reproduce
dat$lifespan <- as.numeric(dat$date_last_repro - dat$date_birth + 1)
dat$lifespan_ffr <- as.numeric(dat$date_last_repro - dat$date_first_repro + 1)
dat$latency <- as.numeric(dat$date_first_repro - dat$date_birth + 1)
dat$total_offspring <- rowSums(dat_raw_repro, na.rm = T)
dat$fecund_mean <- dat$total_offspring / dat$lifespan
dat$fecund_mean_ffr <- dat$total_offspring / dat$lifespan_ffr




#### sample size
dat_raw %>%    # initial
  group_by(light, temp) %>% 
  summarize(n = n()) %>% ungroup()

dat %>%        # final
  group_by(light, temp) %>% 
  summarize(n = n()) %>% ungroup()




#### histogram of lifespan in hot vs. cold treatments
dat_plot <- dat %>%
  mutate(temp_lab = ifelse(temp == 'hot', 'Hot (28°C)', 'Cold (22°C)'))


tt <- theme(axis.title = element_text(size = 18),
            axis.text = element_text(size = 15),
            strip.text = element_text(size = 15),
            panel.grid = element_blank())

p1 <- ggplot(dat_plot) +
  geom_histogram(aes(lifespan), bins = 10, fill = '#0570b0') +
  facet_wrap(~ temp_lab, ncol = 1) +
  xlab('Lifespan (days)') +
  ylab('Count') +
  tt

dev.off()
quartz(height = 8, width = 8)
print(p1)

# ggsave('img/fig1.png', height = 8, width = 8, units = 'in', dpi = 300)


# summary of lifespan (mean, sd, se) by temperature treatment
dat %>%
  group_by(temp) %>% 
  summarize(life_mean = mean(lifespan),
            life_sd = sd(lifespan),
            life_se = sd(lifespan) / sqrt(n()))




#### bootstrap shape_mortality
BootShape <- function(data, n_rep) {
  t <- data$lifespan
  tResampled <- replicate(n_rep, sample(t, length(t), replace = T), simplify = F)
  Shape <- function(y) 1 - (sd(y) / mean(y))
  bootShape <- sapply(tResampled, Shape)
  return(data.frame(bootShape))
}

shapeTreat <- dat %>%
  group_by(light, temp) %>%
  do(BootShape(., n_rep = 10000)) %>%
  ungroup() %>% 
  mutate(treat = paste(light, temp, sep = '-'))



#### flatten dat
dat_flat <- dat %>% 
  dplyr::select(id, light, temp, date_birth, date_last_repro, fecund_mean, fecund_mean_ffr, matches('2013')) %>% 
  gather(date, fecund, matches('2013')) %>%
  arrange(id) %>% 
  mutate(date = as.Date(date, format = '%b_%d_%Y')) %>% 
  filter(date >= date_birth & date <= date_last_repro) %>%
  mutate(fecund = as.numeric(fecund),
         age = as.numeric(date - date_birth + 1),
         died = ifelse(date == date_last_repro, 1, 0),
         surv = ifelse(date == date_last_repro, 0, 1)) %>% 
  dplyr::select(-contains('date'))




#### intrinsic rates of increase
IntrinsicIncrease <- function (age, fecund) {       # function to calculate r from flat data
  n <- max(age)                                     # get lifespan
  lesMat <- matrix(0, nrow = n, ncol = n)           # create empty lifespan x lifespan Leslie matrix
  lesMat[1,] <- fecund                              # fill first row with reproduction data
  lesMat[2:n, 1:(n-1)] <- diag(rep(1, (n-1)))       # fill subdiagonal with survival data
  lambda <- lambda(lesMat)
  return(data.frame(r = log(lambda)))               # return ln of leading eigenvalue of Leslie matrix
}

dat_r <- dat_flat %>% 
  group_by(id) %>% 
  do(IntrinsicIncrease(age = .$age, fecund = .$fecund)) %>% 
  ungroup()

dat$r <- dat_r$r



#### group life expectancy
life_exp <- dat %>%
  group_by(light, temp) %>% 
  summarize(life_exp = mean(lifespan, na.rm = T),
            life_exp_ffr = mean(lifespan_ffr, na.rm = T)) %>% 
  ungroup()



#### mean-standardize age and fecundity
ScaleFn <- function(data, life_exp) {
  data %>%
    left_join(life_exp, by = c('light', 'temp')) %>% 
    mutate(age_std = age / life_exp,
           fecund_std = fecund / fecund_mean)
}

ScaleFnFfr <- function(data, life_exp) {
  data %>%
    group_by(id, light, temp) %>% 
    mutate(first_repro = which(fecund > 0)[1]) %>% 
    filter(age >= first_repro) %>% 
    ungroup() %>%
    left_join(life_exp, by = c('light', 'temp')) %>% 
    mutate(age_std = age / life_exp_ffr,
           fecund_std = fecund / fecund_mean_ffr)
}

dat_scale <- ScaleFn(dat_flat, life_exp) %>%
  mutate(treat = paste(light, temp, sep = '-'))

dat_scale_ffr <- ScaleFnFfr(dat_flat, life_exp) %>%
  mutate(treat = paste(light, temp, sep = '-'))



#### shape fecundity
ShapeFecundityLme <- function(data) {
  mod_std <- lme(fecund_std ~ age_std, random = ~ age_std|id, data = data, method = 'ML')
  return(as.data.frame(t(intervals(mod_std, which = 'fixed')$fixed[2,])))
}

shape_fecund <- dat_scale_ffr %>% 
  group_by(light, temp) %>% 
  do(ShapeFecundityLme(data = .)) %>% 
  ungroup() %>% 
  rename(med = est., low = lower, upp = upper)

dat_scale_ffr %>% 
  group_by(temp) %>% 
  do(ShapeFecundityLme(data = .)) %>% 
  ungroup() %>% 
  rename(med = est., low = lower, upp = upper)


# test
mod1 <- lme(fecund_std ~ age_std * treat, random = ~ age_std|id, data = dat_scale_ffr, method = 'ML')
mod2 <- lme(fecund_std ~ age_std * temp, random = ~ age_std|id, data = dat_scale_ffr, method = 'ML')
mod3 <- lme(fecund_std ~ age_std * light, random = ~ age_std|id, data = dat_scale_ffr, method = 'ML')
mod4 <- lme(fecund_std ~ age_std * temp * light, random = ~ age_std|id, data = dat_scale_ffr, method = 'ML')

summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)

anova(mod1)
anova(mod2)
anova(mod3)
anova(mod4)

# note we effectively have finer temporal resolution for high temp... problem?



#### survivorship and probability of survival by age and subsite
Survship <- function(age_vec, df, std = F) {
  n = nrow(df)
  t <- df$lifespan
  if (std) t <- df$lifespan / mean(df$lifespan)
  surv = sapply(age_vec, function (x) length(which(t >= x)) / n)
  data.frame(age = age_vec, psurv = surv)
}

survTreat <- group_by(dat, light, temp) %>%
  do(Survship(age_vec = 1:40, df = .)) %>% 
  mutate(treat = paste(light, temp, sep = '-')) %>% 
  mutate(treat = factor(treat, levels = c('bright-hot', 'dim-hot', 'bright-cold', 'dim-cold')))

survTreatStd <- group_by(dat, light, temp) %>%
  do(Survship(age_vec = seq(0, 1.8, 0.02), df = ., std = T)) %>% 
  mutate(treat = paste(light, temp, sep = '-')) %>% 
  mutate(treat = factor(treat, levels = c('bright-hot', 'dim-hot', 'bright-cold', 'dim-cold')))



#### hazard curves by treatment
GetHazardCurve <- function(data, k) {
  life_mean <- mean(filter(data, died == 1)$age)
  
  mod <- gam(surv ~ s(age, k = k), family = 'binomial', data = data)
  mod_std <- gam(surv ~ s(age_std, k = k), family = 'binomial', data = data)
  
  newage <- seq(0, max(data$age), length = 100)
  newagestd <- seq(0, max(data$age_std), length = 100)
  
  pred_y <- predict(mod, newdata = data.frame(age = newage), type = 'response')
  pred_y_se <- predict(mod, newdata = data.frame(age = newage), type = 'response', se = T)$se.fit
  
  pred_y_std <- predict(mod_std, newdata = data.frame(age_std = newagestd), type = 'response')
  pred_y_std_se <- predict(mod_std, newdata = data.frame(age_std = newagestd), type = 'response', se = T)$se.fit
  
  pred_y_std_life <- as.numeric(predict(mod_std, newdata = data.frame(age_std = 1), type = 'response'))
  
  df <- data.frame(
    newage,
    newagestd,
    pred_y_std_life,
    haz = -log(pred_y),
    haz_upper = -log(pred_y + pred_y_se),
    haz_lower = -log(pred_y - pred_y_se),
    haz_std = (-log(pred_y_std)) * life_mean,
    haz_std_upper = (-log(pred_y_std + pred_y_std_se)) * life_mean,
    haz_std_lower = (-log(pred_y_std - pred_y_std_se)) * life_mean
  )
  
  return(df)
}

hazard_curves <- dat_scale %>%
  group_by(treat) %>%
  do(GetHazardCurve(., k = 3)) %>% 
  ungroup()



#### plot survivorship, mortality and fecundity trajectories by treatment
# ggplot theme
tt <- theme_bw() +
  theme(axis.title = element_text(size = 16),
        plot.title = element_text(size = 16, hjust = 0.5),
        text = element_text(size = 15),
        panel.grid = element_blank(),
        axis.title.x = element_text(margin = margin(.4, 0, 0, 0, unit = 'cm')),
        axis.title.y = element_text(margin = margin(0, .4, 0, 0, unit = 'cm')))

# function for labelling axes
LogLab <- function(x) ifelse(x %in% c(0.1, 1), as.character(x), '')

s_bri <- 2.3
s_dim <- 1

br <- c(0.008, 0.009, (1:10)/100, (1:10)/10)
treat_lab <- c('Hot (28°C), Bright', 'Hot (28°C), Dim', 'Cold (22°C), Bright', 'Cold (22°C), Dim')

# plot
p1_1 <- ggplot(survTreat, aes(age, psurv, size = treat, col = treat)) +
  geom_step(alpha = 0.85, direction = 'hv') +
  coord_cartesian(xlim = c(0, 40), ylim = c(0.02, 1)) +
  scale_size_manual(values = c(s_bri, s_dim, s_bri, s_dim), name = 'Treatment', guide = 'legend', labels = treat_lab) +
  scale_color_manual(values = c('darkred', 'darkred', 'darkblue', 'darkblue'), name = 'Treatment', guide = 'legend', labels = treat_lab) +
  scale_y_log10(breaks = br, labels = LogLab) +
  xlab(NULL) +
  ylab('Survivorship') +
  ggtitle('Non-standardized') +
  annotate('text', x = Inf, y = Inf, label = '(a)', hjust = 1.4, vjust = 1.7, size = 4.6) +
  tt + theme(legend.position = c(.24, .27))

p1_2 <- ggplot(survTreatStd, aes(age, psurv, size = treat, col = treat)) +
  geom_step(alpha = 0.85, direction = 'hv') +
  coord_cartesian(xlim = c(0, 1.5), ylim = c(0.02, 1)) +
  scale_size_manual(values = c(s_bri, s_dim, s_bri, s_dim), name = 'Treatment', guide = F) +
  scale_color_manual(values = c('darkred', 'darkred', 'darkblue', 'darkblue'), name = 'Treatment', guide = F) +
  scale_x_continuous(labels = function(x) formatC(x, width = 1)) +
  scale_y_log10(breaks = br, labels = LogLab) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle('Standardized') +
  annotate('text', x = Inf, y = Inf, label = '(b)', hjust = 1.4, vjust = 1.6, size = 4.5) +
  tt

p1_3 <- ggplot(hazard_curves, aes(x = newage, y = haz, size = treat, col = treat)) +
  geom_ribbon(inherit.aes = F, aes(x = newage, ymin = haz_lower, ymax = haz_upper, group = treat), fill = 'grey40', alpha = 0.2)+
  geom_line() +
  coord_cartesian(xlim = c(0, 40), ylim = c(0, 1)) +
  scale_size_manual(values = c(s_bri, s_bri, s_dim, s_dim), name = 'Treatment', guide = F) +
  scale_color_manual(values = c('darkblue', 'darkred', 'darkblue', 'darkred'), name = 'Treatment', guide = F) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = function(x) formatC(x, width = 1)) +
  xlab(NULL) +
  ylab('Mortality') +
  annotate('text', x = Inf, y = Inf, label = '(c)', hjust = 1.4, vjust = 1.7, size = 4.6) +
  tt

p1_4 <- ggplot(hazard_curves, aes(x = newagestd, y = haz_std, size = treat, col = treat)) +
  geom_ribbon(inherit.aes = F, aes(x = newagestd, ymin = haz_std_lower, ymax = haz_std_upper, group = treat), fill = 'grey40', alpha = 0.2)+
  geom_line() +
  scale_size_manual(values = c(s_bri, s_bri, s_dim, s_dim), name = 'Treatment', guide = F) +
  scale_color_manual(values = c('darkblue', 'darkred', 'darkblue', 'darkred'), name = 'Treatment', guide = F) +
  scale_x_continuous(labels = function(x) formatC(x, width = 1)) +
  scale_y_continuous(breaks = c(0, 10, 20)) +
  coord_cartesian(xlim = c(0, 1.5), ylim = c(0, 25)) +
  xlab(NULL) +
  ylab(NULL) +
  annotate('text', x = Inf, y = Inf, label = '(d)', hjust = 1.4, vjust = 1.7, size = 4.6) +
  tt

p1_5 <- ggplot(dat_scale, aes(age, fecund, size = treat, col = treat)) +
  geom_smooth(inherit.aes = F, aes(age, fecund, group = treat), method = 'loess', se = T, size = 0, fill = 'grey40', alpha = 0.2) +
  geom_smooth(method = 'loess', se = F) +
  coord_cartesian(xlim = c(0, 40), ylim = c(0, 0.8)) +
  scale_size_manual(values = c(s_bri, s_bri, s_dim, s_dim), name = 'Treatment', guide = F) +
  scale_color_manual(values = c('darkblue', 'darkred', 'darkblue', 'darkred'), name = 'Treatment', guide = F) +
  scale_y_continuous(breaks = seq(0, 1, 0.5), labels = function(x) formatC(x, width = 1)) +
  xlab('Age (days)') +
  ylab('Fecundity') +
  annotate('text', x = Inf, y = Inf, label = '(e)', hjust = 1.4, vjust = 1.7, size = 4.6) +
  tt

p1_6 <- ggplot(dat_scale, aes(age_std, fecund_std, size = treat, col = treat)) +
  geom_smooth(inherit.aes = F, aes(age_std, fecund_std, group = treat), method = 'loess', se = T, size = 0, fill = 'grey40', alpha = 0.2) +
  geom_smooth(method = 'loess', se = F) +
  coord_cartesian(xlim = c(0, 1.5), ylim = c(0, 1.5)) +
  scale_size_manual(values = c(s_bri, s_bri, s_dim, s_dim), name = 'Treatment', guide = F) +
  scale_color_manual(values = c('darkblue', 'darkred', 'darkblue', 'darkred'), name = 'Treatment', guide = F) +
  scale_x_continuous(labels = function(x) formatC(x, width = 1)) +
  scale_y_continuous(labels = function(x) formatC(x, width = 1)) +
  xlab('Age (life expectancies)') +
  ylab(NULL) +
  annotate('text', x = Inf, y = Inf, label = '(f)', hjust = 1.4, vjust = 1.7, size = 4.6) +
  tt

g1_1 <- ggplotGrob(p1_1)
g1_2 <- ggplotGrob(p1_2)
g1_3 <- ggplotGrob(p1_3)
g1_4 <- ggplotGrob(p1_4)
g1_5 <- ggplotGrob(p1_5)
g1_6 <- ggplotGrob(p1_6)

g1_3$widths <- g1_1$widths
g1_5$widths <- g1_1$widths
g1_2$widths <- g1_6$widths
g1_4$widths <- g1_6$widths

f2 <- arrangeGrob(g1_1, g1_2, g1_3, g1_4, g1_5, g1_6, nrow = 3,
                  heights = c(1.09, 0.98, 1.09), widths = c(1.062, 1))

dev.off()
quartz(height = 10, width = 10)
grid.arrange(f2)

# ggsave('img/fig2.png', f2, height = 10, width = 10, units = 'in', dpi = 300)






#### life history traits by treatment
dat_plot <- dat %>%
  dplyr::select(id, treat, lifespan, total_offspring, fecund_mean, r) %>% 
  gather(variable, value, lifespan:r)

ggplot(dat_plot) +
  geom_boxplot(aes(treat, value)) +
  facet_wrap(~ variable, scales = 'free_y') +
  theme_bw()

# stats
summary(lm(scale(lifespan) ~ temp * light, dat))
summary(lm(scale(fecund_mean) ~ temp * light, dat))
summary(lm(scale(total_offspring) ~ temp * light, dat))

anova(lm(scale(lifespan) ~ temp * light, dat))
anova(lm(scale(fecund_mean) ~ temp * light, dat))
anova(lm(scale(total_offspring) ~ temp * light, dat))




#### stats for shape_mortality
# main effect
aov_main <- shapeTreat %>%
  group_by(treat) %>% 
  summarize(mean = mean(bootShape),
            var = var(bootShape))

var_g <- sum((aov_main$mean - mean(aov_main$mean))^2) / (nrow(aov_main) - 1)
var_w <- mean(aov_main$var)
var_g/var_w
pf(var_g/var_w, 3, nrow(dat) - 3, lower.tail = F)

# temp
aov_temp <- shapeTreat %>%
  group_by(temp) %>% 
  summarize(mean_temp = mean(bootShape),
            var_temp = var(bootShape))

var_g <- sum((aov_temp$mean_temp - mean(aov_temp$mean_temp))^2) / (nrow(aov_temp) - 1)
var_w <- mean(aov_temp$var_temp)
var_g/var_w
pf(var_g/var_w, 1, nrow(dat) - 3, lower.tail = F)

# light
aov_light <- shapeTreat %>%
  group_by(light) %>% 
  summarize(mean_light = mean(bootShape),
            var_light = var(bootShape))

var_g <- sum((aov_light$mean_light - mean(aov_light$mean_light))^2) / (nrow(aov_light) - 1)
var_w <- mean(aov_light$var_light)
var_g/var_w
pf(var_g/var_w, 1, nrow(dat) - 3, lower.tail = F)


# interaction
aov_interact <- aov_main %>% 
  left_join(unique(dplyr::select(dat, treat, light, temp)), by = 'treat') %>% 
  left_join(aov_temp, by = 'temp') %>% 
  left_join(aov_light, by = 'light')

var_g <- sum((aov_interact$mean - aov_interact$mean_temp - aov_interact$mean_light + mean(aov_interact$mean))^2) / 1
var_w <- mean(aov_interact$var)
var_g/var_w
pf(var_g/var_w, 1, nrow(dat) - 3, lower.tail = F)


