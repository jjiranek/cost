### Cost of Resistance Analyses
# Load required packages
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(lme4)
library(MASS)
library(emmeans)
source("glmm_funs.R")

## Fecundity Assay
# Read in data
fecundity <- read.csv("fecundityFinal.csv", header = T)

# Code variables as factors
fecundity$Food <- factor(fecundity$Food, levels = c("DA837", "OP50","HB101"))
fecundity$Day <- factor(fecundity$Day)
fecundity$Genotype <- factor(fecundity$Genotype, levels = c("N2", "ERT250"))

# Subset sheet for daily fecundity
fecundity <- fecundity[fecundity$ID != 130,] # remove outlier (no offspring) 
fec_day<- fecundity[!is.na(fecundity$Fecundity),] # remove missing observations

# Subset sheet for lifetime fecundity
fec_lt <- fecundity[fecundity$Censored == 0,] # remove censored worms

# Sum up total offspring to calculate lifetime fecundity
ltfecund <- fec_lt %>%
  group_by(ID) %>%
  summarise(ID = first(ID),
            Fecundity = sum(Fecundity),
            Food = first(Food),
            Genotype = first(Genotype)
  )

# LT Fecundity Means by Food
ltfood <- ltfecund %>%
  group_by(Food) %>%
  summarise(Food = first(Food),
            n = n(),
            Counts = mean(Fecundity),
            SEM = sd(Fecundity)/sqrt(n)
  )

# Daily Fecundity
dayfecund <- fec_day %>%
  group_by(Day, Food) %>%
  summarise(Day = first(Day),
            Food = first(Food),
            n = n(),
            Mean = mean(Fecundity),
            SEM = sd(Fecundity)/sqrt(n)
  )

# Share of Reproduction Broken down by Day
repro_share <- fec_day %>%
  group_by(ID)%>%
  mutate(percentage = Fecundity/sum(Fecundity))%>% 
  group_by(Day, Food) %>%
  summarise(Day = first(Day),
            Food = first(Food),
            n = n(),
            Mean = mean(percentage),
            SEM = sd(percentage)/sqrt(n)
  )

# Generate summary statistics
meanfecund <- ltfecund %>%
  group_by(Food, Genotype) %>%
  summarise(n = n(),
            Mean_Fecundity = mean(Fecundity),
            SEM = sd(Fecundity)/sqrt(n), 
  )

# Model Lifetime Fecundity
null <- glm.nb(Fecundity ~ 1, data = ltfecund) # null model
mod1 <- glm.nb(Fecundity ~ Genotype*Food, data = ltfecund) # interaction
mod2 <- glm.nb(Fecundity ~ Genotype + Food , data = ltfecund) # no interaction
mod3 <- glm.nb(Fecundity ~ Food, data = ltfecund) # only food is significant

# Lifetime Fecundity Model Comparison
f.aic<-c(AIC(null),AIC(mod1),AIC(mod2),AIC(mod3))
delAIC<-f.aic-min(f.aic)
relLik <- exp(-0.5 * delAIC)
aicweight <- relLik/sum(relLik)
aic.table<-data.frame(AIC = f.aic, delAIC = delAIC, relLik = relLik,
                      weight = aicweight) 
rownames(aic.table) <- c("Null","Genotype*Food","Genotype + Food","Food")
aic.table  %>% arrange(desc(weight))

# Summary of Best Models
summary(null)
summary(mod3)

# Model Daily Fecundity
dailynull <- glmer.nb(Fecundity ~ Day + (1|ID), 
                      data = fec_day)
daily1 <- glmer(Fecundity ~ Day*Food + Genotype*Food + (1|ID), 
                data = fec_day, family = poisson ) 
overdisp_fun(daily1) # very overdispersed

daily2 <- glmer.nb(Fecundity ~ Day*Food + Genotype*Food + (1|ID), 
                   data = fec_day)
daily3 <- glmer.nb(Fecundity ~ Day*Food + Genotype + (1|ID), 
                   data = fec_day)
daily4 <- glmer.nb(Fecundity ~ Day*Food + (1|ID), 
                   data = fec_day) # best model
daily5 <- glmer.nb(Fecundity ~ Day + Food*Genotype + (1|ID), 
                   data = fec_day)

# Daily Offspring Production Model Comparison
f.aic<-c(AIC(dailynull),AIC(daily2),AIC(daily3), AIC(daily4), AIC(daily5))
delAIC<-f.aic-min(f.aic)
relLik <- exp(-0.5 * delAIC)
aicweight <- relLik/sum(relLik)
aic.table<-data.frame(AIC = f.aic, delAIC = delAIC, weight = aicweight)
rownames(aic.table) <- c("Day","Day*Food + Genotype*Food",
                         "Day*Food + Genotype", "Day*Food",
                         "Day + Food*Genotype")
aic.table %>% arrange(desc(weight))

# Summary of Best Model
summary(daily4)

#Summary of Second Best Model
summary(daily3)

# Tukey Test
marginal = emmeans(daily4, ~ Day|Food)
pairs(marginal, adjust = "tukey")

# Generate Figures
# Lifetime Fecundity Plot
p<-ggplot(aes(y=Fecundity,x=Genotype,fill=Genotype),data=ltfecund) +
  geom_point(pch=21,cex=3, alpha = 0.3)+
  stat_summary(fun.data=mean_se,geom="errorbar",size=1,width=0.35)+
  stat_summary(fun=mean,
               geom="point",cex=5,pch=21,stroke=1.5)+ylim(-1,405)+
  scale_fill_manual(values=c("gray50","gray10")) +
  xlab("Genotype")+ylab("Total Number of Offspring")+facet_wrap(ltfecund$Food)+
  theme_bw() +
  theme(axis.text.x = element_text(color = "grey20", size = 13),
        axis.text.y = element_text(color = "grey20", size = 15),
        strip.text.x = element_text(size = 15),
        legend.position = "none",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(color = "grey20", size = 22),
        axis.title.x = element_text(color = "grey20", size = 22))
p

# Daily Fecundity Plot
pFS<-ggplot(fec_day, aes(x=Day, y=Fecundity,group=Genotype,color=Genotype)) +
  facet_wrap(~Food)+
  scale_colour_manual(values=c("gray50","gray10"))+
  theme_bw()+
  stat_summary(fun = mean,
               geom = "line",size=1.5,aes(group=Genotype,color=Genotype))+
  stat_summary(fun = mean,
               geom = "linerange",color="black",size=1,
               fun.max = function(x) mean(x) +
                 qt(.975, df = length(x)) * sd(x) / sqrt(length(x)),
               fun.min = function(x) mean(x) -
                 qt(.975, df = length(x)) * sd(x) / sqrt(length(x)),
               aes(group=Genotype,color=Genotype)) +
  xlab("Day")+ylab("Daily Number of Offspring") +
  theme(axis.text.x = element_text(color = "grey20", size = 13),
        axis.text.y = element_text(color = "grey20", size = 15),
        strip.text.x = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(color = "grey20", size = 22),
        axis.title.x = element_text(color = "grey20", size = 22))

figure <- ggarrange(p, pFS,
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1,
                    font.label = list(size = 25, color = "black", face = "bold"),
                    common.legend = TRUE, legend = "right")
figure

## Population Growth Assay
# Read in data
pg <- read.csv("pg.csv", header = T)

# Code variables as factors
pg$Food <- factor(pg$Food)
pg$Genotype <- pg$Genotype <- factor(pg$Genotype, levels = c("N2", "ERT250"))
pg$Counter <- factor(pg$Counter)

# Remove incomplete observations
pg <- pg[pg$Died == 0,] # remove worms that died early
pg <- pg[pg$Censored == 0,] # remove censored worms

# summarize the data
pops <- pg %>%
  group_by(ID) %>%
  summarise(ID = first(ID),
            Mean_Count = mean(Count),
            Est_Pop = Mean_Count*725,
            Food = first(Food),
            Genotype = first(Genotype),
            Counter = first(Counter),
            Contaminated = as.factor(first(Contaminated)),
  )

# mean populations
meanpops <- pops %>%
  group_by(Food, Genotype) %>%
  summarise(n = n(),
            Mean_Pop = mean(Est_Pop),
            SEM = sd(Est_Pop)/sqrt(n) 
  )

# Generate Figures
p<-ggplot(aes(y=Est_Pop,x = Genotype, fill = Genotype),data = pops)
p+ geom_point(pch=21,cex=3, alpha = 0.3)+
  geom_signif(comparisons = list(c("N2", "ERT250")), 
              map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "NS" = .5), tip_length = 0) +
  stat_summary(fun.data=mean_se,geom="errorbar",size=1,width=0.35)+
  ylim(c(0,35000))+
  stat_summary(fun=mean,
               geom="point",cex=5,pch=21,stroke=1.5)+
  scale_fill_manual(values=c("gray50","gray10"))+
  xlab("Genotype")+ylab("Population Size")+facet_grid(cols=vars(Food))+
  theme_bw()+ 
  theme(axis.text.x = element_text(color = "grey20", size = 13),
        axis.text.y = element_text(color = "grey20", size = 15),
        strip.text.x = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(color = "grey20", size = 25),
        axis.title.x = element_text(color = "grey20", size = 25))

# Model Selection and Analysis
null <- glmer(Count ~ (1|ID),family = poisson, data = pg)
mod1 <- glmer(Count ~ Genotype*Food + (1|ID),family = poisson, data = pg)
mod2 <- glmer(Count ~ Genotype+Food + (1|ID),family = poisson, data = pg)
mod3 <- glmer(Count ~ Genotype + (1|ID),family = poisson, data = pg)
mod4 <- glmer(Count ~ Food + (1|ID),family = poisson, data = pg)

# Model Comparison
f.aic<-c(AIC(null),AIC(mod1),AIC(mod2), AIC(mod3), AIC(mod4))
delAIC<-f.aic-min(f.aic)
relLik <- exp(-0.5 * delAIC)
aicweight <- relLik/sum(relLik)
aic.table<-data.frame(AIC = f.aic, delAIC = delAIC, weight = aicweight)
rownames(aic.table) <- c("Null","Genotype*Food","Genotype + Food", 
                         "Genotype", "Food")
aic.table # mod1 is the best model

# 61 is an outlier, so let's see if it is driving trends in the best model
no61 <- pg[pg$ID!= 61,]
mod5 <- glmer(Count ~ Genotype*Food + (1|ID),family = poisson, data = no61)
summary(mod5) # nothing changes about the main conclusions, so best to keep it in

# Summary of best model
summary(mod1)

# Pairwise comparisons from best model
marginal = emmeans(mod1, ~ Genotype*Food)
pairs(marginal, adjust = "tukey")

