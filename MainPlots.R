# load dependencies
library(ggplot2)
library(ggpubr)
library(viridis)
library(RColorBrewer)
library(Metrics)
library(scales)
library(tidyr)




################################################
#### FIGURE 2: SIMULATED ANNEALING SCHEDULE ####
################################################
MIN_DIMER=1
DIMER_ADJ=0.1
pal <- RColorBrewer::brewer.pal(5,"Reds")


# ratio-based normalization
#T_INIT = -2*math.log10(MEAN_DIMER/MAX_DIMER)+2
standard <- expand.grid(MaxCost=1:50, MeanCost=c(20,10,5,1,0.5))
standard$T_INIT <- -2*log10(standard$MeanCost/standard$MaxCost)+2
standard$T_INIT[standard$T_INIT<0.3] <- 0.3
p1 <- ggplot(standard[standard$MeanCost>0.5,])+
  geom_line(aes(x=MaxCost, y=T_INIT, col=factor(MeanCost), group=MeanCost), lwd=1.5)+
  xlab(expression(paste("Max (", Delta, "Cost)")))+
  ylab("T_INIT")+#bquote(T[I]))+
  ggtitle("A")+
  scale_color_manual(values=pal)+
  theme_classic()+
  labs(col=expression(paste("Mean (",Delta, "Cost)")))+
  theme(text=element_text(size=14))

#T_INIT = (-2*math.log10(MEAN_DIMER/MAX_DIMER)+1)/10
deltaG <- expand.grid(MaxCost=1:20, MeanCost=c(10,4,2,1,0.2))
deltaG$T_INIT <- (-2*log10(deltaG$MeanCost / deltaG$MaxCost)+1) / 10
deltaG$T_INIT[deltaG$T_INIT<0.03] <- 0.03
ggplot(deltaG)+
  geom_line(aes(x=MaxCost, y=T_INIT, col=factor(MeanCost), group=MeanCost), lwd=2)+
  xlab(expression(paste("Max (",Delta,"Cost)")))+
  ylab("T_INIT")+
  ggtitle("A")+
  theme_classic()+
  labs(col=expression(paste("Mean (",Delta,"Cost)")))+
  scale_color_manual(values=pal)


## Temperature vs. Acceptance probability
prob <- expand.grid(seq(5,0,-0.01), seq(1,9))
names(prob) <- c("TEMP","Dimers")
prob$P <- exp(-2*prob$Dimers/prob$TEMP)
#math.exp(-PROB_ADJ*comparison/Temp)

p3 <- ggplot(prob[prob$Dimers>0,])+
  geom_line(aes(x=TEMP, y=P, col=as.factor(Dimers)), lwd=1.5)+
  scale_color_manual(values=RColorBrewer::brewer.pal(9,"Reds"))+
  theme_classic()+
  ggtitle("C")+
  xlab("Temp")+
  ylab("P(Acceptance)")+
  scale_x_reverse()+
  labs(col=expression(paste(Delta, "Cost")))+
  theme(text=element_text(size=14), legend.position="none")


## Temperature vs. ITERATIONS
T <- expand.grid(iter=1:100, rate=c(0.9,0.95,0.98,0.99))
T$T <- (2-0.1)*(T$rate^T$iter)+0.1

T$rate <- as.factor(T$rate)
p2 <- ggplot(T)+
  geom_line(aes(x=iter, y=T, group=rate, col=rate, lwd=rate))+
  scale_color_manual(values=blues9[c(3,5,7,9)])+
  scale_linewidth_manual(values=c(1,2,1,1))+
  theme_classic()+
  xlab("% Iterations")+
  ylab("Temp")+
  scale_y_continuous(limits=c(0,2))+
  ggtitle("B")+
  labs(col=expression(alpha), lwd=expression(alpha))+
  theme(text=element_text(size=14))#, legend.position="none")


## Acceptance probability vs. ITERATIONS
T$P_1 <- exp(-2*1/T$T)
T$P_2 <- exp(-2*2/T$T)
T$P_3 <- exp(-2*3/T$T)
T$P_4 <- exp(-2*4/T$T)
T$P_5 <- exp(-2*5/T$T)
T$P_6 <- exp(-2*6/T$T)
T$P_7 <- exp(-2*7/T$T)
T$P_8 <- exp(-2*8/T$T)
T$P_9 <- exp(-2*9/T$T)
ts <- tidyr::pivot_longer(T, c("P_1","P_2","P_3","P_4","P_5",
                               "P_6","P_7","P_8","P_9"))
names(ts)[4:5] <- c("Dimers","P")
p4 <- ggplot(ts[ts$rate==0.95,])+
  geom_line(aes(x=iter, y=P, col=Dimers, group=Dimers), lwd=1.5)+
  scale_color_manual(values=RColorBrewer::brewer.pal(9,"Reds"))+
  theme_classic()+
  xlab("% Iterations")+
  ylab("P(Acceptance)")+
  labs(col=expression(paste(Delta,"Cost")))+
  scale_y_log10(labels=scales::label_scientific(digits=2), breaks=c(0,1e-50,1e-20,1e-10,1e-5,1e-1))+
  #theme(legend.position="inside", legend.position.inside=c(0.9,0.7))+
  ggtitle("D")+
  theme(text=element_text(size=14))


fig2ab <- ggpubr::ggarrange(p1+theme(legend.position="bottom"), #+theme(plot.margin=margin(0.5,0.5,0.5,1.5,"cm")), 
                            p2+theme(legend.position="bottom"))
fig2cd_plots <- ggpubr::ggarrange(p3, 
                                  p4+theme(legend.position="none"), nrow=1)
p2cd_leg <- ggpubr::get_legend(p4+theme(legend.position="bottom"))
fig2cd <- ggpubr::ggarrange(fig2cd_plots, p2cd_leg, nrow=2, heights=c(1,0.2))
ggpubr::ggarrange(fig2ab, fig2cd, nrow=1)
ggsave("Figure2_TempPlots.pdf", width=11.5, height=3.4, dpi=300)





#############################################################
#### FIGURE 3: SIMULATED ANNEALING ALGORITHM TRACE PLOTS ####
#############################################################
trace01 <- read.csv("MW_example/3_OptimizedMultiplexes/Trace_Dimer_Load/Run01_costsTrace.csv")
trace02 <- read.csv("MW_example/3_OptimizedMultiplexes/Trace_Dimer_Load/Run02_costsTrace.csv")
trace03 <- read.csv("MW_example/3_OptimizedMultiplexes/Trace_Dimer_Load/Run03_costsTrace.csv")
trace04 <- read.csv("MW_example/3_OptimizedMultiplexes/Trace_Dimer_Load/Run04_costsTrace.csv")
trace05 <- read.csv("MW_example/3_OptimizedMultiplexes/Trace_Dimer_Load/Run05_costsTrace.csv")
trace06 <- read.csv("MW_example/3_OptimizedMultiplexes/Trace_Dimer_Load/Run06_costsTrace.csv")
trace07 <- read.csv("MW_example/3_OptimizedMultiplexes/Trace_Dimer_Load/Run07_costsTrace.csv")
trace08 <- read.csv("MW_example/3_OptimizedMultiplexes/Trace_Dimer_Load/Run08_costsTrace.csv")
trace09 <- read.csv("MW_example/3_OptimizedMultiplexes/Trace_Dimer_Load/Run09_costsTrace.csv")
trace10 <- read.csv("MW_example/3_OptimizedMultiplexes/Trace_Dimer_Load/Run10_costsTrace.csv")

# combine
trace01$Run <- 1
trace02$Run <- 2
trace03$Run <- 3
trace04$Run <- 4
trace05$Run <- 5
trace06$Run <- 6
trace07$Run <- 7
trace08$Run <- 8
trace09$Run <- 9
trace10$Run <- 10
traces <- rbind(trace01, trace02, trace03, trace04, trace05, trace06,
                trace07, trace08, trace09, trace10)
traces$Run <- as.factor(traces$Run)

## example plot showing multiple runs & GLS vs. ASA
p3a <- ggplot(traces)+
  geom_line(aes(x=Iterations, y=TotalDimers, group=Run, col=Run))+
  scale_color_manual(values=pals::kelly(10))+
  theme_classic()+
  ylab("Total Dimers")+
  ylim(0,550)+
  theme(legend.position="none", text=element_text(size=12))+
  geom_bracket(xmin=0, xmax=5000, y.position=500, label="GLS", tip.length=c(0.05,0.05))+
  geom_bracket(xmin=5000,xmax=15000, y.position=220, label="ASA", tip.length=c(0.05,0.05))+
  ggtitle("A")

## example showing hill-climbing for one run (massive improvement)
min1 <- trace08[which(trace08$TotalDimers==min(trace08$TotalDimers)),]
p3b <- ggplot(trace08[trace08$Iterations>5000,])+
  geom_line(aes(x=Iterations, y=TotalDimers),lwd=0.5,col="gray60",lty="dotted")+
  geom_point(aes(x=Iterations, col=ASA_Temp, y=TotalDimers),size=1,alpha=0.8)+
  geom_text(data=min1, aes(x=Iterations, y=TotalDimers, label="★",
                           family="HiraKakuPro-W3"), size=3, col="cyan")+
  geom_bracket(xmin=7000, xmax=8000, y=170, label="cycle", label.size=3, tip.length=c(0.05,0.05))+
  scale_color_viridis(option="inferno",end=0.8)+
  theme_bw()+
  ylab("Total Dimers")+
  labs(col="Temp")+
  ggtitle("B")+
  scale_x_continuous(breaks=seq(5000,15000,1000), 
                     labels=c(5000,rep("",4),10000,rep("",4),15000))+
  theme(#legend.position="inside", legend.position.inside=c(0.8,0.7),
        panel.grid.minor=element_blank())#, axis.text.x=element_text(angle=45,hjust=1))
  
## example showing hill-climbing where not much improvement is possible
trace <- read.csv("MW_example2/3_OptimizedMultiplexes/Trace_Dimer_Load/Run01_costsTrace.csv")
min2 <- trace[which(trace$TotalDimers==min(trace$TotalDimers)),]
p3c <- ggplot(trace[trace$Iterations>5000,])+
  geom_line(aes(x=Iterations, y=TotalDimers),lwd=0.5,col="gray60",lty="dotted")+
  geom_point(aes(x=Iterations, col=ASA_Temp, y=TotalDimers),size=1,alpha=0.8)+
  geom_text(data=min2, aes(x=Iterations, y=TotalDimers, label="★",
                           family="HiraKakuPro-W3"), size=3, col="cyan")+
  #geom_bracket(xmin=7000, xmax=8000, y=56, label="cycle")+
  scale_color_viridis(option="inferno",end=0.8)+
  theme_bw()+
  ylab("Total Dimers")+
  labs(col="Temp")+
  ggtitle("C")+
  scale_x_continuous(breaks=seq(5000,15000,1000), 
                     labels=c(5000,rep("",4),10000,rep("",4),15000))+
  theme(#legend.position="inside", legend.position.inside=c(0.8,0.7),
        panel.grid.minor=element_blank())#, axis.text.x=element_text(angle=45,hjust=1))

fig3 <- ggpubr::ggarrange(p3a+ylab(element_blank()), 
                  ggpubr::ggarrange(p3b+ylab(element_blank())+xlab(element_blank()), 
                                    p3c+ylab(element_blank())+xlab(element_blank()), 
                                    common.legend=TRUE, nrow=1, legend="bottom"), 
                  nrow=2)
ggpubr::annotate_figure(fig3, left=ggpubr::text_grob("Total Dimers", rot=90))
ggsave("Figure3_TracePlots.pdf", width=5, height=5, dpi=300)





#############################################
#### FIGURE 4: ASSESSING PREVIOUS PANELS ####
#############################################
# read in & clean up assessment data
assess <- read.csv("/Users/maggiehallerud/Desktop/testing_SNP_panels/SNP_panel_assessments.csv")
assess$Pairs <- as.numeric(assess$Pairs)
assess$Dimers <- as.numeric(assess$Dimers)
assess$`Dimers_Opt80%` <- as.numeric(assess$`Dimers_Opt80%`)
names(assess)[which(names(assess)=="Optimized?")] <- "Optimized."
assess$Optimized. <- as.factor(assess$Optimized.)
assess <- assess[-which(assess$Animal=="Marten150"),]

# check shape of pattern
assess$Pairs_sq <- assess$Pairs^2
int <- lm(Dimers~1, data=assess)
lin <- lm(Dimers~Pairs, data=assess); rmse(assess$Dimers, int$fitted.values)
quad1 <- lm(Dimers~Pairs_sq+Pairs, data=assess)
log <- lm(log(Dimers+0.1)~Pairs, data=assess)
power <- lm(log(Dimers+0.1)~log(Pairs), data=assess)
rmse(assess$Dimers[!is.na(assess$Dimers)], int$fitted.values)
rmse(assess$Dimers[!is.na(assess$Dimers)], lin$fitted.values)
rmse(assess$Dimers[!is.na(assess$Dimers)], quad1$fitted.values)
rmse(assess$Dimers[!is.na(assess$Dimers)], log$fitted.values)
rmse(assess$Dimers[!is.na(assess$Dimers)], power$fitted.values)
predsA <- predict(quad1, data.frame(Pairs_sq=c(20:400)^2, Pairs=20:400))

## panel size vs. dimer load
p4a <- ggplot()+
  geom_abline(aes(slope=1, intercept=0), col="red")+
  #geom_abline(aes(slope=0.5, intercept=0), col="red", lty="dashed")+
  geom_line(aes(x=20:400, y=predsA), lty="dashed")+
  annotate(geom="label", x=150, y=1150, 
  label=expression('y=0.009x'^2*'-0.5x+28'), size=9, size.unit="pt", border.color="white")+
  #geom_line(data=assess_long, aes(x=Pairs, y=Dimers, group=Animal), lty="dotted", col="gray10")+
  geom_point(data=assess, aes(x=Pairs, y=Dimers, col=Optimized.), alpha=0.8, size=3)+#col=DimersPerPair, 
  scale_color_manual(values=c("gray30","skyblue"))+
  xlab(element_blank())+
  ylab("Predicted dimers")+
  ggtitle("A")+
  theme_bw()+
  theme(legend.position="none", panel.grid.minor=element_blank(), 
        text=element_text(size=11)) #plot.margin=margin(0.1,0,0,0.1,"cm"),
        #axis.title.y=element_text(vjust=-1))
p4a#check

## dimers per pair vs. primer pairs
int <- lm(`Dimers/Total Pairs`~1, assess)
lin <- lm(`Dimers/Total Pairs`~Pairs, assess)
quad <- lm(`Dimers/Total Pairs`~Pairs_sq+Pairs, assess)
log <- lm(log(`Dimers/Total Pairs`+0.1)~Pairs, assess)
pow <- lm(log(`Dimers/Total Pairs`+0.1)~log(Pairs), assess)
rmse(assess$`Dimers/Total Pairs`[!is.na(assess$`Dimers/Total Pairs`)], 
     int$fitted.values)#1.09
rmse(assess$`Dimers/Total Pairs`[!is.na(assess$`Dimers/Total Pairs`)], 
     lin$fitted.values)#0.760
rmse(assess$`Dimers/Total Pairs`[!is.na(assess$`Dimers/Total Pairs`)], 
     quad$fitted.values)#0.760
rmse(assess$`Dimers/Total Pairs`[!is.na(assess$`Dimers/Total Pairs`)], 
     log$fitted.values)#1.59
rmse(assess$`Dimers/Total Pairs`[!is.na(assess$`Dimers/Total Pairs`)],
     power$fitted.values)#3.68
#linpreds <- predict(lin, data.frame(Pairs=1:400))
predsB <- predict(lin, data.frame(Pairs=1:400))

p4b <- ggplot()+
  geom_hline(aes(yintercept=1), col="red")+
  #geom_line(aes(x=1:400, y=linpreds), lty="dotted")+
  geom_line(aes(x=1:400, y=predsB), lty="dashed")+
  geom_point(data=assess, aes(x=Pairs, y=`Dimers/Total Pairs`, col=Optimized.), alpha=0.8, size=3)+
  geom_label(aes(x=130, y=3.9, label="y=0.007x-0.019"), size=9, border.color="white", label.padding=unit(0,"cm"), size.unit="pt")+
  scale_color_manual(values=c("gray30","skyblue"))+
  xlab(element_blank())+
  ylab("Dimers / Pair")+
  ggtitle("B")+
  theme_bw()+
  theme(legend.position="none", panel.grid.minor=element_blank(),
        text=element_text(size=11)) #plot.margin=margin(0.1,0,0,0.1,"cm"),
        #axis.title.y=element_text(vjust=-2))
p4b

## pairs in dimers vs. panel size
int <- lm(`% in Dimers`~Pairs, assess)
lin <- lm(`% in Dimers`~Pairs, assess)
quad <- lm(`% in Dimers`~Pairs_sq+Pairs, assess)
log <- lm(`% in Dimers`~log(Pairs), assess)
pow <- lm(log(`% in Dimers`+0.01)~log(Pairs), assess)
rmse(assess$`% in Dimers`[!is.na(assess$`% in Dimers`)], predict(int))#0.866
rmse(assess$`% in Dimers`[!is.na(assess$`% in Dimers`)], predict(lin))#1.013
rmse(assess$`% in Dimers`[!is.na(assess$`% in Dimers`)], predict(quad))#1.0158
rmse(assess$`% in Dimers`[!is.na(assess$`% in Dimers`)], predict(log))#0.7905
rmse(assess$`% in Dimers`[!is.na(assess$`% in Dimers`)], predict(pow))#1.45
predsC <- predict(log, data.frame(Pairs=1:400))

p4c <- ggplot()+
  geom_line(aes(x=1:400, y=predsC), lty="dashed")+
  geom_point(data=assess, aes(x=Pairs, y=`% in Dimers`, col=Optimized.),size=3,alpha=0.8)+
  geom_label(aes(x=250, y=0.05, label="y=0.30*ln(x)-0.81"), size=9, size.unit="pt", border.color="white")+
  scale_color_manual(values=c("gray30","skyblue"))+
  xlab(element_blank())+
  ylab("Pairs in Dimers")+
  ggtitle("C")+
  theme_bw()+
  theme(legend.position="none", panel.grid.minor=element_blank(),
        text=element_text(size=11))+ #plot.margin=margin(0.1,0,0,0.1,"cm"))+
  scale_y_continuous(labels=scales::percent)+
  coord_cartesian(ylim=c(0,1))
p4c

## "bad" dimers vs panel size
int <- lm(BadDimers~1, assess)
lin <- lm(BadDimers~Pairs, assess)
quad <- lm(BadDimers~Pairs_sq, assess)
log <- lm(log(BadDimers+0.1)~Pairs, assess)
pow <- lm(log(BadDimers+0.1)~log(Pairs), assess)
rmse(assess$BadDimers[!is.na(assess$BadDimers)], predict(int))#19.63
rmse(assess$BadDimers[!is.na(assess$BadDimers)], predict(lin))#14.58
rmse(assess$BadDimers[!is.na(assess$BadDimers)], predict(quad))#14.95
rmse(assess$BadDimers[!is.na(assess$BadDimers)], predict(log))#26.18
rmse(assess$BadDimers[!is.na(assess$BadDimers)], predict(pow))#26.19
predsD <- predict(lin, data.frame(Pairs=1:400))

p4d <- ggplot()+
  geom_line(aes(x=1:400, y=predsD), lty="dashed")+
  geom_point(data=assess, aes(x=Pairs, y=BadDimers, col=Optimized.), size=3, alpha=0.8)+
  geom_label(aes(x=100, y=57, label="y=0.13x-5.7"), size=9, size.unit="pt", border.color="white", label.padding=unit(0,"cm"))+
  scale_color_manual(values=c("gray30","skyblue"))+
  xlab(element_blank())+
  ylab(expression(paste("Dimers with ",Delta,"G < -10")))+
  ggtitle("D")+
  theme_bw()+
  theme(legend.position="none", panel.grid.minor=element_blank(),
        text=element_text(size=11)) #plot.margin=margin(0.1,0,0,0.1,"cm"),
        #axis.title.y=element_text(vjust=-2))
p4d

## mean delta G
int <- lm(MeanDeltaG~1, assess)
lin <- lm(MeanDeltaG~Pairs, assess)
quad <- lm(MeanDeltaG~Pairs_sq+Pairs, assess)
log <- lm(MeanDeltaG~log(Pairs), assess)
pow <- lm(log(MeanDeltaG+0.001)~log(Pairs), assess)
rmse(assess$MeanDeltaG[!is.na(assess$MeanDeltaG)], predict(int))#0.0286
rmse(assess$MeanDeltaG[!is.na(assess$MeanDeltaG)], predict(lin))#0.0286
rmse(assess$MeanDeltaG[!is.na(assess$MeanDeltaG)], predict(quad))#0.0286
rmse(assess$MeanDeltaG[!is.na(assess$MeanDeltaG)], predict(log))#0.0286
rmse(assess$MeanDeltaG[!is.na(assess$MeanDeltaG)], predict(pow))#6.86

p4e <- ggplot()+
  geom_abline(aes(slope=0, intercept=-0.049346), lty="dashed")+
  geom_point(data=assess, aes(x=Pairs, y=MeanDeltaG, col=Optimized.), size=3, alpha=0.8)+
  scale_color_manual(values=c("gray30","skyblue"))+
  geom_label(aes(x=100, y=-0.075, label="y=-0.0493"), size=9, size.unit="pt", border.color="white", label.padding=unit(0,"cm"))+
  xlab(element_blank())+
  ylab(expression(paste("Mean ",Delta,"G")))+
  ggtitle("E")+
  theme_bw()+
  theme(legend.position="none", panel.grid.minor=element_blank(),
        text=element_text(size=11)) #plot.margin=margin(0.1,0,0,0.1,"cm"),
        #axis.text.y=element_text(vjust=-1))
p4e


### demonstrating improvement with 20% panel size reduction
assess_long <- assess %>% 
  tidyr::pivot_longer(cols=c("Dimers","Dimers_Opt80%"), names_to="Type", 
                      values_to="Dimers")
assess_long$Type[assess_long$Type=="Dimers"] <- "Original"
assess_long$Type[assess_long$Type=="Dimers_Opt80%"] <- "20% Reduction"
assess_long$Pairs[assess_long$Type=="20% Reduction"] <- assess_long$`80%Pairs`[assess_long$Type=="20% Reduction"]
p4f <- ggplot()+
  geom_abline(aes(slope=1, intercept=0), col="red")+
  geom_line(data=assess_long, aes(x=Pairs, y=Dimers, group=Animal), lty="dotted", col="gray10")+
  geom_point(data=assess_long, aes(x=Pairs, y=Dimers, shape=Optimized., col=Type, size=Type), alpha=0.8)+
  scale_size_manual(values=c(3,2))+
  scale_color_manual(values=c("blue","gray30"))+
  xlab(element_blank())+
  ylab("Predicted dimers")+
  ggtitle("F")+
  theme_bw()+
  theme(legend.position="none", panel.grid.minor=element_blank(), 
        text=element_text(size=11))# plot.margin=margin(0.1,0,0,0.1,"cm"),
        #axis.title.y=element_text(vjust=-1))
p4f

fig4 <- ggarrange(p4a, p4b, p4c, p4d, p4e, p4f, align="hv", nrow=2, ncol=3)
annotate_figure(fig4, bottom="Primer pairs")#, top=text_grob("Dimer Loads in Existing Panels",face="bold"))
ggsave("Fig4_panel_assessment.pdf", height=4, width=7, dpi=300)





#######################################################
##### FIGURE 6: SIMULATED PANELS BY SIZE / INPUTS #####
#######################################################
# read in simulations data
summ <- read.csv("MW_SIMULATIONS_MERGED.csv")

# calculate medians, range for each species-input-output-algorithm combination
summ_median <- aggregate(TotalDimers~IN_LOCI+Species+N_LOCI+DeltaG, summ, FUN=median)
summ_min <- aggregate(TotalDimers~IN_LOCI+Species+N_LOCI+DeltaG, summ, FUN=function(x) min(x))
summ_max <- aggregate(TotalDimers~IN_LOCI+Species+N_LOCI+DeltaG, summ, FUN=function(x) max(x))
dg_summ_median <- aggregate(MeanDeltaG~IN_LOCI+Species+N_LOCI+DeltaG, summ, FUN=median)
dg_summ_min <- aggregate(MeanDeltaG~IN_LOCI+Species+N_LOCI+DeltaG, summ, FUN=function(x) min(x))
dg_summ_max <- aggregate(MeanDeltaG~IN_LOCI+Species+N_LOCI+DeltaG, summ, FUN=function(x) max(x))
bad_summ_median <- aggregate(N_BadDimers~IN_LOCI+Species+N_LOCI+DeltaG, summ, FUN=median)
bad_summ_min <- aggregate(N_BadDimers~IN_LOCI+Species+N_LOCI+DeltaG, summ, FUN=function(x) min(x))
bad_summ_max <- aggregate(N_BadDimers~IN_LOCI+Species+N_LOCI+DeltaG, summ, FUN=function(x) max(x))

summs <- cbind(summ_median, summ_min$TotalDimers, summ_max$TotalDimers, 
               dg_summ_median$MeanDeltaG, dg_summ_min$MeanDeltaG, dg_summ_max$MeanDeltaG,
               bad_summ_median$N_BadDimers, bad_summ_min$N_BadDimers, bad_summ_max$N_BadDimers)
names(summs)[5:13] <- c("Median","Min","Max","MedianDeltaG","MinDeltaG","MaxDeltaG","MedianBad","MinBad","MaxBad")
summs$Species <- factor(summs$Species, levels=c("Marten","Fisher","GrayFox","SeaOtter"))
summs$N_LOCI <- as.numeric(as.character(summs$N_LOCI))

summs$DeltaG[summs$DeltaG=="True"] <- "Delta G"
summs$DeltaG[summs$DeltaG%in%c("False","FALSE")] <- "Standard"
names(summs)[4] <- "Algorithm"
unique(summs$Algorithm)
summs$Algorithm <- factor(summs$Algorithm, levels=c("Standard","Delta G"))

#jitterer <- position_jitter(width = 10, height=0, seed=123)

## plot for marten
base_plotA <- ggplot(mapping=aes(x=N_LOCI, y=Median, col=Algorithm))+
  xlab(element_blank())+
  ylab("Predicted dimers")+
  ggtitle("A")+
  scale_color_manual(values=c("green","purple"))+
  theme_classic()+
  facet_wrap(~IN_LOCI, nrow=1)+
  theme(panel.grid.minor = element_blank(),
        text=element_text(size=12),
        axis.text.x=element_text(angle=45, hjust=1),
        panel.background=element_rect(fill=scales::alpha("green",0.1)))
  #coord_cartesian(expand=c(0,0))
base_plotB <- ggplot(mapping=aes(x=N_LOCI, y=MedianDeltaG, col=Algorithm))+
  facet_wrap(~IN_LOCI, nrow=1)+
  xlab(element_blank())+
  ylab(expression(paste("Mean ",Delta,"G")))+
  scale_color_manual(values=c("green","purple"))+
  ggtitle("B")+
  theme_classic()+
  theme(panel.grid.minor = element_blank(),
        text=element_text(size=12),
        axis.text.x=element_text(angle=45,hjust=1),
        panel.background = element_rect(fill=scales::alpha("purple",0.1)))
base_plotC <- ggplot(mapping=aes(x=N_LOCI, y=MedianBad, col=Algorithm))+
  facet_wrap(~IN_LOCI, nrow=1)+
  xlab(element_blank())+
  ylab(expression(paste("Dimers with ",Delta,"G< -10")))+
  scale_color_manual(values=c("green","purple"))+
  ggtitle("C")+
  theme_classic()+
  theme(panel.grid.minor = element_blank(),
        text=element_text(size=12),
        axis.text.x=element_text(angle=45,hjust=1),
        panel.background = element_rect(fill=scales::alpha("purple",0.1)))

p5a <- base_plotA+
  geom_hline(data=summs[which(summs$Species=="Marten"),],
             aes(yintercept=0), lty="dashed", lwd=0.3)+
  geom_abline(data=summs[which(summs$Species=="Marten"),], 
              aes(slope=1, intercept=0), col="black", lwd=0.3)+
  geom_errorbar(data=summs[which(summs$Species=="Marten"),], 
                aes(ymin=Min, ymax=Max), alpha=0.7)+
  geom_point(data=summs[which(summs$Species=="Marten"),], 
             alpha=0.7, size=1)
p5b <- base_plotB+
  geom_hline(data=summs[which(summs$Species=="Marten"),],
             aes(yintercept=0), lty="dashed", lwd=0.3)+
  geom_errorbar(data=summs[which(summs$Species=="Marten"),],
                aes(ymin=MinDeltaG, ymax=MaxDeltaG), 
                alpha=0.7)+
  geom_point(data=summs[which(summs$Species=="Marten"),],
             alpha=0.7, size=1)
p5c <- base_plotC+
  geom_hline(data=summs[which(summs$Species=="Marten"),],
             aes(yintercept=0), lty="dashed", lwd=0.3)+
  geom_errorbar(data=summs[which(summs$Species=="Marten"),],
                aes(ymin=MinBad, ymax=MaxBad), 
                alpha=0.7)+
  geom_point(data=summs[which(summs$Species=="Marten"),],
             alpha=0.7, size=1)

fig5 <- ggarrange(p5a, p5b, p5c, common.legend=TRUE, nrow=3, legend="right")
annotate_figure(fig5, bottom="Primer pairs in optimized panel", top="Dimer Loads x Input Templates")
ggsave("Fig5_Marten_simulations.pdf", width=9, height=6, dpi=300)


## plot for fisher
p5a <- base_plotA+
  geom_hline(data=summs[which(summs$Species=="Fisher"),],
             aes(yintercept=0), lty="dashed", lwd=0.3)+
  geom_abline(data=summs[which(summs$Species=="Fisher"),], 
              aes(slope=1, intercept=0), col="black", lwd=0.3)+
  geom_errorbar(data=summs[which(summs$Species=="Fisher"),], 
                aes(ymin=Min, ymax=Max), alpha=0.7)+
  geom_point(data=summs[which(summs$Species=="Fisher"),], 
             alpha=0.7, size=1)
p5b <- base_plotB+
  geom_hline(data=summs[which(summs$Species=="Fisher"),],
             aes(yintercept=0), lty="dashed", lwd=0.3)+
  geom_errorbar(data=summs[which(summs$Species=="Fisher"),],
                aes(ymin=MinDeltaG, ymax=MaxDeltaG), 
                alpha=0.7)+
  geom_point(data=summs[which(summs$Species=="Fisher"),],
             alpha=0.7, size=1)
p5c <- base_plotC+
  geom_hline(data=summs[which(summs$Species=="Fisher"),],
             aes(yintercept=0), lty="dashed", lwd=0.3)+
  geom_errorbar(data=summs[which(summs$Species=="Fisher"),],
                aes(ymin=MinBad, ymax=MaxBad), 
                alpha=0.7)+
  geom_point(data=summs[which(summs$Species=="Fisher"),],
             alpha=0.7, size=1)


fisher <- ggarrange(p5a, p5b, p5c, nrow=3, legend="right", common.legend=TRUE)
annotate_figure(fisher, bottom="Primer pairs in optimized panels",
                top="Dimer Loads x Input Templates")
ggsave("FigS1_Fisher_simulations.pdf", width=9, height=6, dpi=300)


## plot for gray fox
p5a <- base_plotA+
  geom_hline(data=summs[which(summs$Species=="GrayFox"),],
             aes(yintercept=0), lty="dashed", lwd=0.3)+
  geom_abline(data=summs[which(summs$Species=="GrayFox"),],
              aes(slope=1, intercept=0), lwd=0.3)+
  geom_errorbar(data=summs[which(summs$Species=="GrayFox"),], 
                aes(ymin=Min, ymax=Max), alpha=0.7)+
  geom_point(data=summs[which(summs$Species=="GrayFox"),], 
             alpha=0.7, size=1)
p5b <- base_plotB+
  geom_hline(data=summs[which(summs$Species=="GrayFox"),],
             aes(yintercept=0), lty="dashed", lwd=0.3)+
  geom_errorbar(data=summs[which(summs$Species=="GrayFox"),],
                aes(ymin=MinDeltaG, ymax=MaxDeltaG), 
                alpha=0.7)+
  geom_point(data=summs[which(summs$Species=="GrayFox"),],
             alpha=0.7, size=1)
p5c <- base_plotC+
  geom_hline(data=summs[which(summs$Species=="GrayFox"),],
             aes(yintercept=0), lty="dashed", lwd=0.3)+
  geom_errorbar(data=summs[which(summs$Species=="GrayFox"),],
                aes(ymin=MinBad, ymax=MaxBad), 
                alpha=0.7)+
  geom_point(data=summs[which(summs$Species=="GrayFox"),],
             alpha=0.7, size=1)

grayfox <- ggarrange(p5a, p5b, p5c, common.legend=TRUE, legend="right", 
                     nrow=3, align="hv")
annotate_figure(grayfox,  bottom="Primer pairs in optimized panels",
                top="Dimer Loads x Input Templates")
ggarrange("FigS2_GrayFox_simulations.pdf", height=6, width=7, dpi=300)


## plot for sea otter
p5a <- base_plotA+
  geom_hline(data=summs[which(summs$Species=="SeaOtter"),],
             aes(yintercept=0), lty="dashed", lwd=0.3)+
  geom_errorbar(data=summs[which(summs$Species=="SeaOtter"),], 
                aes(ymin=Min, ymax=Max), alpha=0.7)+
  geom_point(data=summs[which(summs$Species=="SeaOtter"),], 
             size=1, alpha=0.7)
p5b <- base_plotB+
  geom_hline(data=summs[which(summs$Species=="SeaOtter"),],
             aes(yintercept=0), lty="dashed", lwd=0.3)+
  geom_errorbar(data=summs[which(summs$Species=="SeaOtter"),],
                aes(ymin=MinDeltaG, ymax=MaxDeltaG), 
                alpha=0.7)+
  geom_point(data=summs[which(summs$Species=="SeaOtter"),],
             alpha=0.7, size=1)
p5c <- base_plotC+
  geom_hline(data=summs[which(summs$Species=="SeaOtter"),],
             aes(yintercept=0), lty="dashed", lwd=0.3)+
  geom_errorbar(data=summs[which(summs$Species=="SeaOtter"),],
                aes(ymin=MinBad, ymax=MaxBad), 
                alpha=0.7)+
  geom_point(data=summs[which(summs$Species=="SeaOtter"),],
             alpha=0.7, size=1)

otter <- ggarrange(p5a, p5b, p5c, common.legend=TRUE, legend="right", nrow=3)
annotate_figure(otter,  bottom="Primer pairs in optimized panels",
                top="Dimer Loads x Input Templates")
ggarrange("SeaOtter_simulations.pdf", height=6, width=6, dpi=300)





##########################################################
##### FIGURE 7: MULTI-LEVEL OPTIMIZATION SIMULATIONS #####
##########################################################
orders <- read.csv("ORDER_SUMMARY_ALL.csv")

# reformat raw summaries
orders$Species <- NA
orders$Species[startsWith(orders$DataSub,"Fisher")] <- "Fisher"
orders$Species[startsWith(orders$DataSub,"Marten")] <- "Marten"
orders$Species[startsWith(orders$DataSub,"GrayFox")] <- "GrayFox"

# split into level1, level2, level3
order1 <- orders[orders$Run==1,]
order2 <- orders[orders$Run==2,]
order3 <- orders[orders$Run==3,]
names(order1)[8] <- c("Run1_Dimers")
names(order2)[8] <- c("Run2_Dimers")
names(order3)[8] <- c("Run3_Dimers")

# make unique identifiers for each dataset-level combo
order3$Order1 <- paste(order3$DataSub, order3$IN1)
order1$Order1 <- paste(order1$DataSub, order1$IN1)
order3$Order2 <- paste(order3$DataSub, order3$IN1, order3$IN2, sep="-")
order2$Order2 <- paste(order2$DataSub, order2$IN1, order2$IN2, sep="-")

dimers_ordersim <- merge(order3[,c("Species","DataSub","Run3_Dimers","Order1","Order2","IN1","IN2","IN3")], 
                         order1[,c("Run1_Dimers","Order1")],
                         by="Order1", all.x=TRUE)
dimers_ordersim <- merge(dimers_ordersim, order2[,c("Run2_Dimers","Order2")], by="Order2", all.x=TRUE)
dimers_ordersim <- dimers_ordersim[,c("Species","DataSub","IN1","IN2","IN3","Run1_Dimers","Run2_Dimers","Run3_Dimers")]
dimers_ordersim$Order <- paste(dimers_ordersim$IN1, dimers_ordersim$IN2, dimers_ordersim$IN3,sep="-")
dimers_ordersim$Order <- factor(dimers_ordersim$Order, levels=c("N20-N50-N1000",
                                                                "N50-N20-N1000",
                                                                "N20-N1000-N50",
                                                                "N50-N1000-N20",
                                                                "N1000-N20-N50",
                                                                "N1000-N50-N20"))

# convert to long format
order_long <- dimers_ordersim %>% pivot_longer(cols=Run1_Dimers:Run3_Dimers, 
                                 names_to="Run", values_to="Dimers")
order_long$Run[order_long$Run=="Run1_Dimers"] <- "1"
order_long$Run[order_long$Run=="Run2_Dimers"] <- "2"
order_long$Run[order_long$Run=="Run3_Dimers"] <- "3"
order_long$Order <- paste0(order_long$IN1, order_long$IN2, order_long$IN3)
unique(order_long$Order)#check - should be 6 of these

# calculate median & range for each level across runs
orders_mean <- aggregate(Dimers~Species+IN1+IN2+IN3+Run+Order, order_long, FUN=median)
orders_min <- aggregate(Dimers~Species+IN1+IN2+IN3+Run+Order, order_long, FUN=min)
orders_max <- aggregate(Dimers~Species+IN1+IN2+IN3+Run+Order, order_long, FUN=max)
names(orders_mean)[7]<-"Mean"; names(orders_min)[7]<-"Min"; names(orders_max)[7]<-"Max"
orders_agg <- cbind(orders_mean, orders_min$Min, orders_max$Max)
names(orders_agg)[8:9] <- c("Min","Max")

orders_agg$Order[orders_agg$Order=="N1000N20N50"] <- "N100-N10-N20"
orders_agg$Order[orders_agg$Order=="N1000N50N20"] <- "N100-N20-N10"
orders_agg$Order[orders_agg$Order=="N20N50N1000"] <- "N10-N20-N100"
orders_agg$Order[orders_agg$Order=="N20N1000N50"] <- "N10-N100-N20"
orders_agg$Order[orders_agg$Order=="N50N1000N20"] <- "N20-N100-N10"
orders_agg$Order[orders_agg$Order=="N50N20N1000"] <- "N20-N10-N100"
orders_agg$Order <- factor(orders_agg$Order, c("N10-N20-N100","N20-N10-N100",
                                               "N10-N100-N20","N20-N100-N10",
                                               "N100-N10-N20","N100-N20-N10"))
ggplot(orders_agg)+
  #geom_abline(aes(slope=1,intercept=0), col="gray", lty=2)+
  geom_line(aes(x=Run, y=Mean, col=Species, group=Species), alpha=0.5)+
  geom_point(aes(x=Run, y=Mean, col=Species), alpha=0.7)+
  geom_linerange(aes(x=Run, ymin=Min, ymax=Max, col=Species), alpha=0.5)+
  facet_wrap(~Order, nrow=3)+
  xlab("Optimization Level")+
  ylab("Dimer Load")+
  theme_classic()+
  ggtitle("Multi-Level Optimization")
ggsave("Fig7_multilevel_optimization.pdf")

