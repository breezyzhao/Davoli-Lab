#analysis HPV-negative TCGA tumors

####################################################################################################
## Aneuploidy score defined based on thersholds of -0.3 and 0.3 and is SCALED (NORMALIZED)
## aneuploidy score : Unweighted.SCNA.0.3.chrom.arm.scale
####################################################################################################


####################################################################################################
## Define deletions in TCGA
#chr9p21.loss using (-0.3) on 9p21.3
#chr3p14.loss using (-0.3) on 3p14.2
#chr17p13.loss using (-0.3) on 17p13.1
#Chrom7.gain using 0.3 on Chrom7
####################################################################################################

data.keep<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/0-others/CCLE/Scotte-3.4-others/Teresa_data/final_use/TCGA.HNSC.alldata-Feb27.2021-Xin.txt",
                      sep="\t",header = T)
#################################################################################################### Fig1
if(fig1=T){
##---Fig1A
library(readxl)
sta.data<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/0-others/CCLE/Scotte-3.4-others/summary_table/MS/Report_of_extra_analyses.txt",
                     sep="\t",header = T)
data1<-read_excel("/Users/zhaox12/Desktop/Teresas_lab/project/0-others/CCLE/Scotte-3.4-others/summary_table/MS/CD_rawdata/EPOC CD3 CD3CD8 CD68 LOH_for dr Davoli_mean score_03Aug2020.xls")
data1<-as.data.frame(data1)
data1$is.LOH3p<-ifelse(data1$Locus_3p_POS==1,"Loss",
                       ifelse(data1$Locus_3p_POS==0,"No_loss",NA))

data1$is.LOH9p<-ifelse(data1$Locus_9p_POS==1,"Loss",
                       ifelse(data1$Locus_9p_POS==0,"No_loss",NA))
data1$is.LOH17p<-ifelse(data1$Locus_17p_POS==1,"Loss",
                        ifelse(data1$Locus_17p_POS==0,"No_loss",NA))
data1$is.POL.chr7<-ifelse(data1$polysomy7==1,"Gain",
                          ifelse(data1$polysomy7==0,"No_Gain",NA))

data1$is.LOH3p <- factor(data1$is.LOH3p, levels=c("No_loss", "Loss"))
data1$is.LOH9p <- factor(data1$is.LOH9p, levels=c("No_loss", "Loss"))
data1$is.LOH17p <- factor(data1$is.LOH17p, levels=c("No_loss", "Loss"))
data1$is.POL.chr7 <- factor(data1$is.POL.chr7, levels=c("No_Gain", "Gain"))

data1$total_CD3_mean<-ifelse(data1$total_CD3_mean==0,0.01,data1$total_CD3_mean)
data1$total_CD3CD8_mean<-ifelse(data1$total_CD3CD8_mean==0,0.01,data1$total_CD3CD8_mean)
data1$total_CD68_mean<-ifelse(data1$total_CD68_mean==0,0.01,data1$total_CD68_mean)
data1$log.cd3.mean<-log(data1$total_CD3_mean,exp(1))
data1$log.cd3cd8.mean<-log(data1$total_CD3CD8_mean,exp(1))
data1$log.cd68.mean<-log(data1$total_CD68_mean,exp(1))
seed(0)
library(ggstatsplot)
library(ggplot2)
library(ggplot2)
source("/Users/zhaox12/Desktop/Teresas_lab/project/cmd/0_others/Scott/specific_ggstats_plot.R")
g1.1<-data1[data1$is.LOH3p=="Loss",]
g1.2<-data1[data1$is.LOH3p=="No_loss",]
g1.test<-wilcox.test(g1.1$total_CD3_mean,g1.2$total_CD3_mean,use = "pairwise.complete.obs")
g1<-ggbetweenstats(
  data = data1,
  x = is.LOH3p,
  y = log.cd3.mean,
  notch = F, # show notched box plot
  mean.plotting = TRUE, # whether mean for each group is to be displayed
  mean.ci = F, # whether to display confidence interval for means
  type = "Parametric", # which type of test is to be run
  k = 2, # number of decimal places for statistical results
  outlier.tagging = F, # whether outliers need to be tagged
  #outlier.label = REACTOME_GAP_JUNCTION_ASSEMBLY, # variable to be used for the outlier tag
  #outlier.label.args = list(color = "darkgreen"), # changing the color for the text label
  xlab = "", # label for the x-axis variable
  ylab = expression("Num of cells (log)/"~mm^2), # label for the y-axis variable
  title = "CD3+ \n3p14_Loss", # title text for the plot
  ggtheme = ggplot2::theme_classic(), # choosing a different theme
  ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
  package = "ggsci", # package from which color palette is to be taken
  palette = "nrc_npg", # choosing a different color palette
  messages = FALSE,
  pairwise.comparisons = F, # display significant pairwise comparisons
  pairwise.annotation = "p.value", # how do you want to annotate the pairwise comparisons
  p.adjust.method = "bonferroni",
  ggplot.component = list(theme(text = element_text(size = 15)))
)+ggplot2::theme(axis.title.y = element_text( size = 20, face = "plain"),
                 axis.text.x =element_text( size = 30, face = "plain"),
                 plot.title = element_text(size=30))+
  ggplot2::annotate("text", x=1.5, y=0.95*10, label= paste0("p=",sta.data$p.Value.from.mixed.effect.models[1]),size=10)+
  scale_color_manual(values=c("#FF0033","#3300FF"))+scale_x_discrete(expand=c(1,0))+coord_cartesian(ylim = c(0, 10))
#y try 0.95

g2.1<-data1[data1$is.LOH3p=="Loss",]
g2.2<-data1[data1$is.LOH3p=="No_loss",]
g2.test<-wilcox.test(g2.1$total_CD3CD8_mean,g2.2$total_CD3CD8_mean,use = "pairwise.complete.obs")
g2<-ggbetweenstats(
  data = data1,
  x = is.LOH3p,
  y = log.cd3cd8.mean,
  notch = F, # show notched box plot
  mean.plotting = TRUE, # whether mean for each group is to be displayed
  mean.ci = F, # whether to display confidence interval for means
  type = "Parametric", # which type of test is to be run
  k = 2, # number of decimal places for statistical results
  outlier.tagging = F, # whether outliers need to be tagged
  #outlier.label = REACTOME_GAP_JUNCTION_ASSEMBLY, # variable to be used for the outlier tag
  #outlier.label.args = list(color = "darkgreen"), # changing the color for the text label
  xlab = "", # label for the x-axis variable
  ylab = expression("Num of cells (log)/"~mm^2), # label for the y-axis variable
  title = "CD8+ \n3p14_Loss", # title text for the plot
  ggtheme = ggplot2::theme_classic(), # choosing a different theme
  ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
  package = "ggsci", # package from which color palette is to be taken
  palette = "nrc_npg", # choosing a different color palette
  messages = FALSE,
  pairwise.comparisons = F, # display significant pairwise comparisons
  pairwise.annotation = "p.value", # how do you want to annotate the pairwise comparisons
  p.adjust.method = "bonferroni",
  ggplot.component = list(theme(text = element_text(size = 15)))
)+ggplot2::theme(axis.title.y = element_text( size = 20, face = "plain"),
                 axis.text.x =element_text( size = 30, face = "plain"),
                 plot.title = element_text(size=30))+
  ggplot2::annotate("text", x=1.5, y=9.5, label= paste0("p=",sta.data$p.Value.from.mixed.effect.models[9]),size=10)+
  scale_color_manual(values=c("#FF0033","#3300FF"))+scale_x_discrete(expand=c(1,0))+coord_cartesian(ylim = c(0, 10))

g3.1<-data1[data1$is.LOH9p=="Loss",]
g3.2<-data1[data1$is.LOH9p=="No_loss",]
g3.test<-wilcox.test(g3.1$total_CD3_mean,g3.2$total_CD3_mean,use = "pairwise.complete.obs")
g3<-ggbetweenstats(
  data = data1,
  x = is.LOH9p,
  y = log.cd3.mean,
  notch = F, # show notched box plot
  mean.plotting = TRUE, # whether mean for each group is to be displayed
  mean.ci = F, # whether to display confidence interval for means
  type = "Parametric", # which type of test is to be run
  k = 2, # number of decimal places for statistical results
  outlier.tagging = F, # whether outliers need to be tagged
  #outlier.label = REACTOME_GAP_JUNCTION_ASSEMBLY, # variable to be used for the outlier tag
  #outlier.label.args = list(color = "darkgreen"), # changing the color for the text label
  xlab = "", # label for the x-axis variable
  ylab = expression("Num of cells (log)/"~mm^2), # label for the y-axis variable
  title = "9p21.3_Loss", # title text for the plot
  ggtheme = ggplot2::theme_classic(), # choosing a different theme
  ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
  package = "ggsci", # package from which color palette is to be taken
  palette = "nrc_npg", # choosing a different color palette
  messages = FALSE,
  pairwise.comparisons = F, # display significant pairwise comparisons
  pairwise.annotation = "p.value", # how do you want to annotate the pairwise comparisons
  p.adjust.method = "bonferroni",
  ggplot.component = list(theme(text = element_text(size = 15)))
)+ggplot2::theme(axis.title.y = element_text( size = 20, face = "plain"),
                 axis.text.x =element_text( size = 30, face = "plain"),
                 plot.title = element_text(size=30))+
  ggplot2::annotate("text", x=1.5, y=0.95*10, label= paste0("p=",sta.data$p.Value.from.mixed.effect.models[3]),size=10)+
  scale_color_manual(values=c("#FF0033","#3300FF"))+scale_x_discrete(expand=c(1,0))+coord_cartesian(ylim = c(0, 10))


g4.1<-data1[data1$is.LOH9p=="Loss",]
g4.2<-data1[data1$is.LOH9p=="No_loss",]
g4.test<-wilcox.test(g4.1$total_CD3CD8_mean,g4.2$total_CD3CD8_mean,use = "pairwise.complete.obs")
g4<-ggbetweenstats(
  data = data1,
  x = is.LOH9p,
  y = log.cd3cd8.mean,
  notch = F, # show notched box plot
  mean.plotting = TRUE, # whether mean for each group is to be displayed
  mean.ci = F, # whether to display confidence interval for means
  type = "Parametric", # which type of test is to be run
  k = 2, # number of decimal places for statistical results
  outlier.tagging = F, # whether outliers need to be tagged
  #outlier.label = REACTOME_GAP_JUNCTION_ASSEMBLY, # variable to be used for the outlier tag
  #outlier.label.args = list(color = "darkgreen"), # changing the color for the text label
  xlab = "", # label for the x-axis variable
  ylab = expression("Num of cells (log)/"~mm^2), # label for the y-axis variable
  title = "9p21.3_Loss", # title text for the plot
  ggtheme = ggplot2::theme_classic(), # choosing a different theme
  ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
  package = "ggsci", # package from which color palette is to be taken
  palette = "nrc_npg", # choosing a different color palette
  messages = FALSE,
  pairwise.comparisons = F, # display significant pairwise comparisons
  pairwise.annotation = "p.value", # how do you want to annotate the pairwise comparisons
  p.adjust.method = "bonferroni",
  ggplot.component = list(theme(text = element_text(size = 15)))
)+ggplot2::theme(axis.title.y = element_text( size = 20, face = "plain"),
                 axis.text.x =element_text( size = 30, face = "plain"),
                 plot.title = element_text(size=30))+
  ggplot2::annotate("text", x=1.5, y=9.5, label= paste0("p=",sta.data$p.Value.from.mixed.effect.models[11]),size=10)+
  scale_color_manual(values=c("#FF0033","#3300FF"))+scale_x_discrete(expand=c(1,0))+coord_cartesian(ylim = c(0, 10))

g5.1<-data1[data1$is.LOH17p=="Loss",]
g5.2<-data1[data1$is.LOH17p=="No_loss",]
g5.test<-wilcox.test(g5.1$total_CD3_mean,g5.2$total_CD3_mean,use = "pairwise.complete.obs")
g5<-ggbetweenstats(
  data = data1,
  x = is.LOH17p,
  y = log.cd3.mean,
  notch = F, # show notched box plot
  mean.plotting = TRUE, # whether mean for each group is to be displayed
  mean.ci = F, # whether to display confidence interval for means
  type = "Parametric", # which type of test is to be run
  k = 2, # number of decimal places for statistical results
  outlier.tagging = F, # whether outliers need to be tagged
  #outlier.label = REACTOME_GAP_JUNCTION_ASSEMBLY, # variable to be used for the outlier tag
  #outlier.label.args = list(color = "darkgreen"), # changing the color for the text label
  xlab = "", # label for the x-axis variable
  ylab = expression("Num of cells (log)/"~mm^2), # label for the y-axis variable
  title = "17p13.1_Loss", # title text for the plot
  ggtheme = ggplot2::theme_classic(), # choosing a different theme
  ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
  package = "ggsci", # package from which color palette is to be taken
  palette = "nrc_npg", # choosing a different color palette
  messages = FALSE,
  pairwise.comparisons = F, # display significant pairwise comparisons
  pairwise.annotation = "p.value", # how do you want to annotate the pairwise comparisons
  p.adjust.method = "bonferroni",
  ggplot.component = list(theme(text = element_text(size = 15)))
)+ggplot2::theme(axis.title.y = element_text( size = 20, face = "plain"),
                 axis.text.x =element_text( size = 30, face = "plain"),
                 plot.title = element_text(size=30))+
  ggplot2::annotate("text", x=1.5, y=0.95*10, label= paste0("p=",sta.data$p.Value.from.mixed.effect.models[5]),size=10)+
  scale_color_manual(values=c("#FF0033","#3300FF"))+scale_x_discrete(expand=c(1,0))+coord_cartesian(ylim = c(0, 10))


g6.1<-data1[data1$is.LOH17p=="Loss",]
g6.2<-data1[data1$is.LOH17p=="No_loss",]
g6.test<-wilcox.test(g6.1$total_CD3CD8_mean,g6.2$total_CD3CD8_mean,use = "pairwise.complete.obs")
g6<-ggbetweenstats(
  data = data1,
  x = is.LOH17p,
  y = log.cd3cd8.mean,
  notch = F, # show notched box plot
  mean.plotting = TRUE, # whether mean for each group is to be displayed
  mean.ci = F, # whether to display confidence interval for means
  type = "Parametric", # which type of test is to be run
  k = 2, # number of decimal places for statistical results
  outlier.tagging = F, # whether outliers need to be tagged
  #outlier.label = REACTOME_GAP_JUNCTION_ASSEMBLY, # variable to be used for the outlier tag
  #outlier.label.args = list(color = "darkgreen"), # changing the color for the text label
  xlab = "", # label for the x-axis variable
  ylab = expression("Num of cells (log)/"~mm^2), # label for the y-axis variable
  title = "17p13.1_Loss", # title text for the plot
  ggtheme = ggplot2::theme_classic(), # choosing a different theme
  ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
  package = "ggsci", # package from which color palette is to be taken
  palette = "nrc_npg", # choosing a different color palette
  messages = FALSE,
  pairwise.comparisons = F, # display significant pairwise comparisons
  pairwise.annotation = "p.value", # how do you want to annotate the pairwise comparisons
  p.adjust.method = "bonferroni",
  ggplot.component = list(theme(text = element_text(size = 15)))
)+ggplot2::theme(axis.title.y = element_text( size = 20, face = "plain"),
                 axis.text.x =element_text( size = 30, face = "plain"),
                 plot.title = element_text(size=30))+
  ggplot2::annotate("text", x=1.5, y=9.5, label= paste0("p=",sta.data$p.Value.from.mixed.effect.models[13]),size=10)+
  scale_color_manual(values=c("#FF0033","#3300FF"))+scale_x_discrete(expand=c(1,0))+coord_cartesian(ylim = c(0, 10))

g7.1<-data1[data1$is.POL.chr7=="Gain",]
g7.2<-data1[data1$is.POL.chr7=="No_Gain",]
g7.test<-wilcox.test(g7.1$total_CD3_mean,g7.2$total_CD3_mean,use = "pairwise.complete.obs")
g7<-ggbetweenstats(
  data = data1,
  x = is.POL.chr7,
  y = log.cd3.mean,
  notch = F, # show notched box plot
  mean.plotting = TRUE, # whether mean for each group is to be displayed
  mean.ci = F, # whether to display confidence interval for means
  type = "Parametric", # which type of test is to be run
  k = 2, # number of decimal places for statistical results
  outlier.tagging = F, # whether outliers need to be tagged
  #outlier.label = REACTOME_GAP_JUNCTION_ASSEMBLY, # variable to be used for the outlier tag
  #outlier.label.args = list(color = "darkgreen"), # changing the color for the text label
  xlab = "", # label for the x-axis variable
  ylab = expression("Num of cells (log)/"~mm^2), # label for the y-axis variable
  title = "Chr7 gain", # title text for the plot
  ggtheme = ggplot2::theme_classic(), # choosing a different theme
  ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
  package = "ggsci", # package from which color palette is to be taken
  palette = "nrc_npg", # choosing a different color palette
  messages = FALSE,
  pairwise.comparisons = F, # display significant pairwise comparisons
  pairwise.annotation = "p.value", # how do you want to annotate the pairwise comparisons
  p.adjust.method = "bonferroni",
  ggplot.component = list(theme(text = element_text(size = 15)))
)+ggplot2::theme(axis.title.y = element_text( size = 20, face = "plain"),
                 axis.text.x =element_text( size = 30, face = "plain"),
                 plot.title = element_text(size=30))+
  ggplot2::annotate("text", x=1.5, y=0.95*10, label= paste0("p=",sta.data$p.Value.from.mixed.effect.models[7]),size=10)+
  scale_color_manual(values=c("#FF0033","#3300FF"))+scale_x_discrete(expand=c(1,0))+coord_cartesian(ylim = c(0, 10))


g8.1<-data1[data1$is.POL.chr7=="Gain",]
g8.2<-data1[data1$is.POL.chr7=="No_Gain",]
g8.test<-wilcox.test(g8.1$total_CD3CD8_mean,g8.2$total_CD3CD8_mean,use = "pairwise.complete.obs")
g8<-ggbetweenstats(
  data = data1,
  x = is.POL.chr7,
  y = log.cd3cd8.mean,
  notch = F, # show notched box plot
  mean.plotting = TRUE, # whether mean for each group is to be displayed
  mean.ci = F, # whether to display confidence interval for means
  type = "Parametric", # which type of test is to be run
  k = 2, # number of decimal places for statistical results
  outlier.tagging = F, # whether outliers need to be tagged
  #outlier.label = REACTOME_GAP_JUNCTION_ASSEMBLY, # variable to be used for the outlier tag
  #outlier.label.args = list(color = "darkgreen"), # changing the color for the text label
  xlab = "", # label for the x-axis variable
  ylab = expression("Num of cells (log)/"~mm^2), # label for the y-axis variable
  title = "Chr7 gain", # title text for the plot
  ggtheme = ggplot2::theme_classic(), # choosing a different theme
  ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
  package = "ggsci", # package from which color palette is to be taken
  palette = "nrc_npg", # choosing a different color palette
  messages = FALSE,
  pairwise.comparisons = F, # display significant pairwise comparisons
  pairwise.annotation = "p.value", # how do you want to annotate the pairwise comparisons
  p.adjust.method = "bonferroni",
  ggplot.component = list(theme(text = element_text(size = 15)))
)+ggplot2::theme(axis.title.y = element_text( size = 20, face = "plain"),
                 axis.text.x =element_text( size = 30, face = "plain"),
                 plot.title = element_text(size=30))+
  ggplot2::annotate("text", x=1.5, y=9.5, label= paste0("p=",sta.data$p.Value.from.mixed.effect.models[15]),size=10)+
  scale_color_manual(values=c("#FF0033","#3300FF"))+scale_x_discrete(expand=c(1,0))+coord_cartesian(ylim = c(0, 10))

library(ggpubr)
figure<-ggarrange(g1,g2,g3,g4,g5,g6,
                  labels = c("", "", "","","",""),
                  ncol = 2, nrow = 3,
                  font.label = list(size = 30, color = "black"))
annotate_figure(figure,
                top = text_grob("", color = "black", face = "bold", size = 35))
##---Fig1B
library(ggstatsplot)
library(ggplot2)
library(ggpubr)
library(lmerTest)
source("/Users/zhaox12/Desktop/Teresas_lab/project/cmd/0_others/Scott/specific_ggstats_plot.R")
library(readxl)
set.seed(0)
data.clin<-read_excel("/Users/zhaox12/Desktop/Teresas_lab/project/0-others/CCLE/Scotte-3.4-others/summary_table/MS/CD_rawdata/EPOC CD3 CD3CD8 CD68 LOH_for dr Davoli_mean score_03Aug2020.xls")
data.clin<-as.data.frame(data.clin)
data.clin$sum<-apply(data.clin[,2:4],1,sum) #not include chr7 gain
data.clin$score<-ifelse(data.clin$sum>0,"SCNA","No_SCNA")
data.clin$total_CD3_mean<-ifelse(data.clin$total_CD3_mean==0,0.01,data.clin$total_CD3_mean)
data.clin$total_CD3CD8_mean<-ifelse(data.clin$total_CD3CD8_mean==0,0.01,data.clin$total_CD3CD8_mean)
data.clin$total_CD68_mean<-ifelse(data.clin$total_CD68_mean==0,0.01,data.clin$total_CD68_mean)
data.clin$log.cd3.mean<-log(data.clin$total_CD3_mean,exp(1))
data.clin$log.cd3cd8.mean<-log(data.clin$total_CD3CD8_mean,exp(1))
data.clin$log.cd68.mean<-log(data.clin$total_CD68_mean,exp(1))
D<-data.clin
D$score <- factor(D$score, levels=c("No_SCNA", "SCNA"))
data.mix<-read_excel("/Users/zhaox12/Desktop/Teresas_lab/project/0-others/CCLE/Scotte-3.4-others/summary_table/MS/CD_rawdata/EPOC CD3 CD3CD8 LOH_for dr Davoli_multi_records_03Aug2020.xls")
data.mix$total_CD3<-ifelse(data.mix$total_CD3==0,0.01,data.mix$total_CD3)
data.mix$total_CD3CD8<-ifelse(data.mix$total_CD3CD8==0,0.01,data.mix$total_CD3CD8)
data.mix$total_CD68<-ifelse(data.mix$total_CD68==0,0.01,data.mix$total_CD68)
data.mix$total_CD3<-log(data.mix$total_CD3,exp(1))
data.mix$total_CD3CD8<-log(data.mix$total_CD3CD8,exp(1))
data.mix$total_CD68<-log(data.mix$total_CD68,exp(1))
data.mix$sum<-apply(data.mix[,2:4],1,sum) #no
data.mix$score<-ifelse(data.mix$sum>0,1,0)
data.cd3<-summary(lmer(total_CD3~score+(1|rand_num) , data=data.mix))
data.cd3cd8<-summary(lmer(total_CD3CD8~score+(1|rand_num) , data=data.mix))
data.cd68<-summary(lmer(total_CD68~score+(1|rand_num) , data=data.mix))
g1<-ggbetweenstats(
  data = D,
  x = score,
  y = log.cd3.mean,
  notch = F, # show notched box plot
  mean.plotting = TRUE, # whether mean for each group is to be displayed
  mean.ci = F, # whether to display confidence interval for means
  type = "Parametric", # which type of test is to be run
  k = 2, # number of decimal places for statistical results
  outlier.tagging = F, # whether outliers need to be tagged
  xlab = "", # label for the x-axis variable
  ylab = expression("Num of cells (log)/"~mm^2), # label for the y-axis variable
  title = "CD3+", # title text for the plot
  ggtheme = ggplot2::theme_classic(), # choosing a different theme
  ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
  messages = FALSE,
  pairwise.comparisons = F, # display significant pairwise comparisons
  ggplot.component = list(theme(text = element_text(size = 15)))
)+ggplot2::theme(axis.title.y = element_text( size = 15, face = "plain"),
                 axis.text.x =element_text( size = 20, face = "plain"),
                 plot.title = element_text(size=23))+
  ggplot2::annotate("text", x=1.5, y=0.95*11, label= paste0("p=",signif(data.cd3$coefficients[2,5],3)),size=8)+
  scale_color_manual(values=c("#FF0033","#3300FF"))+scale_x_discrete(expand=c(1.5,0))+coord_cartesian(ylim = c(0, 11))
# t2<-wilcox.test(t1.1$log.cd3cd8.mean,t1.2$log.cd3cd8.mean)
g2<-ggbetweenstats(
  data = D,
  x = score,
  y = log.cd3cd8.mean,
  notch = F, # show notched box plot
  mean.plotting = TRUE, # whether mean for each group is to be displayed
  mean.ci = F, # whether to display confidence interval for means
  type = "Parametric", # which type of test is to be run
  k = 2, # number of decimal places for statistical results
  outlier.tagging = F, # whether outliers need to be tagged
  xlab = "", # label for the x-axis variable
  ylab = expression("Num of cells (log)/"~mm^2), # label for the y-axis variable
  title = "CD8+", # title text for the plot
  ggtheme = ggplot2::theme_classic(), # choosing a different theme
  ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
  messages = FALSE,
  pairwise.comparisons = F, # display significant pairwise comparisons
  ggplot.component = list(theme(text = element_text(size = 15)))
)+ggplot2::theme(axis.title.y = element_text( size = 15, face = "plain"),
                 axis.text.x =element_text( size = 20, face = "plain"),
                 plot.title = element_text(size=23))+
  ggplot2::annotate("text", x=1.5, y=0.95*11, label= paste0("p=",signif(data.cd3cd8$coefficients[2,5],3)),size=8)+
  scale_color_manual(values=c("#FF0033","#3300FF"))+scale_x_discrete(expand=c(1.5,0))+coord_cartesian(ylim = c(0, 11))
g3<-ggbetweenstats(
  data = D,
  x = score,
  y = log.cd68.mean,
  notch = F, # show notched box plot
  mean.plotting = TRUE, # whether mean for each group is to be displayed
  mean.ci = F, # whether to display confidence interval for means
  type = "Parametric", # which type of test is to be run
  k = 2, # number of decimal places for statistical results
  outlier.tagging = F, # whether outliers need to be tagged
  xlab = "", # label for the x-axis variable
  ylab = expression("Num of cells (log)/"~mm^2), # label for the y-axis variable
  title = "CD68+", # title text for the plot
  ggtheme = ggplot2::theme_classic(), # choosing a different theme
  ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
  messages = FALSE,
  pairwise.comparisons = F, # display significant pairwise comparisons
  ggplot.component = list(theme(text = element_text(size = 15)))
)+ggplot2::theme(axis.title.y = element_text( size = 15, face = "plain"),
                 axis.text.x =element_text( size = 20, face = "plain"),
                 plot.title = element_text(size=23))+
  ggplot2::annotate("text", x=1.5, y=0.95*11, label= paste0("p=",signif(data.cd68$coefficients[2,5],3)),size=8)+
  scale_color_manual(values=c("#FF0033","#3300FF"))+scale_x_discrete(expand=c(1.5,0))+coord_cartesian(ylim = c(0, 11))
figure<-ggarrange(g1,g2,g3,
                  labels = c("",""),
                  ncol = 1, nrow = 3,
                  font.label = list(size = 30, color = "black"))
annotate_figure(figure,
                top = text_grob("", color = "black", face = "bold", size = 24))
##---Fig 1C
library(readxl)
library(survival)
data.clin<-read_excel("/Users/zhaox12/Desktop/Teresas_lab/project/0-others/CCLE/Scotte-3.4-others/summary_table/MS/CD_rawdata/EPOC data for dr Davoli_09Oct2020_dealed.xls")
data.clin1<-as.data.frame(data.clin)
loh.data<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/0-others/CCLE/Scotte-3.4-others/summary_table/MS/CD_rawdata/dealed/LOH-transfom-transpose.txt",
                     sep="\t",header = T)
loh.use1<-loh.data
loh.use1$D3S1234.scale<-scale(loh.use1$D3S1234)
loh.use1$D9S1748.scale<-scale(loh.use1$D9S1748)
loh.use1$D17S786.scale<-scale(loh.use1$D17S786)
loh.use1$IFNA.scale<-scale(loh.use1$IFNA)
loh.use1$TP53.scale<-scale(loh.use1$TP53)
loh.use1[,2:20]<-apply(loh.use1[,2:20],2,function(x) ifelse(x>1.43 | x<0.7,1,0)) 
loh.use1$LOH.3p<-ifelse(loh.use1$D3S1228==1 & loh.use1$D3S1234==1 & loh.use1$D3S1300==1,1,0)
loh.use1$LOH.9p<-ifelse(loh.use1$D9S171==1 & loh.use1$D9S1748==1 & loh.use1$D9S1751==1 & loh.use1$IFNA==1,1,0)
table(loh.use1$LOH.9p)
loh.use1$LOH.17p<-ifelse(loh.use1$TP53==1 & loh.use1$CHRNB1==1 & loh.use1$D17S786==1,1,0)
loh.use1$LOH.11p<-ifelse(loh.use1$D11S1778==1 & loh.use1$INT2==1,1,0)
loh.use1$LOH.13p<-ifelse(loh.use1$D13S133==1 & loh.use1$D13S170 ==1,1,0)
loh.use1$LOH.4p<-ifelse(loh.use1$D4S243==1 & loh.use1$FABP2==1,1,0)
loh.use1$LOH.8p<-ifelse(loh.use1$D8S261==1 & loh.use1$D8S262==1 & loh.use1$D8S264==1 ,1,0)
merge.clin<-merge(data.clin1,loh.use1,by.x="Alt_id",by.y="sample")
D<-merge.clin
D$FISH.test<-ifelse(D$FISH=="DS",0,1)
D$Histology<-ifelse(D$HISTBASELINER3_2=="1:Dysplasia",1,0)
D$loss.events<-apply(D[,c("LOH.3p","LOH.9p","LOH.17p","LOH.4p","LOH.11p","LOH.13p","LOH.8p")],1,sum)
D$major.loss.events<-apply(D[,c("LOH.3p","LOH.9p","LOH.17p")],1,sum)
D$minor.loss.events<-apply(D[,c("LOH.4p","LOH.11p","LOH.13p","LOH.8p")],1,sum)
D$SCNA.events<-ifelse(D$loss.events>0 & D$polysomy7>0,2,
                      ifelse(D$loss.events>0 | D$polysomy7>0,1,0))
D$major.loss<-ifelse(D$major.loss.events>0,1,0)
D$minor.loss<-ifelse(D$minor.loss.events>0,1,0)

D$cd3.group<-ifelse(D$total_CD3_mean>quantile(D$total_CD3_mean,0.65,na.rm=T),1,
                    ifelse(D$total_CD3_mean<quantile(D$total_CD3_mean,0.35,na.rm=T),0,NA))
D$cd3cd8.group<-ifelse(D$total_CD3CD8_mean>quantile(D$total_CD3CD8_mean,0.65,na.rm=T),1,
                       ifelse(D$total_CD3CD8_mean<quantile(D$total_CD3CD8_mean,0.35,na.rm=T),0,NA))
D$cd68.group<-ifelse(D$total_CD68_mean>quantile(D$total_CD68_mean,0.65,na.rm=T),1,
                     ifelse(D$total_CD68_mean<quantile(D$total_CD68_mean,0.35,na.rm=T),0,NA))
D$pdl1.group<-ifelse(D$PDL1_pct_mean>=quantile(D$PDL1_pct_mean,0.65,na.rm=T),1,
                     ifelse(D$PDL1_pct_mean<=quantile(D$PDL1_pct_mean,0.35,na.rm=T),0,NA))

D.cd3.use<-D[!is.na(D$cd3.group),]
D.cd3cd8.use<-D[!is.na(D$cd3cd8.group),]
D.cd68.use<-D[!is.na(D$cd68.group),]
D.pdl1.use<-D[!is.na(D$pdl1.group),]
#binary
logit.cd3 <- glm(cd3.group ~ Locus_3p_POS + Locus_9p_POS + Locus_17p_POS+ SCNA.events , data=D.cd3.use, family=binomial(link = "logit"))
summary(logit.cd3)$coefficient
logit.cd3cd8 <- glm(cd3cd8.group ~ Locus_3p_POS + Locus_9p_POS + Locus_17p_POS + SCNA.events, data=D.cd3cd8.use, family=binomial(link="logit"))
summary(logit.cd3cd8)$coefficient
logit.cd68 <- glm(cd68.group ~ Locus_3p_POS + Locus_9p_POS + Locus_17p_POS + SCNA.events, data=D.cd68.use, family=binomial(link="logit"))
summary(logit.cd68)$coefficient
#continuous
logit.cd3 <- glm(cd3.group ~ D3S1234.scale + D9S1748.scale + D17S786.scale  +SCNA.events , data=D.cd3.use, family=binomial(link = "logit"))
summary(logit.cd3)$coefficient
logit.cd3cd8 <- glm(cd3cd8.group ~ D3S1234.scale + D9S1748.scale + D17S786.scale  +SCNA.events, data=D.cd3cd8.use, family=binomial(link="logit"))
summary(logit.cd3cd8)$coefficient
logit.cd68 <- glm(cd68.group ~ D3S1234.scale + D9S1748.scale + D17S786.scale  +SCNA.events, data=D.cd68.use, family=binomial(link="logit"))
summary(logit.cd68)$coefficient
}
#################################################################################################### Fig2
if(fig2==T){
 ##---Fig2A and B
library("survival")
library("survminer")
library(readxl)
data<-read_excel("/Users/zhaox12/Desktop/Teresas_lab/project/0-others/CCLE/Scotte-3.4-others/summary_table/MS/CD_rawdata/EPOC data for dr Davoli_09Oct2020.xls")
D<-as.data.frame(data)
D$total_CD3_mean<-ifelse(D$total_CD3_mean==0,0.01,D$total_CD3_mean)
D$total_CD3CD8_mean<-ifelse(D$total_CD3CD8_mean==0,0.01,D$total_CD3CD8_mean)
D$total_CD68_mean<-ifelse(D$total_CD68_mean==0,0.01,D$total_CD68_mean)
D$log.cd3.mean<-log(D$total_CD3_mean,exp(1))
D$log.cd3cd8.mean<-log(D$total_CD3CD8_mean,exp(1))
D$log.cd68.mean<-log(D$total_CD68_mean,exp(1))
D$CD3.level<-ifelse(D$log.cd3.mean>quantile(D$log.cd3.mean,0.5,na.rm=T),1,0)
D$group<-ifelse(D$HISTBASELINER3_2=="1:Dysplasia" & D$Locus_9p_POS==1,"Dysplasia and 9p21.3 loss","others")
D$group.non9p.dys<-ifelse(D$HISTBASELINER3_2=="1:Dysplasia" & D$Locus_9p_POS==1,"Dysplasia and 9p21.3 loss",
                          ifelse(D$HISTBASELINER3_2=="1:Dysplasia" & D$Locus_9p_POS==0,"Dysplasia only",NA))
D$group.9p.nondys<-ifelse(D$HISTBASELINER3_2=="1:Dysplasia" & D$Locus_9p_POS==1,"Dysplasia and 9p21.3 loss",
                          ifelse(D$HISTBASELINER3_2=="0:hyperkeratosis and hyperplasia" & D$Locus_9p_POS==1,"9p21.3 loss",NA))
D$gender.binary<-ifelse(D$gender=="M",1,0)
D$smoke.cat<-ifelse(D$Smk=="Never",0,
                    ifelse(D$Smk=="Current" | D$Smk=="Former",1,"NA"))
D$Alcohol.cat<-ifelse(D$Alcohol=="No",0,
                      ifelse(D$Alcohol=="Yes",1,"NA"))
D$Histology<-ifelse(D$HISTBASELINER3_2=="1:Dysplasia",1,0)
D$Arm9p_only_loss<-ifelse(D$Locus_9p_POS==1,1,
                          ifelse(D$Locus_3p_POS==0 & (D$Locus_9p_POS==0 | is.na(D$Locus_9p_POS)) & D$Locus_17p_POS==0,0,NA))
D$Chr9p21.3_loss<-ifelse(D$Locus_9p_POS==1,"Loss","No_loss")
D$Group<-ifelse(D$Locus_9p_POS==1,"9p21.3Loss",
                ifelse((D$Locus_3p_POS==1 |D$Locus_17p_POS==1) & D$Locus_9p_POS==0, "3p14 loss or 17p13.1 loss",
                       ifelse(D$Locus_3p_POS==0 & D$Locus_17p_POS==0 & D$Locus_9p_POS==0,"No_loss",NA)))

D$Arm3p17p_only_loss<-ifelse(D$Locus_3p_POS==1 | D$Locus_17p_POS==1,1,0)


fit1 <- survfit(Surv(time_CFS ,event_CFS) ~ group , data = D)
fit2 <- survfit(Surv(time_CFS ,event_CFS) ~ Group, data = D, na.action = na.omit)
g1<-ggsurvplot(fit = fit1,
               pval = TRUE, conf.int = F,
               risk.table = F, # Add risk table
               risk.table.col = "", # Change risk table color by groups
               linetype = "strata", # Change line type by groups
               ggtheme = theme_bw(), # Change ggplot2 theme
               palette = c("#FF0033","#3300FF"),
               title=(""),
               font.legend = list(size = 15, color = "black"),
               legend = c(0.3, 0.3),
               pval.coord = c(0, 0.03))
ggpar(g1, 
      font.main = c(15, "bold"),
      font.x = c(15, "bold"),
      font.y = c(15, "bold"),
      font.caption = c(15, "bold"), 
      font.legend = c(15, "bold"), 
      font.tickslab = c(15, "bold"))
g2<-ggsurvplot(fit = fit2,
               pval = TRUE, conf.int = F,
               risk.table = F, # Add risk table
               risk.table.col = "", # Change risk table color by groups
               linetype = "strata", # Change line type by groups
               ggtheme = theme_bw(), # Change ggplot2 theme
               palette = c("green","#FF0033","#3300FF"),
               title=(""),
               font.legend = list(size = 15, color = "black"),
               legend = c(0.3, 0.3),
               pval.coord = c(0, 0.03))
ggpar(g2, 
      font.main = c(15, "bold"),
      font.x = c(15, "bold"),
      font.y = c(15, "bold"),
      font.caption = c(15, "bold"), 
      font.legend = c(15, "bold"), 
      font.tickslab = c(15, "bold"))
##---Fig2D               
res.cox1 <- coxph(Surv(time_CFS ,event_CFS) ~ Locus_3p_POS + Locus_9p_POS  + Locus_17p_POS +Histology , data =  D) 
summary(res.cox1)$conf.int
summary(res.cox1)$coefficient
}
#################################################################################################### Fig3
if(fig3==T){
##---Fig3A
library(ggstatsplot)
library(ggplot2)
library(ggpubr)
source("/Users/zhaox12/Desktop/Teresas_lab/project/cmd/0_others/Scott/specific_ggstats_plot.R")
set.seed(0)
data<-data.keep
D<-data
D$score<-apply(D[,c('chr3p14.loss','chr17p13.1.loss','chr9p21.3.loss')],1,sum)
D$score<-ifelse(D$score>0,"SCNA","No_SCNA")
D$score <- factor(D$score, levels=c("No_SCNA", "SCNA"))
t1.1<-D[D$score=="SCNA",]
t1.2<-D[D$score=="No_SCNA",]
t1<-wilcox.test(t1.1$CD3D,t1.2$CD3D)
g1<-ggbetweenstats(
  data = D,
  x = score,
  y = CD3D,
  notch = F, # show notched box plot
  mean.plotting = TRUE, # whether mean for each group is to be displayed
  mean.ci = F, # whether to display confidence interval for means
  type = "Parametric", # which type of test is to be run
  k = 2, # number of decimal places for statistical results
  outlier.tagging = F, # whether outliers need to be tagged
  xlab = "", # label for the x-axis variable
  ylab = "expression", # label for the y-axis variable
  title = "CD3", # title text for the plot
  ggtheme = ggplot2::theme_classic(), # choosing a different theme
  ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
  messages = FALSE,
  pairwise.comparisons = F, # display significant pairwise comparisons
  ggplot.component = list(theme(text = element_text(size = 15)))
)+ggplot2::theme(axis.title.y = element_text( size = 25, face = "plain"),
                 axis.text.x =element_text( size = 25, face = "plain"),
                 plot.title = element_text(size=25))+
  ggplot2::annotate("text", x=1.5, y=11.5, label= paste0("p=",format(t1$p.value,scientific=T,digits=3)),size=10)+
  scale_color_manual(values=c("#FF0033","#3300FF"))+scale_x_discrete(expand=c(0.8,0))+coord_cartesian(ylim = c(2.5, 12))
t2<-wilcox.test(t1.1$CD8A,t1.2$CD8A)
g2<-ggbetweenstats(
  data = D,
  x = score,
  y = CD8A,
  notch = F, # show notched box plot
  mean.plotting = TRUE, # whether mean for each group is to be displayed
  mean.ci = F, # whether to display confidence interval for means
  type = "Parametric", # which type of test is to be run
  k = 2, # number of decimal places for statistical results
  outlier.tagging = F, # whether outliers need to be tagged
  xlab = "", # label for the x-axis variable
  ylab = "expression", # label for the y-axis variable
  title = "CD8", # title text for the plot
  ggtheme = ggplot2::theme_classic(), # choosing a different theme
  ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
  messages = FALSE,
  pairwise.comparisons = F, # display significant pairwise comparisons
  ggplot.component = list(theme(text = element_text(size = 15)))
)+ggplot2::theme(axis.title.y = element_text( size = 25, face = "plain"),
                 axis.text.x =element_text( size = 25, face = "plain"),
                 plot.title = element_text(size=25))+
  ggplot2::annotate("text", x=1.5, y=11.5, label= paste0("p=",format(t2$p.value,scientific=T,digits=3)),size=10)+
  scale_color_manual(values=c("#FF0033","#3300FF"))+scale_x_discrete(expand=c(0.8,0))+coord_cartesian(ylim = c(2.5, 12))
figure<-ggarrange(g1,g2,
                  labels = c("",""),
                  ncol = 2, nrow = 1,
                  font.label = list(size = 30, color = "black"))
annotate_figure(figure,
                top = text_grob("", color = "black", face = "bold", size = 28))

##---Fig3B and 3c
##CD3
data<-data.keep
D<-data
D$score<-NA
D$score[D$CD3D>quantile(D$CD3D,0.65,na.rm=T)]<-1
D$score[D$CD3D<quantile(D$CD3D,0.35,na.rm=T)]<-0
D<-D[!is.na(D$score),]
covs.in<-c('chr3p14.loss','chr9p21.3.loss','chr17p13.1.loss','Chrom7.gain','Unweighted.SCNA.0.3.chrom.arm.scale') #could include TP53
covs.in.con<-c('X3p14.2.scale','X9p21.3.scale','X17p13.1.scale','Chrom7.scale','Unweighted.SCNA.0.3.chrom.arm.scale')
fmla<- as.formula(paste("score~", paste(covs.in, collapse= "+")))
fmlb<- as.formula(paste("score~", paste(covs.in.con, collapse= "+")))
des.m1<-model.matrix(fmla, D)
des.m2<-model.matrix(fmlb, D)
mymodel1 <- glm(fmla, data.frame(des.m1,'score'=D$score), family = "binomial")#1 binary
mymodel2<-glm(fmlb, data.frame(des.m2,'score'=D$score), family = "binomial")#2 continous
summary(mymodel1)$coefficients
summary(mymodel2)$coefficients
##CD8
data<-data.keep
D<-data
D$score<-NA
D$score[D$CD8A>quantile(D$CD8A,0.65,na.rm=T)]<-1
D$score[D$CD8A<quantile(D$CD8A,0.35,na.rm=T)]<-0
D<-D[!is.na(D$score),]
covs.in<-c('chr3p14.loss','chr9p21.3.loss','chr17p13.1.loss','Chrom7.gain','Unweighted.SCNA.0.3.chrom.arm.scale') #could include TP53
covs.in.con<-c('X3p14.2.scale','X9p21.3.scale','X17p13.1.scale','Chrom7.scale','Unweighted.SCNA.0.3.chrom.arm.scale')
fmla<- as.formula(paste("score~", paste(covs.in, collapse= "+")))
fmlb<- as.formula(paste("score~", paste(covs.in.con, collapse= "+")))
des.m1<-model.matrix(fmla, D)
des.m2<-model.matrix(fmlb, D)
mymodel1 <- glm(fmla, data.frame(des.m1,'score'=D$score), family = "binomial")#1 binary
mymodel2<-glm(fmlb, data.frame(des.m2,'score'=D$score), family = "binomial")#2 continous
summary(mymodel1)$coefficients
summary(mymodel2)$coefficients
##CD68
data<-data.keep
D<-data
D$score<-NA
D$score[D$CD68>quantile(D$CD68,0.65,na.rm=T)]<-1
D$score[D$CD68<quantile(D$CD68,0.35,na.rm=T)]<-0
D<-D[!is.na(D$score),]
covs.in<-c('chr3p14.loss','chr9p21.3.loss','chr17p13.1.loss','Chrom7.gain','Unweighted.SCNA.0.3.chrom.arm.scale') #could include TP53
covs.in.con<-c('X3p14.2.scale','X9p21.3.scale','X17p13.1.scale','Chrom7.scale','Unweighted.SCNA.0.3.chrom.arm.scale')
fmla<- as.formula(paste("score~", paste(covs.in, collapse= "+")))
fmlb<- as.formula(paste("score~", paste(covs.in.con, collapse= "+")))
des.m1<-model.matrix(fmla, D)
des.m2<-model.matrix(fmlb, D)
mymodel1 <- glm(fmla, data.frame(des.m1,'score'=D$score), family = "binomial")#1 binary
mymodel2<-glm(fmlb, data.frame(des.m2,'score'=D$score), family = "binomial")#2 continous
summary(mymodel1)$coefficients
summary(mymodel2)$coefficients
##---Fig3D
data<-data.keep
D<-data
D$chr9p21.loss.arm_ONLY<-0
D$chr9p21.loss.arm_ONLY[D$X9p_cnv<(-0.3) &D$X9p21.3_cnv_focalonly> (-0.1)]<-1
table(D$chr9p21.loss.arm_ONLY)
D$chr9p21.loss.focal_ONLY<-0
D$chr9p21.loss.focal_ONLY[D$X9p21.3_cnv_focalonly<(-0.3)&D$X9p_cnv>(-0.1)]<-1
table(D$chr9p21.loss.focal_ONLY)
D$chr9p21.loss.arm_AND_FOCAL<-0
D$chr9p21.loss.arm_AND_FOCAL[D$X9p21.3_cnv_focalonly<(-0.1)&D$X9p_cnv<(-0.1)]<-1
table(D$chr9p21.loss.arm_AND_FOCAL)
D$score<-NA
D$score[D$CD8A>quantile(D$CD8A,0.65,na.rm=T)]<-1 #CD8A
D$score[D$CD8A<quantile(D$CD8A,0.35,na.rm=T)]<-0 #CD8A
D<-D[!is.na(D$score),]

D1<-D[,c("X9p_cnv","X9p21.3_cnv_focalonly","chr9p21.loss.arm_ONLY","chr9p21.loss.focal_ONLY","chr9p21.loss.arm_AND_FOCAL","score")]

covs.in<-c('chr9p21.loss.arm_ONLY','chr9p21.loss.focal_ONLY','chr9p21.loss.arm_AND_FOCAL','chr3p14.loss','Unweighted.SCNA.0.3.chrom.arm.scale') #could include TP53
covs.in.con<-c('X9p_cnv.scale','X9p21.3_cnv_focalonly.scale','X3p14.2.scale','Unweighted.SCNA.0.3.chrom.arm.scale')
fmla<- as.formula(paste("score~", paste(covs.in, collapse= "+")))
fmlb<- as.formula(paste("score~", paste(covs.in.con, collapse= "+")))
des.m1<-model.matrix(fmla, D)
des.m2<-model.matrix(fmlb, D)
mymodel1 <- glm(fmla, data.frame(des.m1,'score'=D$score), family = "binomial")#1 binary
mymodel2<-glm(fmlb, data.frame(des.m2,'score'=D$score), family = "binomial")#2 continuous
summary(mymodel1)$coefficients
summary(mymodel2)$coefficients
}
#################################################################################################### Fig4
if(fig4==T){
##---Fig4A-D
library(ggpubr)
library(rstatix)
data.keep<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/0-others/CCLE/Scotte-3.4-others/Teresa_data/final_use/TCGA.HNSC.ABSOLUTE-Feb27.2021-Xin.txt",
                        sep="\t",header = T)
set.seed(1234)
data<-merge(data.keep,teresa.data[,c("sample","TP53_Variant_Classification")],by="sample")
data1<-data[,c("stage.use","CD8A","CD3D","chr9p21.3.loss","Unweighted.SCNA.0.3.chrom.arm.scale","X9p.loss","TP53_Variant_Classification","SASP_epithelial_signature",
               "chr9p21.3_focalonly.loss","CDKN2A.focal.loss")]
data1$chr9p21.3.loss<-ifelse(data1$chr9p21.3.loss==1,"Loss","No_loss")
data1$chr9p21.3_focalonly.loss<-ifelse(data1$chr9p21.3_focalonly.loss==1,"Loss","No_loss")
data1$CDKN2A.focal.loss<-ifelse(data1$CDKN2A.focal.loss==1,"Loss","No_loss")
data1<-data1[data1$stage.use!=0,]
data1<-data1[!is.na(data1$stage.use),]
data1$chr9p21.3.loss<-factor(data1$chr9p21.3.loss,levels = c("No_loss","Loss"))
data1$chr9p21.3_focalonly.loss<-factor(data1$chr9p21.3_focalonly.loss,levels = c("No_loss","Loss"))
data.rep<-data1
data.rep<-data.rep[!is.na(data1$X9p.loss),]
data.rep$Arm9p.loss<-ifelse(data.rep$X9p.loss==1,"Loss","No_loss")
data.rep$Arm9p.loss<-factor(data.rep$Arm9p.loss,levels = c("No_loss","Loss"))

data1$stage.use <- as.factor(data1$stage.use)
#Fig4B
stat.test1 <- data1 %>%
  group_by(stage.use) %>%
  wilcox_test(CD8A ~ chr9p21.3.loss) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test1 <- stat.test1 %>%
  add_xy_position(x = "stage.use", dodge = 0.8)
stat.test1$y.position<-12

n_fun <- function(x){
  return(data.frame(y = 0.95*15.2,
                    label = paste0("N=",length(x))))
}
g1<-ggplot(data1, aes(x=stage.use, y=CD8A)) + 
  geom_point(aes(fill = chr9p21.3.loss, color=chr9p21.3.loss), alpha=0.5,
             size = 2, shape = 21, position = position_jitterdodge()) +
  geom_boxplot(aes(fill = chr9p21.3.loss),alpha = 0.5,outlier.colour = NA)+
  stat_summary(fun.data = n_fun, geom = "text", 
               aes(group=chr9p21.3.loss),
               hjust = 0.5, position = position_dodge(0.8),size=5) +
  scale_color_manual(values=c("#FF0033","#3300FF"))+
  scale_fill_manual(values=c("#FF0033","#3300FF")) +
  theme_classic2()+
  ggtitle("9p21.3 loss (arm or focal) and CD8 by Stage (TCGA)")+
  ylab("CD8")+xlab("")+scale_y_continuous(limits = c(-1, 16))+
  theme(axis.title.y = element_text( size = 15, face = "plain"),
        axis.text.x =element_text( size = 15, face = "plain"),
        plot.title = element_text(size=15))+
  theme(legend.position = c(0.09, 0.2))+
  stat_pvalue_manual(
    stat.test1,  label = "p.adj",tip.length = 0,size=5
  )
#Fig4A
stat.test3 <- data.rep %>%
  group_by(stage.use) %>%
  wilcox_test(CD8A ~ Arm9p.loss) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test3 <- stat.test3 %>%
  add_xy_position(x = "stage.use", dodge = 0.8)
stat.test3$y.position<-12
g3<-ggplot(data.rep, aes(x=stage.use, y=CD8A)) + 
  geom_point(aes(fill = Arm9p.loss, color=Arm9p.loss), alpha=0.5,
             size = 2, shape = 21, position = position_jitterdodge()) +
  geom_boxplot(aes(fill=Arm9p.loss),alpha = 0.5,outlier.colour = NA)+
  stat_summary(fun.data = n_fun, geom = "text", 
               aes(group=Arm9p.loss),
               hjust = 0.5, position = position_dodge(0.8),size=5) +
  scale_color_manual(values=c("#FF0033","#3300FF"))+
  scale_fill_manual(values=c("#FF0033","#3300FF")) +
  theme_classic2()+
  ggtitle("9p arm loss and CD8 by Stage (TCGA)")+
  ylab("CD8")+xlab("")+scale_y_continuous(limits = c(-1, 16))+
  theme(axis.title.y = element_text( size = 15, face = "plain"),
        axis.text.x =element_text( size = 15, face = "plain"),
        plot.title = element_text(size=15))+
  theme(legend.position = c(0.08, 0.2))+
  stat_pvalue_manual(
    stat.test3,  label = "p.adj",tip.length = 0,size=5
  )
group<-c(0,1)
used.stage<-c("Stage I","Stage II")

data2<-data.rep
data2[, 7][is.na(data2[, 7])] <- 0
data2<-data2[data2$TP53_Variant_Classification %in% group,]
data2<-data2[data2$stage.use %in% used.stage,]
table(data2$stage.use,exclude = NULL)
data2$TP53_mut<-ifelse(data2$TP53_Variant_Classification==1,"Mut","WT")
data2$TP53_mut<-factor(data2$TP53_mut,levels = c("WT","Mut"))
data2$SCNA.level<-ifelse(data2$Unweighted.SCNA.0.3.chrom.arm.scale>quantile(data2$Unweighted.SCNA.0.3.chrom.arm.scale,0.5,na.rm=T),"High","Low")
data2$SCNA.level <- factor(data2$SCNA.level, levels=c("Low", "High"))

stat.test4 <- data2 %>%
  group_by(TP53_mut) %>%
  wilcox_test(CD8A ~ Arm9p.loss) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test4 <- stat.test4 %>%
  add_xy_position(x = "TP53_mut", dodge = 0.8)
stat.test4$y.position<-12
#fig 4C
g4<-ggplot(data2, aes(x=TP53_mut, y=CD8A)) + 
  geom_point(aes(fill = Arm9p.loss, color=Arm9p.loss), alpha=0.5,
             size = 2, shape = 21, position = position_jitterdodge()) +
  geom_boxplot(aes(fill = Arm9p.loss),alpha = 0.5,outlier.colour = NA)+
  stat_summary(fun.data = n_fun, geom = "text", 
               aes(group=Arm9p.loss),
               hjust = 0.5, position = position_dodge(0.8),size=5) +
  scale_color_manual(values=c("#FF0033","#3300FF"))+
  scale_fill_manual(values=c("#FF0033","#3300FF")) +
  theme_classic2()+
  ggtitle("9p arm and CD8 by TP53")+
  ylab("CD8")+xlab("")+scale_y_continuous(limits = c(-1, 16))+
  theme(axis.title.y = element_text( size = 15, face = "plain"),
        axis.text.x =element_text( size = 15, face = "plain"),
        plot.title = element_text(size=15))+
  theme(legend.position = c(0.168, 0.2))+
  stat_pvalue_manual(
    stat.test4,  label = "p",tip.length = 0,size=5
  )

data2<-data1
data2[, 7][is.na(data2[, 7])] <- 0
data2<-data2[data2$TP53_Variant_Classification %in% group,]
data2<-data2[data2$stage.use %in% used.stage,]
table(data2$stage.use,exclude = NULL)
data2$TP53_mut<-ifelse(data2$TP53_Variant_Classification==1,"Mut","WT")
data2$TP53_mut<-factor(data2$TP53_mut,levels = c("WT","Mut"))
data2$SCNA.level<-ifelse(data2$Unweighted.SCNA.0.3.chrom.arm.scale>quantile(data2$Unweighted.SCNA.0.3.chrom.arm.scale,0.5,na.rm=T),"High","Low")
data2$SCNA.level <- factor(data2$SCNA.level, levels=c("Low", "High"))
stat.test6 <- data2 %>%
  group_by(TP53_mut) %>%
  wilcox_test(CD8A ~ chr9p21.3.loss) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test6 <- stat.test6 %>%
  add_xy_position(x = "TP53_mut", dodge = 0.8)
stat.test6$y.position<-12
#Fig4D
g6<-ggplot(data2, aes(x=TP53_mut, y=CD8A)) + 
  geom_point(aes(fill = chr9p21.3.loss, color=chr9p21.3.loss), alpha=0.5,
             size = 2, shape = 21, position = position_jitterdodge()) +
  geom_boxplot(aes(fill = chr9p21.3.loss),alpha = 0.5,outlier.colour = NA)+
  scale_color_manual(values=c("#FF0033","#3300FF"))+
  stat_summary(fun.data = n_fun, geom = "text", 
               aes(group=chr9p21.3.loss),
               hjust = 0.5, position = position_dodge(0.8),size=5) +
  scale_fill_manual(values=c("#FF0033","#3300FF")) +
  theme_classic2()+
  ggtitle("9p21.3 (arm or focal) and CD8 by TP53")+
  ylab("CD8")+xlab("")+scale_y_continuous(limits = c(-1, 16))+
  theme(axis.title.y = element_text( size = 15, face = "plain"),
        axis.text.x =element_text( size = 15, face = "plain"),
        plot.title = element_text(size=15))+
  theme(legend.position = c(0.168, 0.2))+
  stat_pvalue_manual(
    stat.test6,  label = "p",tip.length = 0,size=5
  )
figure<-ggarrange(
  g3,g1,
  ggarrange(g4,g6, ncol = 2, labels = "", align = "h",widths = c(1,1)),
  nrow = 3, font.label = list(size = 30, color = "black"),heights = c(1,1,1))
}
#################################################################################################### Fig5
if(fig5==T){
##---Fig5A from GSEA attachement
##---Fig5B~E
# CCLE
library(ggpubr)
library(rstatix)
data<-read.delim('/Users/zhaox12/Desktop/Teresas_lab/project/0-others/CCLE/Scotte-3.4-others/Teresa_data/final_use/CCLE.HNSC.alldata-Feb27.2021-Xin.txt',as.is=T)
datakeep<-data
###Fig5B~E
  D<-datakeep
  data1<-D[,c("CD274","CD8A","KEGG_JAK_STAT_signature","KEGG_Cytokine_Cytokine_receptor_interaction_signature","chr9p21.3.loss","chr3p14.loss","X9p_cnv.loss","X3p_cnv.loss",
              "TP53_mut_2019","chr9p24.1.loss","SASP_epithelial_signature","chr9p21.3_focal_only.loss","TNFA_SIGNALING_VIA_NFKB_signature","CDKN2A.loss","cell_cycle_signature")]
  colnames(data1)[7]<-"Chr9p_cnv.loss"
  colnames(data1)[8]<-"Chr3p_cnv.loss"
  data1$chr9p21.3.loss<-ifelse(data1$chr9p21.3.loss==1,"Loss","No_loss")
  data1$chr9p21.3.loss<-factor(data1$chr9p21.3.loss,levels = c("No_loss","Loss"))
  data1$chr3p14.loss<-ifelse(data1$chr3p14.loss==1,"Loss","No_loss")
  data1$chr3p14.loss<-factor(data1$chr3p14.loss,levels = c("No_loss","Loss"))
  data1$Chr9p_cnv.loss<-ifelse(data1$Chr9p_cnv.loss==1,"Loss","No_loss")
  data1$Chr9p_cnv.loss<-factor(data1$Chr9p_cnv.loss,levels = c("No_loss","Loss"))
  data1$Chr3p_cnv.loss<-ifelse(data1$Chr3p_cnv.loss==1,"Loss","No_loss")
  data1$Chr3p_cnv.loss<-factor(data1$Chr3p_cnv.loss,levels = c("No_loss","Loss"))
  data1$chr9p24.1.loss<-ifelse(data1$chr9p24.1.loss==1,"Loss","No_loss")
  data1$chr9p24.1.loss<-factor(data1$chr9p24.1.loss,levels = c("No_loss","Loss"))
  data1$CDKN2A.loss<-ifelse(data1$CDKN2A.loss==1,"Loss","No_loss")
  data1$CDKN2A.loss<-factor(data1$chr9p24.1.loss,levels = c("No_loss","Loss"))
  data.rep<-data1
  data.rep<-data.rep[!is.na(data.rep$chr9p21.3_focal_only.loss),]
  data.rep$chr9p21.3.focal.loss<-ifelse(data.rep$chr9p21.3_focal_only.loss==1,"Loss","No_loss")
  data.rep$chr9p21.3.focal.loss<-factor(data.rep$chr9p21.3.focal.loss,levels = c("No_loss","Loss"))
  
  
  stat.test1 <- data1 %>%
    wilcox_test(SASP_epithelial_signature ~ chr9p21.3.loss,alternative = "greater") %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  stat.test1$y.position<-0.4
  stat.test2 <- data1 %>%
    wilcox_test(SASP_epithelial_signature ~ Chr9p_cnv.loss,alternative = "greater") %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  stat.test2$y.position<-0.4
  
  n_fun <- function(x){
    return(data.frame(y = 0.95*0.55,
                      label = paste0("N=",length(x))))
  }
  g1<-ggplot(data1, aes(x=chr9p21.3.loss, y=SASP_epithelial_signature)) + 
    geom_point(aes(fill = chr9p21.3.loss, color=chr9p21.3.loss), alpha=0.5,
               size = 2, shape = 21, position = position_jitterdodge()) +
    geom_boxplot(aes(fill = chr9p21.3.loss),alpha = 0.5)+
    stat_summary(fun.data = n_fun, geom = "text", 
                 aes(group=chr9p21.3.loss),
                 hjust = 0.5, position = position_dodge(0.8),size=7) +
    scale_fill_manual(values=c("#FF0033","#3300FF")) +
    scale_color_manual(values=c("#FF0033","#3300FF"))+
    theme_classic2()+
    ggtitle("9p21.3 loss (arm or focal) and SASP")+
    ylab("SASP_enrichment")+xlab("")+scale_y_continuous(limits = c(-0.5, 0.55))+
    theme(axis.title.y = element_text( size = 25, face = "plain"),
          axis.text.x =element_text( size = 25, face = "plain"),
          plot.title = element_text(size=20),
          legend.title=element_text(size=16), 
          legend.text=element_text(size=15))+
    theme(legend.position = c(0.18, 0.22))+
    stat_pvalue_manual(
      stat.test1,  label = "p.adj",tip.length = 0,size = 7
    )
  
  g2<-ggplot(data1, aes(x=Chr9p_cnv.loss, y=SASP_epithelial_signature)) + 
    geom_point(aes(fill = Chr9p_cnv.loss, color=Chr9p_cnv.loss), alpha=0.5,
               size = 2, shape = 21, position = position_jitterdodge()) +
    geom_boxplot(aes(fill = Chr9p_cnv.loss),alpha = 0.5)+
    stat_summary(fun.data = n_fun, geom = "text", 
                 aes(group=Chr9p_cnv.loss),
                 hjust = 0.5, position = position_dodge(0.8),size=7) +
    scale_fill_manual(values=c("#FF0033","#3300FF")) +
    scale_color_manual(values=c("#FF0033","#3300FF"))+
    theme_classic2()+
    ggtitle("9p arm loss and SASP")+
    ylab("SASP_enrichment")+xlab("")+scale_y_continuous(limits = c(-0.5, 0.55))+
    theme(text = element_text(size = 15),
          axis.title.y = element_text( size = 25, face = "plain"),
          axis.text.x =element_text( size = 25, face = "plain"),
          plot.title = element_text(size=20),
          legend.title=element_text(size=16), 
          legend.text=element_text(size=15))+
    theme(legend.position = c(0.18, 0.22))+
    stat_pvalue_manual(
      stat.test2,  label = "p.adj",tip.length = 0,size=7
    )
  n_fun1 <- function(x){
    return(data.frame(y = 0.95*0.2,
                      label = paste0("N=",length(x))))
  }
  stat.test3 <- data1 %>%
    wilcox_test(KEGG_JAK_STAT_signature ~ chr9p21.3.loss,alternative = "greater") %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  stat.test3$y.position<-0.1
  
  g3<-ggplot(data1, aes(x=chr9p21.3.loss, y=KEGG_JAK_STAT_signature)) + 
    geom_point(aes(fill = chr9p21.3.loss, color=chr9p21.3.loss), alpha=0.5,
               size = 2, shape = 21, position = position_jitterdodge()) +
    geom_boxplot(aes(fill = chr9p21.3.loss),alpha = 0.5)+
    stat_summary(fun.data = n_fun1, geom = "text", 
                 aes(group=chr9p21.3.loss),
                 hjust = 0.5, position = position_dodge(0.8),size=7) +
    scale_fill_manual(values=c("#FF0033","#3300FF")) +
    scale_color_manual(values=c("#FF0033","#3300FF"))+
    theme_classic2()+
    ggtitle("9p21.3 loss (arm or focal) and JAK-STAT")+
    ylab("JAK-STAT_enrichment")+xlab("")+scale_y_continuous(limits = c(-0.3, 0.2))+
    theme(axis.title.y = element_text( size = 25, face = "plain"),
          axis.text.x =element_text(size = 25, face = "plain"),
          plot.title = element_text(size=20),
          legend.title=element_text(size=16), 
          legend.text=element_text(size=15))+
    theme(legend.position = c(0.18, 0.22))+
    ylim(-0.4,0.2)+
    stat_pvalue_manual(
      stat.test3,  label = "p.adj",tip.length = 0,size = 7
    )
  
  stat.test4 <- data1 %>%
    wilcox_test(KEGG_JAK_STAT_signature ~ Chr9p_cnv.loss,alternative = "greater") %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  stat.test4$y.position<-0.1
  g4<-ggplot(data1, aes(x=Chr9p_cnv.loss, y=KEGG_JAK_STAT_signature)) + 
    geom_point(aes(fill = Chr9p_cnv.loss, color=Chr9p_cnv.loss), alpha=0.5,
               size = 2, shape = 21, position = position_jitterdodge()) +
    geom_boxplot(aes(fill = Chr9p_cnv.loss),alpha = 0.5)+
    stat_summary(fun.data = n_fun1, geom = "text", 
                 aes(group=Chr9p_cnv.loss),
                 hjust = 0.5, position = position_dodge(0.8),size=7) +
    scale_fill_manual(values=c("#FF0033","#3300FF")) +
    scale_color_manual(values=c("#FF0033","#3300FF"))+
    theme_classic2()+
    ggtitle("9p arm loss and JAK-STAT")+
    ylab("JAK-STAT_enrichment")+xlab("")+scale_y_continuous(limits = c(-0.4, 0.2))+
    theme(axis.title.y = element_text( size = 25, face = "plain"),
          axis.text.x =element_text(size = 25, face = "plain"),
          plot.title = element_text(size=20),
          legend.title=element_text(size=16), 
          legend.text=element_text(size=15))+
    theme(legend.position = c(0.18, 0.22))+
    ylim(-0.4,0.2)+
    stat_pvalue_manual(
      stat.test4,  label = "p.adj",tip.length = 0,size = 7
    )
  stat.test5 <- data1 %>%
    wilcox_test(KEGG_Cytokine_Cytokine_receptor_interaction_signature ~ chr9p21.3.loss,alternative = "greater") %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  stat.test5$y.position<-0.1
  g5<-ggplot(data1, aes(x=chr9p21.3.loss, y=KEGG_Cytokine_Cytokine_receptor_interaction_signature)) + 
    geom_point(aes(fill = chr9p21.3.loss, color=chr9p21.3.loss), alpha=0.5,
               size = 2, shape = 21, position = position_jitterdodge()) +
    geom_boxplot(aes(fill = chr9p21.3.loss),alpha = 0.5)+
    stat_summary(fun.data = n_fun1, geom = "text", 
                 aes(group=chr9p21.3.loss),
                 hjust = 0.5, position = position_dodge(0.8),size=7) +
    scale_fill_manual(values=c("#FF0033","#3300FF")) +
    scale_color_manual(values=c("#FF0033","#3300FF"))+
    theme_classic2()+
    ggtitle("9p21.3 loss (arm or focal) and Cytokine")+
    ylab("Cytokine_enrichment")+xlab("")+scale_y_continuous(limits = c(-0.4, 0.2))+
    theme(axis.title.y = element_text( size = 25, face = "plain"),
          axis.text.x =element_text(size = 25, face = "plain"),
          plot.title = element_text(size=20),
          legend.title=element_text(size=16), 
          legend.text=element_text(size=15))+
    theme(legend.position = c(0.18, 0.22))+
    stat_pvalue_manual(
      stat.test5,  label = "p.adj",tip.length = 0,size = 7
    )
  
  stat.test6 <- data1 %>%
    wilcox_test(KEGG_Cytokine_Cytokine_receptor_interaction_signature ~ Chr9p_cnv.loss,alternative = "greater") %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  stat.test6$y.position<-0.1
  g6<-ggplot(data1, aes(x=Chr9p_cnv.loss, y=KEGG_Cytokine_Cytokine_receptor_interaction_signature)) + 
    geom_point(aes(fill = Chr9p_cnv.loss, color=Chr9p_cnv.loss), alpha=0.5,
               size = 2, shape = 21, position = position_jitterdodge()) +
    geom_boxplot(aes(fill = Chr9p_cnv.loss),alpha = 0.5)+
    stat_summary(fun.data = n_fun1, geom = "text", 
                 aes(group=Chr9p_cnv.loss),
                 hjust = 0.5, position = position_dodge(0.8),size=7) +
    scale_fill_manual(values=c("#FF0033","#3300FF")) +
    scale_color_manual(values=c("#FF0033","#3300FF"))+
    theme_classic2()+
    ggtitle("9p arm loss and Cytokine")+
    ylab("Cytokine_enrichment")+xlab("")+scale_y_continuous(limits = c(-0.4, 0.2))+
    theme(axis.title.y = element_text( size = 25, face = "plain"),
          axis.text.x =element_text(size = 25, face = "plain"),
          plot.title = element_text(size=20),
          legend.title=element_text(size=16), 
          legend.text=element_text(size=15))+
    theme(legend.position = c(0.18, 0.22))+
    stat_pvalue_manual(
      stat.test6,  label = "p.adj",tip.length = 0,size = 7
    )
  stat.test7 <- data1 %>%
    wilcox_test(TNFA_SIGNALING_VIA_NFKB_signature ~ chr9p21.3.loss,alternative = "greater") %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  stat.test7$y.position<-0.45
  g7<-ggplot(data1, aes(x=chr9p21.3.loss, y=TNFA_SIGNALING_VIA_NFKB_signature)) + 
    geom_point(aes(fill = chr9p21.3.loss, color=chr9p21.3.loss), alpha=0.5,
               size = 2, shape = 21, position = position_jitterdodge()) +
    geom_boxplot(aes(fill = chr9p21.3.loss),alpha = 0.5)+
    stat_summary(fun.data = n_fun, geom = "text", 
                 aes(group=chr9p21.3.loss),
                 hjust = 0.5, position = position_dodge(0.8),size=7) +
    scale_fill_manual(values=c("#FF0033","#3300FF")) +
    scale_color_manual(values=c("#FF0033","#3300FF"))+
    theme_classic2()+
    ggtitle("9p21.3 loss (arm or focal) and TNFA")+
    ylab("TNFA_enrichment")+xlab("")+scale_y_continuous(limits = c(-0.1, 0.55))+
    theme(axis.title.y = element_text( size = 25, face = "plain"),
          axis.text.x =element_text(size = 25, face = "plain"),
          plot.title = element_text(size=20),
          legend.title=element_text(size=16), 
          legend.text=element_text(size=15))+
    theme(legend.position = c(0.18, 0.22))+
    stat_pvalue_manual(
      stat.test7,  label = "p.adj",tip.length = 0,size = 7
    )
  
  stat.test8 <- data1 %>%
    wilcox_test(TNFA_SIGNALING_VIA_NFKB_signature ~ Chr9p_cnv.loss,alternative = "greater") %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  stat.test8$y.position<-0.45
  g8<-ggplot(data1, aes(x=Chr9p_cnv.loss, y=TNFA_SIGNALING_VIA_NFKB_signature)) + 
    geom_point(aes(fill = Chr9p_cnv.loss, color=Chr9p_cnv.loss), alpha=0.5,
               size = 2, shape = 21, position = position_jitterdodge()) +
    geom_boxplot(aes(fill = Chr9p_cnv.loss),alpha = 0.5)+
    stat_summary(fun.data = n_fun, geom = "text", 
                 aes(group=Chr9p_cnv.loss),
                 hjust = 0.5, position = position_dodge(0.8),size=7) +
    scale_fill_manual(values=c("#FF0033","#3300FF")) +
    scale_color_manual(values=c("#FF0033","#3300FF"))+
    theme_classic2()+
    ggtitle("9p arm loss and TNFA")+
    ylab("TNFA_enrichment")+xlab("")+scale_y_continuous(limits = c(-0.1, 0.55))+
    theme(axis.title.y = element_text( size = 25, face = "plain"),
          axis.text.x =element_text(size = 25, face = "plain"),
          plot.title = element_text(size=20),
          legend.title=element_text(size=16), 
          legend.text=element_text(size=15))+
    theme(legend.position = c(0.18, 0.22))+
    stat_pvalue_manual(
      stat.test8,  label = "p.adj",tip.length = 0,size = 7
    )
figure<-ggarrange(
  g2,g4,g6,g8,g1,g3,g5,g7,
  nrow = 2,ncol=4, font.label = list(size = 30, color = "black"))
}
#################################################################################################### Fig6
##---Fig6 all data/plots from Caris
#################################################################################################### Fig7
##---Fig7 from Catherine Eng
#################################################################################################### FigS1
if(figS1==T){
par(mfrow=c(3,1))
plot(data.keep$Immune.Signature,data.keep$CD8A, pch = 19, col = "lightblue",xlab="Immune Score", ylab="CD8A")
abline(lm(data.keep$CD8A~data.keep$Immune.Signature), col = "red", lwd = 3)
text(paste("Correlation:", round(cor(data.keep$Immune.Signature,data.keep$CD8A), 2)), x = 55, y = 10.5)

plot(data.keep$Immune.Signature,data.keep$CD3D, pch = 19, col = "lightblue",xlab="Immune Score", ylab="CD3D")
abline(lm(data.keep$CD3D~data.keep$Immune.Signature), col = "red", lwd = 3)
text(paste("Correlation:", round(cor(data.keep$Immune.Signature,data.keep$CD3D), 2)), x = 55, y = 10)

plot(data.keep$Immune.Signature,data.keep$CD68, pch = 19, col = "lightblue",xlab="Immune Score", ylab="CD68")
abline(lm(data.keep$CD68~data.keep$Immune.Signature), col = "red", lwd = 3)
text(paste("Correlation:", round(cor(data.keep$Immune.Signature,data.keep$CD68), 2)), x = 55, y = 13.5)
}
#################################################################################################### FigS2
##---FigS2 is based on Teresa's scripts
#################################################################################################### FigS3
if(figS3==T){
###---FigS3
library(ggstatsplot)
library(ggplot2)
library(ggpubr)
source("/Users/zhaox12/Desktop/Teresas_lab/project/cmd/0_others/Scott/specific_ggstats_plot.R")

set.seed(1234)
cn_data<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/0-others/CCLE/Scotte-3.4-others/TCGA/hnsc/broad_values_by_arm.txt",
                    sep="\t",header = T)
rownames(cn_data)<-make.names(cn_data$Chromosome.Arm)
cn_data1<-as.data.frame(t(cn_data[,-1]))
rownames(cn_data1)<-substr(rownames(cn_data1),1,12)
cn_Loss<-cn_data1
cn_Loss[-0.3<cn_Loss]<-0
cn_Loss[-0.3>cn_Loss]<-1
cn_Loss$Loss.sum<-apply(cn_Loss,1,sum)
cn.Loss1<-cn_Loss[,c(1,42)]

cn_amp<-cn_data1
cn_amp[0.3<cn_amp]<-1
cn_amp[0.3>cn_amp]<-0
cn_amp$amp.sum<-apply(cn_amp,1,sum)
cn.amp1<-cn_amp[,c(1,42)]

cn_cnv<-cn_data1
cn_cnv[0.3<cn_cnv]<-1
cn_cnv[-0.3>cn_cnv]<-1
cn_cnv[-0.3<cn_cnv & cn_cnv<0.3]<-0
cn_cnv$cnv.sum<-apply(cn_cnv,1,sum)
cn.cnv1<-cn_cnv[,c(1,42)]

data1<-cbind(cn.Loss1,cn.amp1,cn.cnv1)
data2<-data1[,c(-1,-3,-5)]
data2$sample<-rownames(data2)

data3<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/0-others/CCLE/Scotte-3.4-others/Teresa_data/final_use/3.8-tcga.summary_table_TRANS.txt",
                  sep="\t",header = T)
merge.data<-merge(data3,data2,by="sample")
merge.data$region_3p14 <- factor(merge.data$region_3p14, levels=c("No_loss", "Loss"))
merge.data$region_9p21 <- factor(merge.data$region_9p21, levels=c("No_loss", "Loss"))
merge.data$region_17p13.1 <- factor(merge.data$region_17p13.1, levels=c("No_loss", "Loss"))
merge.data$three_region_Loss <- factor(merge.data$three_region_loss, levels=c("No_loss", "Loss"))

t1.1<-merge.data[merge.data$region_3p14=="Loss",]
t1.2<-merge.data[merge.data$region_3p14=="No_loss",]
t1<-wilcox.test(t1.1$Loss.sum,t1.2$Loss.sum)
g1<-ggbetweenstats(
  data = merge.data,
  x = region_3p14,
  y = Loss.sum,
  notch = F, # show notched box plot
  mean.plotting = TRUE, # whether mean for each group is to be displayed
  mean.ci = F, # whether to display confidence interval for means
  type = "Parametric", # which type of test is to be run
  k = 2, # number of decimal places for statistical results
  outlier.tagging = F, # whether outliers need to be tagged
  xlab = "", # label for the x-axis variable
  ylab = "total N of arm losses", # label for the y-axis variable
  title = "3p14", # title text for the plot
  ggtheme = ggplot2::theme_classic(), # choosing a different theme
  ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
  messages = FALSE,
  pairwise.comparisons = F, # display significant pairwise comparisons
  ggplot.component = list(theme(text = element_text(size = 15)))
)+ggplot2::theme(axis.title.y = element_text( size = 30, face = "plain"),
                 axis.text.x =element_text( size = 30, face = "plain"),
                 plot.title = element_text(size=30))+
  ggplot2::annotate("text", x=1.5, y=0.95*max(merge.data$Loss.sum), label= paste0("p=",signif(t1$p.value,3)),size=10)+
  scale_color_manual(values=c("#FF0033","#3300FF"))

t2.1<-merge.data[merge.data$region_3p14=="Loss",]
t2.2<-merge.data[merge.data$region_3p14=="No_loss",]
t2<-wilcox.test(t2.1$amp.sum,t2.2$amp.sum)
g2<-ggbetweenstats(
  data = merge.data,
  x = region_3p14,
  y = amp.sum,
  notch = F, # show notched box plot
  mean.plotting = TRUE, # whether mean for each group is to be displayed
  mean.ci = F, # whether to display confidence interval for means
  type = "Parametric", # which type of test is to be run
  k = 2, # number of decimal places for statistical results
  outlier.tagging = F, # whether outliers need to be tagged
  xlab = "", # label for the x-axis variable
  ylab = "total N of arm gains", # label for the y-axis variable
  title = "", # title text for the plot
  ggtheme = ggplot2::theme_classic(), # choosing a different theme
  ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
  messages = FALSE,
  pairwise.comparisons = F, # display significant pairwise comparisons
  ggplot.component = list(theme(text = element_text(size = 15)))
)+ggplot2::theme(axis.title.y = element_text( size = 30, face = "plain"),
                 axis.text.x =element_text( size = 30, face = "plain"),
                 plot.title = element_text(size=30))+
  ggplot2::annotate("text", x=1.5, y=0.95*max(merge.data$amp.sum), label= paste0("p=",signif(t2$p.value,3)),size=10)+
  scale_color_manual(values=c("#FF0033","#3300FF"))

t3.1<-merge.data[merge.data$region_3p14=="Loss",]
t3.2<-merge.data[merge.data$region_3p14=="No_loss",]
t3<-wilcox.test(t3.1$cnv.sum,t3.2$cnv.sum)
g3<-ggbetweenstats(
  data = merge.data,
  x = region_3p14,
  y = cnv.sum,
  notch = F, # show notched box plot
  mean.plotting = TRUE, # whether mean for each group is to be displayed
  mean.ci = F, # whether to display confidence interval for means
  type = "Parametric", # which type of test is to be run
  k = 2, # number of decimal places for statistical results
  outlier.tagging = F, # whether outliers need to be tagged
  xlab = "", # label for the x-axis variable
  ylab = "total N of arm losses or gains", # label for the y-axis variable
  title = "", # title text for the plot
  ggtheme = ggplot2::theme_classic(), # choosing a different theme
  ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
  messages = FALSE,
  pairwise.comparisons = F, # display significant pairwise comparisons
  ggplot.component = list(theme(text = element_text(size = 15)))
)+ggplot2::theme(axis.title.y = element_text( size = 30, face = "plain"),
                 axis.text.x =element_text( size = 30, face = "plain"),
                 plot.title = element_text(size=30))+
  ggplot2::annotate("text", x=1.5, y=0.95*max(merge.data$cnv.sum), label= paste0("p=",signif(t3$p.value,3)),size=10)+
  scale_color_manual(values=c("#FF0033","#3300FF"))

t4.1<-merge.data[merge.data$region_9p21=="Loss",]
t4.2<-merge.data[merge.data$region_9p21=="No_loss",]
t4<-wilcox.test(t4.1$Loss.sum,t4.2$Loss.sum)
g4<-ggbetweenstats(
  data = merge.data,
  x = region_9p21,
  y = Loss.sum,
  notch = F, # show notched box plot
  mean.plotting = TRUE, # whether mean for each group is to be displayed
  mean.ci = F, # whether to display confidence interval for means
  type = "Parametric", # which type of test is to be run
  k = 2, # number of decimal places for statistical results
  outlier.tagging = F, # whether outliers need to be tagged
  xlab = "", # label for the x-axis variable
  ylab = "total N of arm losses", # label for the y-axis variable
  title = "9p21.3", # title text for the plot
  ggtheme = ggplot2::theme_classic(), # choosing a different theme
  ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
  messages = FALSE,
  pairwise.comparisons = F, # display significant pairwise comparisons
  ggplot.component = list(theme(text = element_text(size = 15)))
)+ggplot2::theme(axis.title.y = element_text( size = 30, face = "plain"),
                 axis.text.x =element_text( size = 30, face = "plain"),
                 plot.title = element_text(size=30))+
  ggplot2::annotate("text", x=1.5, y=0.95*max(merge.data$Loss.sum), label= paste0("p=",signif(t4$p.value,3)),size=10)+
  scale_color_manual(values=c("#FF0033","#3300FF"))

t5.1<-merge.data[merge.data$region_9p21=="Loss",]
t5.2<-merge.data[merge.data$region_9p21=="No_loss",]
t5<-wilcox.test(t5.1$amp.sum,t5.2$amp.sum)
g5<-ggbetweenstats(
  data = merge.data,
  x = region_9p21,
  y = amp.sum,
  notch = F, # show notched box plot
  mean.plotting = TRUE, # whether mean for each group is to be displayed
  mean.ci = F, # whether to display confidence interval for means
  type = "Parametric", # which type of test is to be run
  k = 2, # number of decimal places for statistical results
  outlier.tagging = F, # whether outliers need to be tagged
  xlab = "", # label for the x-axis variable
  ylab = "total N of arm gains", # label for the y-axis variable
  title = "", # title text for the plot
  ggtheme = ggplot2::theme_classic(), # choosing a different theme
  ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
  package = "ggsci", # package from which color palette is to be taken
  palette = "nrc_npg", # choosing a different color palette
  messages = FALSE,
  pairwise.comparisons = F, # display significant pairwise comparisons
  ggplot.component = list(theme(text = element_text(size = 15)))
)+ggplot2::theme(axis.title.y = element_text( size = 30, face = "plain"),
                 axis.text.x =element_text( size = 30, face = "plain"),
                 plot.title = element_text(size=30))+
  ggplot2::annotate("text", x=1.5, y=0.95*max(merge.data$amp.sum), label= paste0("p=",signif(t5$p.value,3)),size=10)+
  scale_color_manual(values=c("#FF0033","#3300FF"))

t6.1<-merge.data[merge.data$region_9p21=="Loss",]
t6.2<-merge.data[merge.data$region_9p21=="No_loss",]
t6<-wilcox.test(t6.1$cnv.sum,t6.2$cnv.sum)
g6<-ggbetweenstats(
  data = merge.data,
  x = region_9p21,
  y = cnv.sum,
  notch = F, # show notched box plot
  mean.plotting = TRUE, # whether mean for each group is to be displayed
  mean.ci = F, # whether to display confidence interval for means
  type = "Parametric", # which type of test is to be run
  k = 2, # number of decimal places for statistical results
  outlier.tagging = F, # whether outliers need to be tagged
  xlab = "", # label for the x-axis variable
  ylab = "total N of arm losses or gains", # label for the y-axis variable
  title = "", # title text for the plot
  ggtheme = ggplot2::theme_classic(), # choosing a different theme
  ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
  messages = FALSE,
  pairwise.comparisons = F, # display significant pairwise comparisons
  ggplot.component = list(theme(text = element_text(size = 15)))
)+ggplot2::theme(axis.title.y = element_text( size = 30, face = "plain"),
                 axis.text.x =element_text( size = 30, face = "plain"),
                 plot.title = element_text(size=30))+
  ggplot2::annotate("text", x=1.5, y=0.95*max(merge.data$cnv.sum), label= paste0("p=",signif(t6$p.value,3)),size=10)+
  scale_color_manual(values=c("#FF0033","#3300FF"))

t7.1<-merge.data[merge.data$region_17p13.1=="Loss",]
t7.2<-merge.data[merge.data$region_17p13.1=="No_loss",]
t7<-wilcox.test(t7.1$Loss.sum,t7.2$Loss.sum)
g7<-ggbetweenstats(
  data = merge.data,
  x = region_17p13.1,
  y = Loss.sum,
  notch = F, # show notched box plot
  mean.plotting = TRUE, # whether mean for each group is to be displayed
  mean.ci = F, # whether to display confidence interval for means
  type = "Parametric", # which type of test is to be run
  k = 2, # number of decimal places for statistical results
  outlier.tagging = F, # whether outliers need to be tagged
  xlab = "", # label for the x-axis variable
  ylab = "total N of arm losses", # label for the y-axis variable
  title = "17p13.1", # title text for the plot
  ggtheme = ggplot2::theme_classic(), # choosing a different theme
  ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
  messages = FALSE,
  pairwise.comparisons = F, # display significant pairwise comparisons
  ggplot.component = list(theme(text = element_text(size = 15)))
)+ggplot2::theme(axis.title.y = element_text( size = 30, face = "plain"),
                 axis.text.x =element_text( size = 30, face = "plain"),
                 plot.title = element_text(size=30))+
  ggplot2::annotate("text", x=1.5, y=0.95*max(merge.data$Loss.sum), label= paste0("p=",signif(t7$p.value,3)),size=10)+
  scale_color_manual(values=c("#FF0033","#3300FF"))

t8.1<-merge.data[merge.data$region_17p13.1=="Loss",]
t8.2<-merge.data[merge.data$region_17p13.1=="No_loss",]
t8<-wilcox.test(t8.1$amp.sum,t8.2$amp.sum)
g8<-ggbetweenstats(
  data = merge.data,
  x = region_17p13.1,
  y = amp.sum,
  notch = F, # show notched box plot
  mean.plotting = TRUE, # whether mean for each group is to be displayed
  mean.ci = F, # whether to display confidence interval for means
  type = "Parametric", # which type of test is to be run
  k = 2, # number of decimal places for statistical results
  outlier.tagging = F, # whether outliers need to be tagged
  xlab = "", # label for the x-axis variable
  ylab = "total N of arm gains", # label for the y-axis variable
  title = "", # title text for the plot
  ggtheme = ggplot2::theme_classic(), # choosing a different theme
  ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
  messages = FALSE,
  pairwise.comparisons = F, # display significant pairwise comparisons
  ggplot.component = list(theme(text = element_text(size = 15)))
)+ggplot2::theme(axis.title.y = element_text( size = 30, face = "plain"),
                 axis.text.x =element_text( size = 30, face = "plain"),
                 plot.title = element_text(size=30))+
  ggplot2::annotate("text", x=1.5, y=0.95*max(merge.data$amp.sum), label= paste0("p=",signif(t8$p.value,3)),size=10)+
  scale_color_manual(values=c("#FF0033","#3300FF"))

t9.1<-merge.data[merge.data$region_17p13.1=="Loss",]
t9.2<-merge.data[merge.data$region_17p13.1=="No_loss",]
t9<-wilcox.test(t9.1$cnv.sum,t9.2$cnv.sum)
g9<-ggbetweenstats(
  data = merge.data,
  x = region_17p13.1,
  y = cnv.sum,
  notch = F, # show notched box plot
  mean.plotting = TRUE, # whether mean for each group is to be displayed
  mean.ci = F, # whether to display confidence interval for means
  type = "Parametric", # which type of test is to be run
  k = 2, # number of decimal places for statistical results
  outlier.tagging = F, # whether outliers need to be tagged
  xlab = "", # label for the x-axis variable
  ylab = "total N of arm losses or gains", # label for the y-axis variable
  title = "", # title text for the plot
  ggtheme = ggplot2::theme_classic(), # choosing a different theme
  ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
  messages = FALSE,
  pairwise.comparisons = F, # display significant pairwise comparisons
  ggplot.component = list(theme(text = element_text(size = 15)))
)+ggplot2::theme(axis.title.y = element_text( size = 30, face = "plain"),
                 axis.text.x =element_text( size = 30, face = "plain"),
                 plot.title = element_text(size=30))+
  ggplot2::annotate("text", x=1.5, y=0.95*max(merge.data$cnv.sum), label= paste0("p=",signif(t9$p.value,3)),size=10)+
  scale_color_manual(values=c("#FF0033","#3300FF"))

t10.1<-merge.data[merge.data$three_region_Loss=="Loss",]
t10.2<-merge.data[merge.data$three_region_Loss=="No_loss",]
t10<-wilcox.test(t10.1$Loss.sum,t10.2$Loss.sum)
g10<-ggbetweenstats(
  data = merge.data,
  x = three_region_Loss,
  y = Loss.sum,
  notch = F, # show notched box plot
  mean.plotting = TRUE, # whether mean for each group is to be displayed
  mean.ci = F, # whether to display confidence interval for means
  type = "Parametric", # which type of test is to be run
  k = 2, # number of decimal places for statistical results
  outlier.tagging = F, # whether outliers need to be tagged
  xlab = "", # label for the x-axis variable
  ylab = "total N of arm losses", # label for the y-axis variable
  title = "17p13.1/9p21.3/3p14", # title text for the plot
  ggtheme = ggplot2::theme_classic(), # choosing a different theme
  ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
  messages = FALSE,
  pairwise.comparisons = F, # display significant pairwise comparisons
  ggplot.component = list(theme(text = element_text(size = 15)))
)+ggplot2::theme(axis.title.y = element_text( size = 30, face = "plain"),
                 axis.text.x =element_text( size = 30, face = "plain"),
                 plot.title = element_text(size=30))+
  ggplot2::annotate("text", x=1.5, y=0.95*max(merge.data$Loss.sum), label= paste0("p=",signif(t10$p.value,3)),size=10)+
  scale_color_manual(values=c("#FF0033","#3300FF"))

t11.1<-merge.data[merge.data$three_region_Loss=="Loss",]
t11.2<-merge.data[merge.data$three_region_Loss=="No_loss",]
t11<-wilcox.test(t11.1$amp.sum,t11.2$amp.sum)
g11<-ggbetweenstats(
  data = merge.data,
  x = three_region_Loss,
  y = amp.sum,
  notch = F, # show notched box plot
  mean.plotting = TRUE, # whether mean for each group is to be displayed
  mean.ci = F, # whether to display confidence interval for means
  type = "Parametric", # which type of test is to be run
  k = 2, # number of decimal places for statistical results
  outlier.tagging = F, # whether outliers need to be tagged
  xlab = "", # label for the x-axis variable
  ylab = "total N of arm gains", # label for the y-axis variable
  title = "", # title text for the plot
  ggtheme = ggplot2::theme_classic(), # choosing a different theme
  ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
  messages = FALSE,
  pairwise.comparisons = F, # display significant pairwise comparisons
  ggplot.component = list(theme(text = element_text(size = 15)))
)+ggplot2::theme(axis.title.y = element_text( size = 30, face = "plain"),
                 axis.text.x =element_text( size = 30, face = "plain"),
                 plot.title = element_text(size=30))+
  ggplot2::annotate("text", x=1.5, y=0.95*max(merge.data$amp.sum), label= paste0("p=",signif(t11$p.value,3)),size=10)+
  scale_color_manual(values=c("#FF0033","#3300FF"))

t12.1<-merge.data[merge.data$three_region_Loss=="Loss",]
t12.2<-merge.data[merge.data$three_region_Loss=="No_loss",]
t12<-wilcox.test(t12.1$cnv.sum,t12.2$cnv.sum)
g12<-ggbetweenstats(
  data = merge.data,
  x = three_region_Loss,
  y = cnv.sum,
  notch = F, # show notched box plot
  mean.plotting = TRUE, # whether mean for each group is to be displayed
  mean.ci = F, # whether to display confidence interval for means
  type = "Parametric", # which type of test is to be run
  k = 2, # number of decimal places for statistical results
  outlier.tagging = F, # whether outliers need to be tagged
  xlab = "", # label for the x-axis variable
  ylab = "total N of arm losses or gains", # label for the y-axis variable
  title = "", # title text for the plot
  ggtheme = ggplot2::theme_classic(), # choosing a different theme
  ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
  messages = FALSE,
  pairwise.comparisons = F, # display significant pairwise comparisons
  ggplot.component = list(theme(text = element_text(size = 15)))
)+ggplot2::theme(axis.title.y = element_text( size = 30, face = "plain"),
                 axis.text.x =element_text( size = 30, face = "plain"),
                 plot.title = element_text(size=30))+
  ggplot2::annotate("text", x=1.5, y=0.95*max(merge.data$cnv.sum), label= paste0("p=",signif(t12$p.value,3)),size=10)+
  scale_color_manual(values=c("#FF0033","#3300FF"))

figure<-ggarrange(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,
                  labels = c("A", "", "","B","","","C","","","D","",""),
                  ncol = 3, nrow = 4,
                  font.label = list(size = 30, color = "black"))
annotate_figure(figure,
                top = text_grob("", color = "black", face = "bold", size = 35))
}
#################################################################################################### FigS4
if(figS4==T){
##---Fig S4
library(Gviz)
#library(IdeoViz) # load then unload
library(RColorBrewer)
source("/Users/zhaox12/Desktop/Teresas_lab/project/cmd/function_cmd/IdeoViz/R/plotChromValuePair.R")
source("/Users/zhaox12/Desktop/Teresas_lab/project/cmd/function_cmd/IdeoViz/R/plotOnIdeo.R")
source("/Users/zhaox12/Desktop/Teresas_lab/project/cmd/function_cmd/IdeoViz/R/getIdeo.R")
data.header<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/0-others/Pan/Biomart.hg19.gene.info.txt",
                        sep="\t",header = T)
data.header.use<-data.header[data.header$Chromosome==9,c(1,2,3,7,8)]

scna<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/2.Nature_2010_CNV_TCGA/DataConcise_2019/HNSC/HNSC_all_data_by_genes_purity_rescale.txt",
                 sep="\t",header = T)
rownames(scna)<-scna$Gene.Symbol
scna.use<-scna[scna$Gene.Symbol %in% data.header.use$Gene.name,]
teresa.data<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/0-others/CCLE/Scotte-3.4-others/Teresa_data/final_use/TCGA.HNSC.alldata-Feb27.2021-Xin.txt",
                        sep="\t",header = T)
scna.use1<-scna.use[,colnames(scna.use) %in% teresa.data$sample]
scna.use1$sum.scna<-apply(scna.use1[,1:343],1,sum)
scna.use1$sum.scna<-(-1)*scna.use1$sum.scna
scna.use1$freq.loss<- rowSums(scna.use1[,1:343] < (-0.3))/343*100
scna.use1$gene<-rownames(scna.use1)

scna.use2<-scna.use1[,c("gene","freq.loss","sum.scna")]
merge.data<-merge(data.header.use,scna.use2,by.x="Gene.name",by.y="gene")
merge.data$Chromosome<-paste0("chr",merge.data$Chromosome)
rownames(merge.data)<-make.names(merge.data$Gene.name,unique = T)

gr1<-makeGRangesFromDataFrame(merge.data,
                              keep.extra.columns=T,
                              ignore.strand=T,
                              seqinfo=NULL,
                              seqnames.field=c("seqnames", "seqname",
                                               "chromosome", "chrom",
                                               "chr", "chromosome_name",
                                               "seqid"),
                              start.field="start_position",
                              end.field=c("end_position", "stop"),
                              starts.in.df.are.0based=FALSE)

library(BSgenome.Hsapiens.UCSC.hg19)
gene<-merge.data[merge.data$Gene.name=="CD274"|merge.data$Gene.name=="CDKN2A"|merge.data$Gene.name=="CDKN2B"|merge.data$Gene.name=="IFNA1"|merge.data$Gene.name=="IFNA2",]

gene$chromosome<-gene$Chromosome
gene$start<-gene$start_position
gene$end<-gene$end_position
gene$width<-gene$end_position-gene$start_position
gene$strand<-"+"
gene_feature<-"protein_coding"
gene$symbol<-gene$Gene.name

gene.use<-gene[,8:13]

chr <- as.character(unique(seqnames(gr1)))
gen<-"hg19"
atrack <- AnnotationTrack(gr1, name = "Gene",stacking="dense",background.title = "darkgrey")


grtrack <- GeneRegionTrack(gene.use, genome = gen, chromosome = chr,
                           name = "Potential target loss ",fill="#3300FF80",
                           transcriptAnnotation = "symbol",background.title = "#3300FF80",cex=10)
dtrack1<- DataTrack(data = merge.data$freq.loss, start = merge.data$start_position,
                    end = merge.data$end_position, chromosome = chr, genome = gen, 
                    name = "Percentage of patients (%)",background.title="#3300FF",col="#3300FF",type="a",cex=100)
dtrack2<- DataTrack(data = merge.data$sum.scna, start = merge.data$start_position,
                    end = merge.data$end_position, chromosome = chr, genome = gen, 
                    name = "gene.sum.loss",background.title="#3300FF80",col="#3300FF80",type="p",cex=1.5)

seg.data<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/0-others/CCLE/Scotte-3.4-others/new_analysis/HNSC-10kb.HPV_neg.mean.sum.txt",
                     sep="\t",header = T)
seg.data$sum<-(-1)*seg.data$sum
dtrack3<- DataTrack(data = seg.data$sum, start = seg.data$start_position,
                    end = seg.data$end_position, chromosome = chr, genome = gen, 
                    name = "segments.loss total",background.title="#3300FF80",col="#3300FF80",type="p",cex=1.5)

itrack <- IdeogramTrack(genome = gen, chromosome = chr)
strack <- SequenceTrack(Hsapiens, chromosome = chr)
plotTracks(list(itrack, atrack, grtrack,dtrack1), 
           from = 0, to = 49000000)
}
#################################################################################################### FigS5
if(figS5==T){
##---FigS5A/B
data.keep<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/0-others/CCLE/Scotte-3.4-others/Teresa_data/final_use/TCGA.HNSC.ABSOLUTE-Feb27.2021-Xin.txt")
data1<-data.keep[,c("stage.use","CD8A","CD3D","chr9p21.3.loss","Unweighted.SCNA.0.3.chrom.arm.scale","X9p.loss","TP53_Variant_Classification","SASP_epithelial_signature",
               "chr9p21.3_focalonly.loss")]
data1<-data1[data1$stage.use!=0,]
data1<-data1[!is.na(data1$stage.use),]
n_fun <- function(x){
  return(data.frame(y = 0.95*15.2,
                    label = paste0("N=",length(x))))
}
data1$stage.use <- as.factor(data1$stage.use)
data1.0<-data1
data1.1<-data1.0[data1.0$stage.use=="Stage I",]
data1.1$SCNA.level<-ifelse(data1.1$Unweighted.SCNA.0.3.chrom.arm.scale>quantile(data1.1$Unweighted.SCNA.0.3.chrom.arm.scale,0.65,na.rm=T),"High",
                           ifelse(data1.1$Unweighted.SCNA.0.3.chrom.arm.scale<quantile(data1.1$Unweighted.SCNA.0.3.chrom.arm.scale,0.35,na.rm = T), "Low",NA))
data1.2<-data1.0[data1.0$stage.use=="Stage II",]
data1.2$SCNA.level<-ifelse(data1.2$Unweighted.SCNA.0.3.chrom.arm.scale>quantile(data1.2$Unweighted.SCNA.0.3.chrom.arm.scale,0.65,na.rm=T),"High",
                           ifelse(data1.2$Unweighted.SCNA.0.3.chrom.arm.scale<quantile(data1.2$Unweighted.SCNA.0.3.chrom.arm.scale,0.35,na.rm = T),"Low",NA))
data1.3<-data1.0[data1.0$stage.use=="Stage III",]
data1.3$SCNA.level<-ifelse(data1.3$Unweighted.SCNA.0.3.chrom.arm.scale>quantile(data1.3$Unweighted.SCNA.0.3.chrom.arm.scale,0.65,na.rm=T),"High",
                           ifelse(data1.3$Unweighted.SCNA.0.3.chrom.arm.scale<quantile(data1.3$Unweighted.SCNA.0.3.chrom.arm.scale,0.35,na.rm = T),"Low",NA))
data1.4<-data1.0[data1.0$stage.use=="Stage IV",]
data1.4$SCNA.level<-ifelse(data1.4$Unweighted.SCNA.0.3.chrom.arm.scale>quantile(data1.4$Unweighted.SCNA.0.3.chrom.arm.scale,0.65,na.rm=T),"High",
                           ifelse(data1.4$Unweighted.SCNA.0.3.chrom.arm.scale<quantile(data1.4$Unweighted.SCNA.0.3.chrom.arm.scale,0.35,na.rm = T),"Low",NA))
data1.5<-rbind(data1.1,data1.2,data1.3,data1.4)
data1.5<-data1.5[!is.na(data1.5$SCNA.level),]
data1.5$SCNA.level <- factor(data1.5$SCNA.level, levels=c("Low", "High"))
stat.test2 <- data1.5 %>%
  group_by(stage.use) %>%
  wilcox_test(CD8A ~ SCNA.level) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test2 <- stat.test2 %>%
  add_xy_position(x = "stage.use", dodge = 0.8)
stat.test2$y.position<-12

g2<-ggplot(data1.5, aes(x=stage.use, y=CD8A)) + 
  geom_point(aes(fill = SCNA.level, color=SCNA.level), alpha=0.5,
             size = 2, shape = 21, position = position_jitterdodge()) +
  geom_boxplot(aes(fill=SCNA.level),alpha = 0.5,outlier.colour = NA)+
  stat_summary(fun.data = n_fun, geom = "text", 
               aes(group=SCNA.level),
               hjust = 0.5, position = position_dodge(0.8),size=5) +
  scale_color_manual(values=c("#FF0033","#3300FF"))+
  scale_fill_manual(values=c("#FF0033","#3300FF")) +
  theme_classic2()+
  ggtitle("Relationship between overall SCNA level and CD8 by Stage (TCGA)")+
  ylab("CD8")+xlab("")+scale_y_continuous(limits = c(-1, 16))+
  theme(axis.title.y = element_text( size = 15, face = "plain"),
        axis.text.x =element_text( size = 15, face = "plain"),
        plot.title = element_text(size=15))+
  theme(legend.position = c(0.08, 0.2))+
  stat_pvalue_manual(
    stat.test2,  label = "p",tip.length = 0,size=5
  )

stat.test2.1 <- data1.5 %>%
  group_by(stage.use) %>%
  wilcox_test(CD3D ~ SCNA.level) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test2.1 <- stat.test2.1 %>%
  add_xy_position(x = "stage.use", dodge = 0.8)
stat.test2.1$y.position<-12

g2.1<-ggplot(data1.5, aes(x=stage.use, y=CD3D)) + 
  geom_point(aes(fill = SCNA.level, color=SCNA.level), alpha=0.5,
             size = 2, shape = 21, position = position_jitterdodge()) +
  geom_boxplot(aes(fill=SCNA.level),alpha = 0.5,outlier.colour = NA)+
  stat_summary(fun.data = n_fun, geom = "text", 
               aes(group=SCNA.level),
               hjust = 0.5, position = position_dodge(0.8),size=5) +
  scale_color_manual(values=c("#FF0033","#3300FF"))+
  scale_fill_manual(values=c("#FF0033","#3300FF")) +
  theme_classic2()+
  ggtitle("Relationship between overall SCNA level and CD3 by Stage (TCGA)")+
  ylab("CD3")+xlab("")+scale_y_continuous(limits = c(-1, 16))+
  theme(axis.title.y = element_text( size = 15, face = "plain"),
        axis.text.x =element_text( size = 15, face = "plain"),
        plot.title = element_text(size=15))+
  theme(legend.position = c(0.08, 0.2))+
  stat_pvalue_manual(
    stat.test2.1,  label = "p",tip.length = 0,size=5
  )
##---FigS5C/D
rm(list=ls())
library(ggstatsplot)
source("/Users/zhaox12/Desktop/Teresas_lab/project/cmd/function_cmd/box_plot/Script.4group.color_specific_ggstats_plot.R")
setwd("/Users/zhaox12/Desktop/Teresas_lab/project/0-others/CCLE/Scotte-3.4-others/new_analysis/cybersort/neg/plot/")
cyber<-read.delim("../CIBERSORT-Results.txt",sep="\t",header = T)
merge.data<-merge(data.keep[,c("sample","stage.use","Unweighted.SCNA.0.3.chrom.arm.scale","TP53_Variant_Classification")],cyber,by.x="sample",by.y="Mixture")
data1<-merge.data
data1<-data1[data1$stage.use!=0,]
data1<-data1[!is.na(data1$stage.use),]
data1$stage.use <- as.factor(data1$stage.use)
data1.0<-data1
data1.1<-data1.0[data1.0$stage.use=="Stage I",]
data1.1$SCNA.level<-ifelse(data1.1$Unweighted.SCNA.0.3.chrom.arm.scale>quantile(data1.1$Unweighted.SCNA.0.3.chrom.arm.scale,0.5,na.rm=T),"High",
                           ifelse(data1.1$Unweighted.SCNA.0.3.chrom.arm.scale<quantile(data1.1$Unweighted.SCNA.0.3.chrom.arm.scale,0.5,na.rm = T), "Low",NA))
data1.2<-data1.0[data1.0$stage.use=="Stage II",]
data1.2$SCNA.level<-ifelse(data1.2$Unweighted.SCNA.0.3.chrom.arm.scale>quantile(data1.2$Unweighted.SCNA.0.3.chrom.arm.scale,0.5,na.rm=T),"High",
                           ifelse(data1.2$Unweighted.SCNA.0.3.chrom.arm.scale<quantile(data1.2$Unweighted.SCNA.0.3.chrom.arm.scale,0.5,na.rm = T),"Low",NA))
data1.3<-data1.0[data1.0$stage.use=="Stage III",]
data1.3$SCNA.level<-ifelse(data1.3$Unweighted.SCNA.0.3.chrom.arm.scale>quantile(data1.3$Unweighted.SCNA.0.3.chrom.arm.scale,0.5,na.rm=T),"High",
                           ifelse(data1.3$Unweighted.SCNA.0.3.chrom.arm.scale<quantile(data1.3$Unweighted.SCNA.0.3.chrom.arm.scale,0.5,na.rm = T),"Low",NA))
data1.4<-data1.0[data1.0$stage.use=="Stage IV",]
data1.4$SCNA.level<-ifelse(data1.4$Unweighted.SCNA.0.3.chrom.arm.scale>quantile(data1.4$Unweighted.SCNA.0.3.chrom.arm.scale,0.5,na.rm=T),"High",
                           ifelse(data1.4$Unweighted.SCNA.0.3.chrom.arm.scale<quantile(data1.4$Unweighted.SCNA.0.3.chrom.arm.scale,0.5,na.rm = T),"Low",NA))
data1.5<-rbind(data1.1,data1.2,data1.3,data1.4)
data1.5<-data1.5[!is.na(data1.5$SCNA.level),]
data1.5$SCNA.level <- factor(data1.5$SCNA.level, levels=c("Low", "High"))
n_fun <- function(x){
  return(data.frame(y = 0.95*0.5,
                    label = paste0("N=",length(x))))
}
for (i in 1:22){
  #i<-13
  D<-data1.5[,c(2,31,i+3)]
  D$events<-D[,3]
  stat.test2 <- D %>%
    group_by(stage.use) %>%
    wilcox_test(events ~ SCNA.level) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  stat.test2 <- stat.test2 %>%
    add_xy_position(x = "stage.use", dodge = 0.8)
  stat.test2$y.position<-0.4
  
  g2<-ggplot(D, aes(x=stage.use, y=events)) + 
    geom_point(aes(fill = SCNA.level, color=SCNA.level), alpha=0.5,
               size = 2, shape = 21, position = position_jitterdodge()) +
    geom_boxplot(aes(fill=SCNA.level),alpha = 0.5,outlier.colour = NA)+
    stat_summary(fun.data = n_fun, geom = "text", 
                 aes(group=SCNA.level),
                 hjust = 0.5, position = position_dodge(0.8),size=5) +
    scale_color_manual(values=c("#FF0033","#3300FF"))+
    scale_fill_manual(values=c("#FF0033","#3300FF")) +
    theme_classic2()+
    ggtitle("Relationship between overall Cell Fraction and SCNA by Stage (TCGA)",subtitle = colnames(D)[3])+
    ylab("Fraction")+xlab("")+scale_y_continuous(limits = c(-0.1, 0.5))+
    theme(axis.title.y = element_text( size = 15, face = "plain"),
          axis.text.x =element_text( size = 15, face = "plain"),
          plot.title = element_text(size=15))+
    theme(legend.position = "bottom")+
    stat_pvalue_manual(
      stat.test2,  label = "p",tip.length = 0,size=5
    )
}
}
#################################################################################################### FigS6
##---FigS6 should be exact same with Fig4, using CD3D instead of CD8A
#################################################################################################### FigS7
if(figS7==T){
##---FigS7A from GSEA
##---FigS7B~E
data.keep<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/0-others/CCLE/Scotte-3.4-others/Teresa_data/final_use/TCGA.HNSC.alldata-Feb27.2021-Xin.txt",
                      sep="\t",header = T)
  D<-data.keep
  data1<-D[,c("CD274","CD8A","KEGG_JAK_STAT_signature","KEGG_Cytokine_Cytokine_receptor_interaction_signature","chr9p21.3.loss","chr3p14.loss","X9p_cnv.loss","X3p_cnv.loss",
              "TP53_Variant_Classification","SASP_epithelial_signature","CDKN2A.loss","chr9p21.3.focal.loss","TNFA_SIGNALING_VIA_NFKB_signature","chr9p24.1.loss","chr9p24.1.focal.loss")]
  colnames(data1)[7]<-"Chr9p_cnv.loss"
  colnames(data1)[8]<-"Chr3p_cnv.loss"
  data1$chr9p21.3.loss<-ifelse(data1$chr9p21.3.loss==1,"Loss","No_loss")
  data1$chr9p21.3.loss<-factor(data1$chr9p21.3.loss,levels = c("No_loss","Loss"))
  data1$chr9p24.1.loss<-ifelse(data1$chr9p24.1.loss==1,"Loss","No_loss")
  data1$chr9p24.1.loss<-factor(data1$chr9p24.1.loss,levels = c("No_loss","Loss"))
  data1$CDKN2A.loss<-ifelse(data1$CDKN2A.loss==1,"Loss","No_loss")
  data1$CDKN2A.loss<-factor(data1$CDKN2A.loss,levels = c("No_loss","Loss"))
  data.rep<-data1
  data.rep<-data.rep[!is.na(data.rep$Chr9p_cnv.loss),]
  data.rep$Chr9p_cnv.loss<-ifelse(data.rep$Chr9p_cnv.loss==1,"Loss","No_loss")
  data.rep$Chr9p_cnv.loss<-factor(data.rep$Chr9p_cnv.loss,levels = c("No_loss","Loss"))
  data.rep1<-data1
  data.rep1<-data.rep1[!is.na(data.rep1$chr9p21.3.focal.loss),]
  data.rep1$chr9p21.3.focal.loss<-ifelse(data.rep1$chr9p21.3.focal.loss==1,"Loss","No_loss")
  data.rep1$chr9p21.3.focal.loss<-factor(data.rep1$chr9p21.3.focal.loss,levels = c("No_loss","Loss"))
  data.rep2<-data1
  data.rep2<-data.rep2[!is.na(data.rep2$chr9p24.1.focal.loss),]
  data.rep2$chr9p24.1.focal.loss<-ifelse(data.rep2$chr9p24.1.focal.loss==1,"Loss","No_loss")
  data.rep2$chr9p24.1.focal.loss<-factor(data.rep2$chr9p24.1.focal.loss,levels = c("No_loss","Loss"))
  
  
  stat.test1 <- data1 %>%
    wilcox_test(SASP_epithelial_signature ~ chr9p21.3.loss) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  stat.test1$y.position<-0.5
  stat.test2 <- data.rep %>%
    wilcox_test(SASP_epithelial_signature ~ Chr9p_cnv.loss) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  stat.test2$y.position<-0.5
  stat.test2.1 <- data.rep1 %>%
    wilcox_test(SASP_epithelial_signature ~ chr9p21.3.focal.loss) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  stat.test2.1$y.position<-0.5
  
  
  n_fun <- function(x){
    return(data.frame(y = 0.95*0.65,
                      label = paste0("N=",length(x))))
  }
  g1<-ggplot(data1, aes(x=chr9p21.3.loss, y=SASP_epithelial_signature)) + 
    geom_point(aes(fill = chr9p21.3.loss, color=chr9p21.3.loss), alpha=0.5,
               size = 2, shape = 21, position = position_jitterdodge()) +
    geom_boxplot(aes(fill = chr9p21.3.loss),alpha = 0.5)+
    stat_summary(fun.data = n_fun, geom = "text", 
                 aes(group=chr9p21.3.loss),
                 hjust = 0.5, position = position_dodge(0.8),size=7) +
    scale_fill_manual(values=c("#FF0033","#3300FF")) +
    scale_color_manual(values=c("#FF0033","#3300FF"))+
    theme_classic2()+
    ggtitle("9p21.3 loss (arm or focal) and SASP")+
    ylab("SASP_enrichment")+xlab("")+scale_y_continuous(limits = c(-0.5, 0.65))+
    theme(axis.title.y = element_text( size = 25, face = "plain"),
          axis.text.x =element_text( size = 25, face = "plain"),
          plot.title = element_text(size=20),
          legend.title=element_text(size=16), 
          legend.text=element_text(size=15))+
    theme(legend.position = c(0.18, 0.22))+
    stat_pvalue_manual(
      stat.test1,  label = "p.adj",tip.length = 0,size = 7
    )
  
  g2<-ggplot(data.rep, aes(x=Chr9p_cnv.loss, y=SASP_epithelial_signature)) + 
    geom_point(aes(fill = Chr9p_cnv.loss, color=Chr9p_cnv.loss), alpha=0.5,
               size = 2, shape = 21, position = position_jitterdodge()) +
    geom_boxplot(aes(fill = Chr9p_cnv.loss),alpha = 0.5)+
    stat_summary(fun.data = n_fun, geom = "text", 
                 aes(group=Chr9p_cnv.loss),
                 hjust = 0.5, position = position_dodge(0.8),size=7) +
    scale_fill_manual(values=c("#FF0033","#3300FF")) +
    scale_color_manual(values=c("#FF0033","#3300FF"))+
    theme_classic2()+
    ggtitle("9p arm loss and SASP")+
    ylab("SASP_enrichment")+xlab("")+scale_y_continuous(limits = c(-0.5, 0.65))+
    theme(text = element_text(size = 15),
          axis.title.y = element_text( size = 25, face = "plain"),
          axis.text.x =element_text( size = 25, face = "plain"),
          plot.title = element_text(size=20),
          legend.title=element_text(size=16), 
          legend.text=element_text(size=15))+
    theme(legend.position = c(0.18, 0.22))+
    stat_pvalue_manual(
      stat.test2,  label = "p.adj",tip.length = 0,size=7
    )
  g2.1<-ggplot(data.rep1, aes(x=chr9p21.3.focal.loss, y=SASP_epithelial_signature)) + 
    geom_point(aes(fill = chr9p21.3.focal.loss, color=chr9p21.3.focal.loss), alpha=0.5,
               size = 2, shape = 21, position = position_jitterdodge()) +
    geom_boxplot(aes(fill = chr9p21.3.focal.loss),alpha = 0.5)+
    stat_summary(fun.data = n_fun, geom = "text", 
                 aes(group=chr9p21.3.focal.loss),
                 hjust = 0.5, position = position_dodge(0.8),size=7) +
    scale_fill_manual(values=c("#FF0033","#3300FF")) +
    scale_color_manual(values=c("#FF0033","#3300FF"))+
    theme_classic2()+
    ggtitle("9p21.3 focal loss and SASP")+
    ylab("SASP_enrichment")+xlab("")+scale_y_continuous(limits = c(-0.5, 0.65))+
    theme(text = element_text(size = 15),
          axis.title.y = element_text( size = 25, face = "plain"),
          axis.text.x =element_text( size = 25, face = "plain"),
          plot.title = element_text(size=20),
          legend.title=element_text(size=16), 
          legend.text=element_text(size=15))+
    theme(legend.position = c(0.22, 0.22))+
    stat_pvalue_manual(
      stat.test2.1,  label = "p.adj",tip.length = 0,size=7
    )
  n_fun1 <- function(x){
    return(data.frame(y = 0.95*0.3,
                      label = paste0("N=",length(x))))
  }
  stat.test3 <- data1 %>%
    wilcox_test(KEGG_JAK_STAT_signature ~ chr9p21.3.loss) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  stat.test3$y.position<-0.2
  
  g3<-ggplot(data1, aes(x=chr9p21.3.loss, y=KEGG_JAK_STAT_signature)) + 
    geom_point(aes(fill = chr9p21.3.loss, color=chr9p21.3.loss), alpha=0.5,
               size = 2, shape = 21, position = position_jitterdodge()) +
    geom_boxplot(aes(fill = chr9p21.3.loss),alpha = 0.5)+
    stat_summary(fun.data = n_fun1, geom = "text", 
                 aes(group=chr9p21.3.loss),
                 hjust = 0.5, position = position_dodge(0.8),size=7) +
    scale_fill_manual(values=c("#FF0033","#3300FF")) +
    scale_color_manual(values=c("#FF0033","#3300FF"))+
    theme_classic2()+
    ggtitle("9p21.3 loss (arm or focal) and JAK-STAT")+
    ylab("JAK-STAT_enrichment")+xlab("")+scale_y_continuous(limits = c(-0.2, 0.3))+
    theme(axis.title.y = element_text( size = 25, face = "plain"),
          axis.text.x =element_text(size = 25, face = "plain"),
          plot.title = element_text(size=20),
          legend.title=element_text(size=16), 
          legend.text=element_text(size=15))+
    theme(legend.position = c(0.18, 0.22))+
    stat_pvalue_manual(
      stat.test3,  label = "p.adj",tip.length = 0,size = 7
    )
  
  stat.test4 <- data.rep %>%
    wilcox_test(KEGG_JAK_STAT_signature ~ Chr9p_cnv.loss) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  stat.test4$y.position<-0.2
  g4<-ggplot(data.rep, aes(x=Chr9p_cnv.loss, y=KEGG_JAK_STAT_signature)) + 
    geom_point(aes(fill = Chr9p_cnv.loss, color=Chr9p_cnv.loss), alpha=0.5,
               size = 2, shape = 21, position = position_jitterdodge()) +
    geom_boxplot(aes(fill = Chr9p_cnv.loss),alpha = 0.5)+
    stat_summary(fun.data = n_fun1, geom = "text", 
                 aes(group=Chr9p_cnv.loss),
                 hjust = 0.5, position = position_dodge(0.8),size=7) +
    scale_fill_manual(values=c("#FF0033","#3300FF")) +
    scale_color_manual(values=c("#FF0033","#3300FF"))+
    theme_classic2()+
    ggtitle("9p arm loss and JAK-STAT")+
    ylab("JAK-STAT_enrichment")+xlab("")+scale_y_continuous(limits = c(-0.2, 0.3))+
    theme(axis.title.y = element_text( size = 25, face = "plain"),
          axis.text.x =element_text(size = 25, face = "plain"),
          plot.title = element_text(size=20),
          legend.title=element_text(size=16), 
          legend.text=element_text(size=15))+
    theme(legend.position = c(0.18, 0.22))+
    stat_pvalue_manual(
      stat.test4,  label = "p.adj",tip.length = 0,size = 7
    )
  
  stat.test4.1 <- data.rep1 %>%
    wilcox_test(KEGG_JAK_STAT_signature ~ chr9p21.3.focal.loss) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  stat.test4.1$y.position<-0.2
  g4.1<-ggplot(data.rep1, aes(x=chr9p21.3.focal.loss, y=KEGG_JAK_STAT_signature)) + 
    geom_point(aes(fill = chr9p21.3.focal.loss, color=chr9p21.3.focal.loss), alpha=0.5,
               size = 2, shape = 21, position = position_jitterdodge()) +
    geom_boxplot(aes(fill = chr9p21.3.focal.loss),alpha = 0.5)+
    stat_summary(fun.data = n_fun1, geom = "text", 
                 aes(group=chr9p21.3.focal.loss),
                 hjust = 0.5, position = position_dodge(0.8),size=7) +
    scale_fill_manual(values=c("#FF0033","#3300FF")) +
    scale_color_manual(values=c("#FF0033","#3300FF"))+
    theme_classic2()+
    ggtitle("9p21.3 focal loss and JAK-STAT")+
    ylab("JAK-STAT_enrichment")+xlab("")+scale_y_continuous(limits = c(-0.2, 0.3))+
    theme(axis.title.y = element_text( size = 25, face = "plain"),
          axis.text.x =element_text(size = 25, face = "plain"),
          plot.title = element_text(size=20),
          legend.title=element_text(size=16), 
          legend.text=element_text(size=15))+
    theme(legend.position = c(0.22, 0.22))+
    stat_pvalue_manual(
      stat.test4.1,  label = "p.adj",tip.length = 0,size = 7
    )
  
  stat.test5 <- data1 %>%
    wilcox_test(KEGG_Cytokine_Cytokine_receptor_interaction_signature ~ chr9p21.3.loss) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  stat.test5$y.position<-0.2
  g5<-ggplot(data1, aes(x=chr9p21.3.loss, y=KEGG_Cytokine_Cytokine_receptor_interaction_signature)) + 
    geom_point(aes(fill = chr9p21.3.loss, color=chr9p21.3.loss), alpha=0.5,
               size = 2, shape = 21, position = position_jitterdodge()) +
    geom_boxplot(aes(fill = chr9p21.3.loss),alpha = 0.5)+
    stat_summary(fun.data = n_fun1, geom = "text", 
                 aes(group=chr9p21.3.loss),
                 hjust = 0.5, position = position_dodge(0.8),size=7) +
    scale_fill_manual(values=c("#FF0033","#3300FF")) +
    scale_color_manual(values=c("#FF0033","#3300FF"))+
    theme_classic2()+
    ggtitle("9p21.3 loss (arm or focal) and Cytokine")+
    ylab("Cytokine_enrichment")+xlab("")+scale_y_continuous(limits = c(-0.4, 0.3))+
    theme(axis.title.y = element_text( size = 25, face = "plain"),
          axis.text.x =element_text(size = 25, face = "plain"),
          plot.title = element_text(size=20),
          legend.title=element_text(size=16), 
          legend.text=element_text(size=15))+
    theme(legend.position = c(0.18, 0.22))+
    stat_pvalue_manual(
      stat.test5,  label = "p.adj",tip.length = 0,size = 7
    )
  
  stat.test6 <- data.rep %>%
    wilcox_test(KEGG_Cytokine_Cytokine_receptor_interaction_signature ~ Chr9p_cnv.loss) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  stat.test6$y.position<-0.2
  g6<-ggplot(data.rep, aes(x=Chr9p_cnv.loss, y=KEGG_Cytokine_Cytokine_receptor_interaction_signature)) + 
    geom_point(aes(fill = Chr9p_cnv.loss, color=Chr9p_cnv.loss), alpha=0.5,
               size = 2, shape = 21, position = position_jitterdodge()) +
    geom_boxplot(aes(fill = Chr9p_cnv.loss),alpha = 0.5)+
    stat_summary(fun.data = n_fun1, geom = "text", 
                 aes(group=Chr9p_cnv.loss),
                 hjust = 0.5, position = position_dodge(0.8),size=7) +
    scale_fill_manual(values=c("#FF0033","#3300FF")) +
    scale_color_manual(values=c("#FF0033","#3300FF"))+
    theme_classic2()+
    ggtitle("9p arm loss and Cytokine")+
    ylab("Cytokine_enrichment")+xlab("")+scale_y_continuous(limits = c(-0.4, 0.3))+
    theme(axis.title.y = element_text( size = 25, face = "plain"),
          axis.text.x =element_text(size = 25, face = "plain"),
          plot.title = element_text(size=20),
          legend.title=element_text(size=16), 
          legend.text=element_text(size=15))+
    theme(legend.position = c(0.18, 0.22))+
    stat_pvalue_manual(
      stat.test6,  label = "p.adj",tip.length = 0,size = 7
    )
  
  stat.test6.1 <- data.rep1 %>%
    wilcox_test(KEGG_Cytokine_Cytokine_receptor_interaction_signature ~ chr9p21.3.focal.loss) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  stat.test6.1$y.position<-0.2
  g6.1<-ggplot(data.rep1, aes(x=chr9p21.3.focal.loss, y=KEGG_Cytokine_Cytokine_receptor_interaction_signature)) + 
    geom_point(aes(fill = chr9p21.3.focal.loss, color=chr9p21.3.focal.loss), alpha=0.5,
               size = 2, shape = 21, position = position_jitterdodge()) +
    geom_boxplot(aes(fill = chr9p21.3.focal.loss),alpha = 0.5)+
    stat_summary(fun.data = n_fun1, geom = "text", 
                 aes(group=chr9p21.3.focal.loss),
                 hjust = 0.5, position = position_dodge(0.8),size=7) +
    scale_fill_manual(values=c("#FF0033","#3300FF")) +
    scale_color_manual(values=c("#FF0033","#3300FF"))+
    theme_classic2()+
    ggtitle("9p21.3 focal loss and Cytokine")+
    ylab("Cytokine_enrichment")+xlab("")+scale_y_continuous(limits = c(-0.4, 0.3))+
    theme(axis.title.y = element_text( size = 25, face = "plain"),
          axis.text.x =element_text(size = 25, face = "plain"),
          plot.title = element_text(size=20),
          legend.title=element_text(size=16), 
          legend.text=element_text(size=15))+
    theme(legend.position = c(0.22, 0.22))+
    stat_pvalue_manual(
      stat.test6.1,  label = "p.adj",tip.length = 0,size = 7
    )
  
  stat.test7 <- data1 %>%
    wilcox_test(TNFA_SIGNALING_VIA_NFKB_signature ~ chr9p21.3.loss) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  stat.test7$y.position<-0.55
  g7<-ggplot(data1, aes(x=chr9p21.3.loss, y=TNFA_SIGNALING_VIA_NFKB_signature)) + 
    geom_point(aes(fill = chr9p21.3.loss, color=chr9p21.3.loss), alpha=0.5,
               size = 2, shape = 21, position = position_jitterdodge()) +
    geom_boxplot(aes(fill = chr9p21.3.loss),alpha = 0.5)+
    stat_summary(fun.data = n_fun, geom = "text", 
                 aes(group=chr9p21.3.loss),
                 hjust = 0.5, position = position_dodge(0.8),size=7) +
    scale_fill_manual(values=c("#FF0033","#3300FF")) +
    scale_color_manual(values=c("#FF0033","#3300FF"))+
    theme_classic2()+
    ggtitle("9p21.3 loss (arm or focal) and TNFA")+
    ylab("TNFA_enrichment")+xlab("")+scale_y_continuous(limits = c(0, 0.65))+
    theme(axis.title.y = element_text( size = 25, face = "plain"),
          axis.text.x =element_text(size = 25, face = "plain"),
          plot.title = element_text(size=20),
          legend.title=element_text(size=16), 
          legend.text=element_text(size=15))+
    theme(legend.position = c(0.18, 0.22))+
    stat_pvalue_manual(
      stat.test7,  label = "p.adj",tip.length = 0,size = 7
    )
  
  stat.test8 <- data.rep %>%
    wilcox_test(TNFA_SIGNALING_VIA_NFKB_signature ~ Chr9p_cnv.loss) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  stat.test8$y.position<-0.55
  g8<-ggplot(data.rep, aes(x=Chr9p_cnv.loss, y=TNFA_SIGNALING_VIA_NFKB_signature)) + 
    geom_point(aes(fill = Chr9p_cnv.loss, color=Chr9p_cnv.loss), alpha=0.5,
               size = 2, shape = 21, position = position_jitterdodge()) +
    geom_boxplot(aes(fill = Chr9p_cnv.loss),alpha = 0.5)+
    stat_summary(fun.data = n_fun, geom = "text", 
                 aes(group=Chr9p_cnv.loss),
                 hjust = 0.5, position = position_dodge(0.8),size=7) +
    scale_fill_manual(values=c("#FF0033","#3300FF")) +
    scale_color_manual(values=c("#FF0033","#3300FF"))+
    theme_classic2()+
    ggtitle("9p arm loss and TNFA")+
    ylab("TNFA_enrichment")+xlab("")+scale_y_continuous(limits = c(0, 0.65))+
    theme(axis.title.y = element_text( size = 25, face = "plain"),
          axis.text.x =element_text(size = 25, face = "plain"),
          plot.title = element_text(size=20),
          legend.title=element_text(size=16), 
          legend.text=element_text(size=15))+
    theme(legend.position = c(0.18, 0.22))+
    stat_pvalue_manual(
      stat.test8,  label = "p.adj",tip.length = 0,size = 7
    )
  
  stat.test8.1 <- data.rep1 %>%
    wilcox_test(TNFA_SIGNALING_VIA_NFKB_signature ~ chr9p21.3.focal.loss) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  stat.test8.1$y.position<-0.55
  g8.1<-ggplot(data.rep1, aes(x=chr9p21.3.focal.loss, y=TNFA_SIGNALING_VIA_NFKB_signature)) + 
    geom_point(aes(fill = chr9p21.3.focal.loss, color=chr9p21.3.focal.loss), alpha=0.5,
               size = 2, shape = 21, position = position_jitterdodge()) +
    geom_boxplot(aes(fill = chr9p21.3.focal.loss),alpha = 0.5)+
    stat_summary(fun.data = n_fun, geom = "text", 
                 aes(group=chr9p21.3.focal.loss),
                 hjust = 0.5, position = position_dodge(0.8),size=7) +
    scale_fill_manual(values=c("#FF0033","#3300FF")) +
    scale_color_manual(values=c("#FF0033","#3300FF"))+
    theme_classic2()+
    ggtitle("9p21.3 focal loss and TNFA")+
    ylab("TNFA_enrichment")+xlab("")+scale_y_continuous(limits = c(0, 0.65))+
    theme(axis.title.y = element_text( size = 25, face = "plain"),
          axis.text.x =element_text(size = 25, face = "plain"),
          plot.title = element_text(size=20),
          legend.title=element_text(size=16), 
          legend.text=element_text(size=15))+
    theme(legend.position = c(0.22, 0.22))+
    stat_pvalue_manual(
      stat.test8.1,  label = "p.adj",tip.length = 0,size = 7
    )
figure<-ggarrange(
  g2,g4,g6,g8,g1,g3,g5,g7,g2.1,g4.1,g6.1,g8.1,
  nrow = 5,ncol=4, font.label = list(size = 30, color = "black"))
}
#################################################################################################### FigS8
if(figS8==T){
##---FigS8
library(ggplot2)
library(ggpubr)
library(rstatix)
data.keep<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/0-others/CCLE/Scotte-3.4-others/Teresa_data/final_use/TCGA.HNSC.alldata-Feb27.2021-Xin.txt",
                      sep="\t",header = T)
#TCGA
data<-data.keep
merge.data<-data[,c("sample","X9p21.3","X9p_cnv","IFNA_gene_family_signature")]
merge.data$chr9p21.3.loss<-ifelse(merge.data$X9p21.3<(-0.3),"Loss","No_loss")
merge.data$chr9p21.3.loss<-factor(merge.data$chr9p21.3.loss,levels = c("No_loss","Loss"))
merge.data$arm9p.loss<-ifelse(merge.data$X9p_cnv<(-0.3),"Loss","No_loss")
merge.data$arm9p.loss<-factor(merge.data$arm9p.loss,levels = c("No_loss","Loss"))

#plot
set.seed(0)
data.use<-merge.data
data.use$event<-data.use[,4]
stat.test1 <- data.use %>%
  wilcox_test(event ~ chr9p21.3.loss) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test1$y.position<--0.2

n_fun <- function(x){
  return(data.frame(y = 0.99*0,
                    label = paste0("N=",length(x))))
}
g1<-ggplot(data.use, aes(x=chr9p21.3.loss, y=event)) + 
  geom_boxplot(aes(fill = chr9p21.3.loss),alpha = 0.5)+
  geom_point(aes(fill = chr9p21.3.loss, color=chr9p21.3.loss), alpha=0.5,
             size = 2, shape = 21, position = position_jitterdodge()) +
  stat_summary(fun.data = n_fun, geom = "text", 
               aes(group=chr9p21.3.loss),
               hjust = 0.5, position = position_dodge(0.8),size=5) +
  scale_fill_manual(values=c("#FF0033","#3300FF")) +
  scale_color_manual(values=c("#FF0033","#3300FF"))+
  theme_classic2()+
  ggtitle("9p21.3 loss and IFNA")+
  ylab("IFNA_enrichment")+xlab("")+scale_y_continuous(limits = c(-1, 0))+
  theme(axis.title.y = element_text( size = 18, face = "plain"),
        axis.text.x =element_text( size = 18, face = "plain"),
        plot.title = element_text(size=18))+
  theme(legend.position = c(0.2, 0.2))+
  stat_pvalue_manual(
    stat.test1,  label = "p",tip.length = 0,size =5
  )

stat.test2 <- data.use %>%
  wilcox_test(event ~ arm9p.loss) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test2$y.position<--0.2

g2<-ggplot(data.use, aes(x=arm9p.loss, y=event)) + 
  geom_boxplot(aes(fill = arm9p.loss),alpha = 0.5)+
  geom_point(aes(fill = arm9p.loss, color=arm9p.loss), alpha=0.5,
             size = 2, shape = 21, position = position_jitterdodge()) +
  stat_summary(fun.data = n_fun, geom = "text", 
               aes(group=arm9p.loss),
               hjust = 0.5, position = position_dodge(0.8),size=5) +
  scale_fill_manual(values=c("#FF0033","#3300FF")) +
  scale_color_manual(values=c("#FF0033","#3300FF"))+
  theme_classic2()+
  ggtitle("9p arm loss and IFNA (TCGA)")+
  ylab("IFNA_enrichment")+xlab("")+scale_y_continuous(limits = c(-1, 0))+
  theme(axis.title.y = element_text( size = 18, face = "plain"),
        axis.text.x =element_text( size = 18, face = "plain"),
        plot.title = element_text(size=18))+
  theme(legend.position = c(0.2, 0.2))+
  stat_pvalue_manual(
    stat.test2,  label = "p",tip.length = 0,size =5
  )
com<-ggarrange(g2,g1,nrow=2,ncol = 1)

##### cell lines
datakeep<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/0-others/CCLE/Scotte-3.4-others/Teresa_data/final_use/CCLE.HNSC.alldata-Feb27.2021-Xin.txt")
data<-datakeep[,c("sample","chr9p21.3.loss","X9p_cnv.loss","IFNA_gene_family_signature")]
data$chr9p21.3.loss<-ifelse(data$chr9p21.3.loss==1,"Loss","No_loss")
data$chr9p21.3.loss<-factor(data$chr9p21.3.loss,levels = c("No_loss","Loss"))
data$arm9p.loss<-ifelse(data$X9p_cnv.loss==1,"Loss","No_loss")
data$arm9p.loss<-factor(data$arm9p.loss,levels = c("No_loss","Loss"))

#plot
set.seed(0)
data.use<-data
data.use$event<-data.use[,4]
stat.test3 <- data.use %>%
  wilcox_test(event ~ chr9p21.3.loss,alternative = "less") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test3$y.position<--0.2

n_fun <- function(x){
  return(data.frame(y = 0.99*0,
                    label = paste0("N=",length(x))))
}
g3<-ggplot(data.use, aes(x=chr9p21.3.loss, y=event)) + 
  geom_boxplot(aes(fill = chr9p21.3.loss),alpha = 0.5)+
  geom_point(aes(fill = chr9p21.3.loss, color=chr9p21.3.loss), alpha=0.5,
             size = 2, shape = 21, position = position_jitterdodge()) +
  stat_summary(fun.data = n_fun, geom = "text", 
               aes(group=chr9p21.3.loss),
               hjust = 0.5, position = position_dodge(0.8),size=5) +
  scale_fill_manual(values=c("#FF0033","#3300FF")) +
  scale_color_manual(values=c("#FF0033","#3300FF"))+
  theme_classic2()+
  ggtitle("9p21.3 loss and IFNA (cell lines)")+
  ylab("IFNA_enrichment")+xlab("")+scale_y_continuous(limits = c(-1, 0))+
  theme(axis.title.y = element_text( size = 18, face = "plain"),
        axis.text.x =element_text( size = 18, face = "plain"),
        plot.title = element_text(size=18))+
  theme(legend.position = c(0.2, 0.2))+
  stat_pvalue_manual(
    stat.test3,  label = "p",tip.length = 0,size =5
  )

stat.test4 <- data.use %>%
  wilcox_test(event ~ arm9p.loss,alternative = "less") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test4$y.position<--0.2

g4<-ggplot(data.use, aes(x=arm9p.loss, y=event)) + 
  geom_boxplot(aes(fill = arm9p.loss),alpha = 0.5)+
  geom_point(aes(fill = arm9p.loss, color=arm9p.loss), alpha=0.5,
             size = 2, shape = 21, position = position_jitterdodge()) +
  stat_summary(fun.data = n_fun, geom = "text", 
               aes(group=arm9p.loss),
               hjust = 0.5, position = position_dodge(0.8),size=5) +
  scale_fill_manual(values=c("#FF0033","#3300FF")) +
  scale_color_manual(values=c("#FF0033","#3300FF"))+
  theme_classic2()+
  ggtitle("9p arm loss and IFNA (cell lines)")+
  ylab("IFNA_enrichment")+xlab("")+scale_y_continuous(limits = c(-1, 0))+
  theme(axis.title.y = element_text( size = 18, face = "plain"),
        axis.text.x =element_text( size = 18, face = "plain"),
        plot.title = element_text(size=18))+
  theme(legend.position = c(0.2, 0.2))+
  stat_pvalue_manual(
    stat.test4,  label = "p",tip.length = 0,size =5
  )
com1<-ggarrange(g4,g3,nrow=2,ncol = 1)
}
#################################################################################################### FigS9
if(figS9==T){
##---FigS9
library(ggplot2)
library(ggpubr)
library(rstatix)
#TCGA
data.keep<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/0-others/CCLE/Scotte-3.4-others/Teresa_data/final_use/TCGA.HNSC.alldata-Feb27.2021-Xin.txt",
                      sep="\t",header = T)
data<-data.keep
gene<-c("CCL2","TNFSF10","CCL24","CCL20","CXCL3","TNFRSF11B")
gene<-factor(gene,levels =gene)

exp<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/2.Nature_2010_CNV_TCGA/DataConcise_2019/HNSC/HNSC_rnaseq_raw_counts.txt",sep="\t",header = T)
exp1<-as.data.frame(t(exp))
colnames(exp1)<-substr(colnames(exp1), start = 1, stop = 12)
colnames(exp1)<-gsub("\\-",".",colnames(exp1))
rownames(exp1)<-make.names(gsub("\\..*","",rownames(exp1)),unique = T)
exp2<-exp1[rownames(exp1) %in% gene, colnames(exp1) %in% data$sample]

exp3<-as.data.frame(t(log2(exp2+1)))
exp3$sample<-rownames(exp3)
merge.data<-merge(data[,c("sample","chr9p21.3.loss","X9p_cnv.loss")],exp3,by="sample")
merge.data$chr9p21.3.loss<-ifelse(merge.data$chr9p21.3.loss==1,"Loss","No_loss")
merge.data$chr9p21.3.loss<-factor(merge.data$chr9p21.3.loss,levels = c("No_loss","Loss"))
merge.data$arm9p.loss<-ifelse(merge.data$X9p_cnv.loss==1,"Loss","No_loss")
merge.data$arm9p.loss<-factor(merge.data$arm9p.loss,levels = c("No_loss","Loss"))
#plot
set.seed(0)
glist<-list()
for (i in 1:nrow(exp2)){
  #i<-1
  data.use<-merge.data[,c(1,2,ncol(merge.data),i+3)]
  data.use$event<-data.use[,4]
  stat.test1 <- data.use %>%
    wilcox_test(event ~ chr9p21.3.loss) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  stat.test1$y.position<-18
  
  n_fun <- function(x){
    return(data.frame(y = 0.99*20,
                      label = paste0("N=",length(x))))
  }
  g1<-ggplot(data.use, aes(x=chr9p21.3.loss, y=event)) + 
    geom_boxplot(aes(fill = chr9p21.3.loss),alpha = 0.5)+
    geom_point(aes(fill = chr9p21.3.loss, color=chr9p21.3.loss), alpha=0.5,
               size = 2, shape = 21, position = position_jitterdodge()) +
    stat_summary(fun.data = n_fun, geom = "text", 
                 aes(group=chr9p21.3.loss),
                 hjust = 0.5, position = position_dodge(0.8),size=5) +
    scale_fill_manual(values=c("#FF0033","#3300FF")) +
    scale_color_manual(values=c("#FF0033","#3300FF"))+
    theme_classic2()+
    ggtitle(colnames(data.use)[4])+
    ylab("log2 read counts")+xlab("")+scale_y_continuous(limits = c(-5, 20))+
    theme(axis.title.y = element_text( size = 18, face = "plain"),
          axis.text.x =element_text( size = 18, face = "plain"),
          plot.title = element_text(size=18))+
    theme(legend.position = c(0.1, 0.2))+
    stat_pvalue_manual(
      stat.test1,  label = "p",tip.length = 0,size =5
    )
  
  stat.test2 <- data.use %>%
    wilcox_test(event ~ arm9p.loss) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  stat.test2$y.position<-16
  
  g2<-ggplot(data.use, aes(x=arm9p.loss, y=event)) + 
    geom_boxplot(aes(fill = arm9p.loss),alpha = 0.5)+
    geom_point(aes(fill = arm9p.loss, color=arm9p.loss), alpha=0.5,
               size = 2, shape = 21, position = position_jitterdodge()) +
    stat_summary(fun.data = n_fun, geom = "text", 
                 aes(group=arm9p.loss),
                 hjust = 0.5, position = position_dodge(0.8),size=5) +
    scale_fill_manual(values=c("#FF0033","#3300FF")) +
    scale_color_manual(values=c("#FF0033","#3300FF"))+
    theme_classic2()+
    ggtitle(colnames(data.use)[4])+
    ylab("log2 read counts")+xlab("")+scale_y_continuous(limits = c(-5, 20))+
    theme(axis.title.y = element_text( size = 15, face = "plain"),
          axis.text.x =element_text( size = 15, face = "plain"),
          plot.title = element_text(size=15))+
    theme(legend.position = "none")+
    stat_pvalue_manual(
      stat.test2,  label = "p",tip.length = 0,size =5
    )
  glist[[i]]<-g2
}
figure<-ggarrange(glist[[2]],glist[[3]],glist[[5]],glist[[1]],glist[[4]],glist[[6]],ncol = 1,nrow = nrow(exp2))

##### cell lines
datakeep<-read.delim("/Users/zhaox12/Desktop/Teresas_lab/project/0-others/CCLE/Scotte-3.4-others/Teresa_data/final_use/CCLE.HNSC.alldata-Feb27.2021-Xin.txt")
data<-datakeep
cts1 <- read.csv("/Users/zhaox12/Desktop/Teresas_lab/project/0-others/CCLE/Updated_data_2019_version/data/CCLE_RNAseq_reads.csv",sep=",",header = T)
rownames(cts1)<-gsub("[-]",".",cts1$X)
cts1.1<-cts1[rownames(cts1) %in%datakeep$sample,]
cts2<-as.data.frame(t(cts1.1[,-1]))
rownames(cts2)<-make.names(gsub("\\.\\..*","",rownames(cts2)),unique = T)
cts<-cts2[rownames(cts2) %in% gene,]
cts.use<-log2(as.data.frame(t(cts[order(rownames(cts)),order(colnames(cts))]))+1)
cts.use$sample<-rownames(cts.use)
merge.data<-merge(data[,c("sample","chr9p21.3.loss","X9p_cnv.loss")],cts.use,by.y="sample",by.x="sample")
merge.data$chr9p21.3.loss<-ifelse(merge.data$chr9p21.3.loss==1,"Loss","No_loss")
merge.data$chr9p21.3.loss<-factor(merge.data$chr9p21.3.loss,levels = c("No_loss","Loss"))
merge.data$arm9p.loss<-ifelse(merge.data$X9p_cnv.loss==1,"Loss","No_loss")
merge.data$arm9p.loss<-factor(merge.data$arm9p.loss,levels = c("No_loss","Loss"))
#plot
set.seed(0)
glist<-list()
for (i in 1:nrow(cts)){
  #i<-1
  data.use<-merge.data[,c(1,2,ncol(merge.data),i+3)]
  data.use$event<-data.use[,4]
  stat.test3 <- data.use %>%
    wilcox_test(event ~ chr9p21.3.loss) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  stat.test3$y.position<-16
  
  n_fun <- function(x){
    return(data.frame(y = 0.99*20,
                      label = paste0("N=",length(x))))
  }
  g3<-ggplot(data.use, aes(x=chr9p21.3.loss, y=event)) + 
    geom_boxplot(aes(fill = chr9p21.3.loss),alpha = 0.5)+
    geom_point(aes(fill = chr9p21.3.loss, color=chr9p21.3.loss), alpha=0.5,
               size = 2, shape = 21, position = position_jitterdodge()) +
    stat_summary(fun.data = n_fun, geom = "text", 
                 aes(group=chr9p21.3.loss),
                 hjust = 0.5, position = position_dodge(0.8),size=5) +
    scale_fill_manual(values=c("#FF0033","#3300FF")) +
    scale_color_manual(values=c("#FF0033","#3300FF"))+
    theme_classic2()+
    ggtitle(colnames(data.use)[4])+
    ylab("log2 read counts")+xlab("")+scale_y_continuous(limits = c(-5, 20))+
    theme(axis.title.y = element_text( size = 18, face = "plain"),
          axis.text.x =element_text( size = 18, face = "plain"),
          plot.title = element_text(size=18))+
    theme(legend.position = c(0.1, 0.2))+
    stat_pvalue_manual(
      stat.test3,  label = "p",tip.length = 0,size =5
    )
  
  stat.test4 <- data.use %>%
    wilcox_test(event ~ arm9p.loss) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  stat.test4$y.position<-16
  
  g4<-ggplot(data.use, aes(x=arm9p.loss, y=event)) + 
    geom_boxplot(aes(fill = arm9p.loss),alpha = 0.5)+
    geom_point(aes(fill = arm9p.loss, color=arm9p.loss), alpha=0.5,
               size = 2, shape = 21, position = position_jitterdodge()) +
    stat_summary(fun.data = n_fun, geom = "text", 
                 aes(group=arm9p.loss),
                 hjust = 0.5, position = position_dodge(0.8),size=5) +
    scale_fill_manual(values=c("#FF0033","#3300FF")) +
    scale_color_manual(values=c("#FF0033","#3300FF"))+
    theme_classic2()+
    ggtitle(colnames(data.use)[4])+
    ylab("log2 read counts")+xlab("")+scale_y_continuous(limits = c(-5, 20))+
    theme(axis.title.y = element_text( size = 15, face = "plain"),
          axis.text.x =element_text( size = 15, face = "plain"),
          plot.title = element_text(size=15))+
    theme(legend.position = "none")+
    stat_pvalue_manual(
      stat.test4,  label = "p",tip.length = 0,size =5
    )
  glist[[i]]<-g4
}
figure<-ggarrange(glist[[3]],glist[[1]],glist[[5]],glist[[2]],glist[[4]],glist[[6]],ncol = 1,nrow = length(unique(gene)))

### same for "CCL19","CCL21","CXCL9","CXCL10"
}



