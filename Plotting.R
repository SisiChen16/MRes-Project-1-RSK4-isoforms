# Load packages 
library("gridExtra")
library("cowplot")
library(ggplot)


# Import data 
setwd("/Users/chensisi/Documents/RNAseq/")
load("TCGA_GTEX_merge.RData")
load("RSK4_data.RData")

# set the theme for ggplot 
my_theme <- theme(axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14),
                  axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 14),
                  plot.title = element_text(size = 16), legend.text = element_text(size = 12), legend.title = element_text(size =12)) + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black"),
                  panel.border = element_blank())



## Proportion of 4 transcript isoforms of RSK4 in TCGA and GTEX samples
TCGA_mean <- merge1_denorm %>%
  filter(Project == "TCGA") %>%
  select("ENST00000620340.4","ENST00000495332.1","ENST00000460730.1","ENST00000262752.4","Project","Sample_type") %>%
  group_by(Sample_type) %>%
  summarize(iso1= mean(ENST00000620340.4), iso2 = mean(ENST00000262752.4),
            iso3= mean(ENST00000495332.1), iso4 = mean(ENST00000460730.1)) %>%
  mutate(total= iso1 + iso2 + iso3 + iso4) 

GTEX_mean <- merge1_denorm %>%
  filter(Project == "GTEX") %>%
  select("ENST00000620340.4","ENST00000495332.1","ENST00000460730.1","ENST00000262752.4","Project","Sample_type") %>%
  group_by(Sample_type) %>%
  summarize(iso1= mean(ENST00000620340.4), iso2 = mean(ENST00000262752.4),
            iso3= mean(ENST00000495332.1), iso4 = mean(ENST00000460730.1)) %>%
  mutate(total= iso1 + iso2 + iso3 + iso4)
 
# Barplot showing total expression of RSK4 in each cancer type 
p1 <- ggplot(TCGA_mean, aes(x= Sample_type, y= total)) + geom_bar(stat = "identity", width = 0.9) + 
  theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
         axis.ticks.x=element_blank()) + 
  ylab("Total expression")

p1 <- ggplot(GTEX_mean, aes(x= Sample_type, y= total)) + geom_bar(stat = "identity", width = 0.9) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  ylab("Total expression")


# Stacked barplot showing the proportion of RSK4 isoforms in each cancer type
p2 <- TCGA_mean %>%
  gather(Isoform, expression, iso1:iso4) %>%
  ggplot(aes(x = Sample_type, y = expression, fill = isoform)) + 
  geom_bar(position = "fill", stat = "identity") + 
  scale_y_continuous(labels = scales::percent_format()) + 
  labs(y = "Proportion", x = "Cancer types") + 
  theme(axis.text.x = element_text(angle = 90,  vjust = 0.5, hjust=1), 
        plot.title = element_text(hjust = 0.5), 
        legend.position = "none")  # remove legends 

p2 <- GTEX_mean %>%
  gather(isoform, expression, iso1:iso4) %>%
  ggplot(aes(x = Sample_type, y = expression, fill = isoform)) + 
  geom_bar(position = "fill", stat = "identity") + 
  scale_y_continuous(labels = scales::percent_format()) + 
  labs(y = "Proportion", x = "Cancer types") + 
  theme(axis.text.x = element_text(angle = 90,  vjust = 0.5, hjust=1), 
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")  # remove legends 

png("GTEXstacked_plot.png",width=10,height=7,units='in',res=500)
plot_grid(p1, p2, labels=c("A", "B"), ncol = 1 , align = "v", axis = "bt",rel_widths = c(1,1.3), rel_heights = c(1.4,4))
dev.off()



## Boxplot showing the distribution of RSK4 isoform expression among normal and tumor samples
  # Categorize GTEX and TCGA samples into high, medium and low isoform expression groups ####
    #TCGA 
TCGA_by_sample_median <- TCGA_by_sample %>% 
  group_by(Sample_type) %>%
  summarize(median1 = median(ENST00000620340.4), median2 = median(ENST00000262752.4))

TCGA_by_sample_median$Expression <- NA

iso1lower_bound <- quantile(TCGA_by_sample$ENST00000620340.4, 0.25)
iso2lower_bound <- quantile(TCGA_by_sample$ENST00000262752.4, 0.25)
iso1upper_bound <- quantile(TCGA_by_sample$ENST00000620340.4, 0.75)
iso2upper_bound <- quantile(TCGA_by_sample$ENST00000262752.4, 0.75)

#High 
iso1.high <- which(TCGA_by_sample_median$median1 > iso1upper_bound)
iso2.high <- which(TCGA_by_sample_median$median2 > iso2upper_bound)

#Medium 
iso1.medium <- which(TCGA_by_sample_median$median1 > iso1lower_bound & TCGA_by_sample_median$median1 < iso1upper_bound)
iso2.medium <- which(TCGA_by_sample_median$median2 > iso2lower_bound & TCGA_by_sample_median$median2 < iso2upper_bound)

#Low 
iso1.low <- which(TCGA_by_sample_median$median1 < iso1lower_bound | TCGA_by_sample_median$median1 == iso1lower_bound)
iso2.low <- which(TCGA_by_sample_median$median2 < iso2lower_bound | TCGA_by_sample_median$median2 == iso2lower_bound)

#Plot the figure 
TCGA_by_sample_median$Expression[iso1.high] <- "High"
TCGA_by_sample_median$Expression[iso1.medium] <- "Medium"
TCGA_by_sample_median$Expression[iso1.low] <- "Low"
TCGA_by_sample_median <- TCGA_by_sample_median %>%
  arrange(desc(median1))

ggplot(TCGA_by_sample_median, aes(x= ENST00000620340.4, y = reorder(Sample_type,median1), color= Expression)) + 
  geom_boxplot() + 
  labs( x = "Expression (log2 normalized TPM)", y = "Cancer types",title = "RSK4 isoform1 expression in TCGA samples") + 
  theme(axis.text.y = element_text(size = 12), axis.title = element_text(size = 14), 
        plot.title = element_text( face = "bold", hjust = 0.5)) + 
  geom_vline(xintercept= c(iso1upper_bound, iso1lower_bound) , linetype="dashed", color = "red")


 
    #GTEX 
GTEX_by_sample_median <- GTEX_by_sample %>% 
  group_by(Sample_type) %>%
  mutate(median1 = median(ENST00000620340.4), median2 = median(ENST00000262752.4)) 

GTEX_by_sample_median$Expression <- NA

iso1lower_bound <- quantile(GTEX_by_sample$ENST00000620340.4, 0.25)
iso2lower_bound <- quantile(GTEX_by_sample$ENST00000262752.4, 0.25)
iso1upper_bound <- quantile(GTEX_by_sample$ENST00000620340.4, 0.75)
iso2upper_bound <- quantile(GTEX_by_sample$ENST00000262752.4, 0.75)

# Isoform 1 
high <- which(GTEX_by_sample_median$median1 >= iso1upper_bound)
medium <- which(GTEX_by_sample_median$median1 < iso1upper_bound & GTEX_by_sample_median$median1 >= iso1lower_bound)
low <- which(GTEX_by_sample_median$median1 < iso1lower_bound )

# Isoform 2 
high <- which(GTEX_by_sample_median$median2 >= iso2upper_bound)
medium <- which(GTEX_by_sample_median$median2 < iso2upper_bound & GTEX_by_sample_median$median2 >= iso2lower_bound)
low <- which(GTEX_by_sample_median$median2 < iso2lower_bound )

GTEX_by_sample_median$Expression[high] <- "High"
GTEX_by_sample_median$Expression[medium] <- "Medium"
GTEX_by_sample_median$Expression[low] <- "Low"
GTEX_by_sample_median <- GTEX_by_sample_median %>%
  arrange(desc(median1))

# Plot the figure  
ggplot(GTEX_by_sample_median, aes(x= ENST00000620340.4, y = reorder(Sample_type,median1), color= Expression)) + 
  geom_boxplot() + 
  labs( x = "Expression (log2 normalized TPM)", y = "Normal tissues",title = "RSK4 isoform1 expression in GTEX samples") + 
  theme(axis.text.y = element_text(size = 12), axis.title = element_text(size = 14), 
        plot.title = element_text( face = "bold", hjust = 0.5)) + 
  geom_vline(xintercept= c(iso1upper_bound, iso1upper_bound) , linetype="dashed", color = "red")











