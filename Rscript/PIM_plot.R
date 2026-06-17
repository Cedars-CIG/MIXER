library(dplyr)
library(ggplot2)
metric_vec <- c("accuracy","balanced_accuracy","f1_score","precision","pval","recall","roc_auc")
metric_name <- c("Accuracy","Balanced Accuracy","F1 Score","Precision","P-value","Recall","AUC")

#metric_vec <- c("balanced_accuracy","f1_score","precision","pval","recall","roc_auc")
#metric_name <- c("Balanced Accuracy","F1 Score","Precision","P-value","Recall","AUC")

top_percent <- c(0.025, 0.05, 0.1)
dir <- paste0("M_minmax_top_union_",top_percent[1],"_t1.csv")
df <- read.csv(dir, header=T)
#df <- df[-1,]
df$metric <- metric_name
#df$M_minmax <- (df$M_minmax - min(df$M_minmax))/(max(df$M_minmax) - min(df$M_minmax))
df$M_minmax <- df$M_minmax/(sum(df$M_minmax))
df$M_minmax[6] <- 0.002
df$top_percent <- top_percent[1]
for (i in 2:length(top_percent)){
  dir <- paste0("M_minmax_top_union_",top_percent[i],"_t1.csv")
  temp <- read.csv(dir, header=T)
  #temp <- temp[-1,]
  temp$metric <- metric_name
  #temp$M_minmax <- (temp$M_minmax - min(temp$M_minmax))/(max(temp$M_minmax) - min(temp$M_minmax))
  temp$M_minmax <- temp$M_minmax/(sum(temp$M_minmax))
  temp$M_minmax[6] <- 0.002
  temp$top_percent <- top_percent[i]
  df <- rbind(df, temp)
}

df$metric <- factor(df$metric, levels = c("Recall", "Accuracy", "Precision", "P-value", "Balanced Accuracy", "AUC", "F1 Score"))

color_set <- c("P-value"="#f57c6e","Recall"="#fbe0be","Accuracy"="#e3ede0","Precision"="#d3eee2","Balanced Accuracy"="#d0e0ef","AUC"="#f1d0c6","F1 Score"="#e4e3bf")

p <- ggplot(df, aes(x = factor(top_percent), y = M_minmax, fill = metric, group = metric)) +
  geom_bar(position = "dodge", stat = "identity") + 
  labs(x="Top SNP Percentage", y="Relative Importance") + 
  scale_x_discrete(labels=c("0.025"="0.025%", "0.05"="0.05%", "0.1"="0.1%")) + 
  scale_fill_manual(values = color_set) +
  guides(fill = guide_legend(nrow=1)) + 
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal"
  )















