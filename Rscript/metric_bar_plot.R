library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

df_blank <- data.frame(c("blank1","blank2"), c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0))

top_percent <- c(0.025, 0.05, 0.1)
dir <- paste0("top_union_",top_percent[1],"_t1_comparison_test.csv")
df <- read.csv(dir, header=T)
colnames(df_blank) <- colnames(df)
df <- rbind(df, df_blank)
df$top_percent <- top_percent[1]
for (i in 2:length(top_percent)){
  dir <- paste0("top_union_",top_percent[i],"_t1_comparison_test.csv")
  temp <- read.csv(dir, header=T)
  temp <- rbind(temp, df_blank)
  temp$top_percent <- top_percent[i]
  df <- rbind(df, temp)
}

df$top_percent <- factor(df$top_percent)

SNP_set <- c("P-value","blank1","blank2","Recall","Accuracy","Precision","Balanced Accuracy","ROC AUC","F1 Score","Adaptive LASSO")
df$SNP.Set <- factor(df$SNP.Set, levels=SNP_set)

#color_set <- c("#d95f02","red", "darkred","#e5f5f9","#deebf7","#c6dbef","#d0d1e6","#9ecae1","#a6bddb","#1b9e77")
#color_set_legend <- c("#d95f02","#e5f5f9","#deebf7","#c6dbef","#d0d1e6","#9ecae1","#a6bddb","#1b9e77")

color_set <- c("P-value"="#d95f02","blank1"="red","blank2"="darkred","Recall"="#e5f5f9","Accuracy"="#deebf7","Precision"="#c6dbef","Balanced Accuracy"="#d0d1e6","ROC AUC"="#9ecae1","F1 Score"="#a6bddb","Adaptive LASSO"="#1b9e77")


metric_name <- c("Balanced Accuracy","F1 Score","Precision","P-value","Recall","ROC AUC")
metric_vec <- c("Balanced.Accuracy","F1.Score","Precision","P.value","Recall","ROC.AUC")

colnames(df) <- c("SNP_set", metric_name, "top_percent")

df_long <- df %>% pivot_longer(cols = -c("SNP_set", "top_percent"), names_to = "variable", values_to = "value")

p_list <- vector(mode = "list", length = length(metric_name))

i=1 # deal with balanced accuracy, cut y-axis at 0.5
df_temp <- df_long[df_long$variable==metric_name[i],]
p_list[[i]] <- ggplot(df_temp, aes(x = top_percent, y = value, fill = SNP_set)) +
  geom_bar(stat = "identity", position = position_dodge2(padding=0.1,width=0.8), width=0.6) +
  scale_x_discrete(labels=c("0.025"="0.025%", "0.05"="0.05%", "0.1"="0.1%")) + 
  facet_wrap(.~ variable, scales="free_y", ncol=1) + 
  coord_cartesian(ylim=c(0.5, NA)) + 
  scale_fill_manual(values = color_set) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "none"
  )

for (i in 2:(length(metric_name)-1)){
  df_temp <- df_long[df_long$variable==metric_name[i],]
  print(sort(df_temp$value))
  p_list[[i]] <- ggplot(df_temp, aes(x = top_percent, y = value, fill = SNP_set)) +
    geom_bar(stat = "identity", position = position_dodge2(padding=0.1,width=0.8), width=0.6) +
    scale_x_discrete(labels=c("0.025"="0.025%", "0.05"="0.05%", "0.1"="0.1%")) + 
    facet_wrap(.~ variable, scales="free_y", ncol=1) + 
    coord_cartesian(ylim=c(0.9*sort(df_temp$value)[7], NA)) + 
    scale_fill_manual(values = color_set) +
    theme(axis.line = element_line(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "none"
    )
}
# deal with roc auc cut y-axis at 0.5
df_temp <- df_long[df_long$variable==metric_name[length(metric_name)],]
p_list[[length(metric_name)]] <- ggplot(df_temp, aes(x = top_percent, y = value, fill = SNP_set)) +
  geom_bar(stat = "identity", position = position_dodge2(padding=0.1,width=0.8), width=0.6) +
  scale_x_discrete(labels=c("0.025"="0.025%", "0.05"="0.05%", "0.1"="0.1%")) + 
  facet_wrap(.~ variable, scales="free_y", ncol=1) + 
  coord_cartesian(ylim=c(0.5, NA)) + 
  labs(x="Top Percentage") + 
  scale_fill_manual(values = color_set,
                    breaks = c("P-value","Recall","Accuracy","Precision","Balanced Accuracy","ROC AUC","F1 Score","Adaptive LASSO")) + 
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom"
  )
p <- p_list[[1]]
for (i in 2:length(p_list)){
  p <- p + p_list[[i]]
}
p <- p + plot_layout(ncol = 1)




df_blank <- data.frame(c("blank1","blank2"), c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0))

top_percent <- c(0.025, 0.05, 0.1)
dir <- paste0("top_union_",top_percent[1],"_comparison_ADSP.csv")
df <- read.csv(dir, header=T)
colnames(df_blank) <- colnames(df)
df <- rbind(df, df_blank)
df$top_percent <- top_percent[1]
for (i in 2:length(top_percent)){
  dir <- paste0("top_union_",top_percent[i],"_comparison_ADSP.csv")
  temp <- read.csv(dir, header=T)
  temp <- rbind(temp, df_blank)
  temp$top_percent <- top_percent[i]
  df <- rbind(df, temp)
}

df$top_percent <- factor(df$top_percent)

SNP_set <- c("P-value","blank1","blank2","Recall","Accuracy","Precision","Balanced Accuracy","ROC AUC","F1 Score","Adaptive LASSO")
df$SNP.Set <- factor(df$SNP.Set, levels=SNP_set)

#color_set <- c("#d95f02","red", "darkred","#e5f5f9","#deebf7","#c6dbef","#d0d1e6","#9ecae1","#a6bddb","#1b9e77")
#color_set_legend <- c("#d95f02","#e5f5f9","#deebf7","#c6dbef","#d0d1e6","#9ecae1","#a6bddb","#1b9e77")

color_set <- c("P-value"="#d95f02","blank1"="red","blank2"="darkred","Recall"="#e5f5f9","Accuracy"="#deebf7","Precision"="#c6dbef","Balanced Accuracy"="#d0d1e6","ROC AUC"="#9ecae1","F1 Score"="#a6bddb","Adaptive LASSO"="#1b9e77")


metric_name <- c("Balanced Accuracy","F1 Score","Precision","P-value","Recall","ROC AUC")
metric_vec <- c("Balanced.Accuracy","F1.Score","Precision","P.value","Recall","ROC.AUC")

colnames(df) <- c("SNP_set", metric_name, "top_percent")

df_long <- df %>% pivot_longer(cols = -c("SNP_set", "top_percent"), names_to = "variable", values_to = "value")

p_list <- vector(mode = "list", length = length(metric_name))

i=1 # deal with balanced accuracy, cut y-axis at 0.5
df_temp <- df_long[df_long$variable==metric_name[i],]
p_list[[i]] <- ggplot(df_temp, aes(x = top_percent, y = value, fill = SNP_set)) +
  geom_bar(stat = "identity", position = position_dodge2(padding=0.1,width=0.8), width=0.6) +
  scale_x_discrete(labels=c("0.025"="0.025%", "0.05"="0.05%", "0.1"="0.1%")) + 
  facet_wrap(.~ variable, scales="free_y", ncol=1) + 
  coord_cartesian(ylim=c(0.5, NA)) + 
  scale_fill_manual(values = color_set) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "none"
  )

for (i in 2:(length(metric_name)-1)){
  df_temp <- df_long[df_long$variable==metric_name[i],]
  print(sort(df_temp$value))
  p_list[[i]] <- ggplot(df_temp, aes(x = top_percent, y = value, fill = SNP_set)) +
    geom_bar(stat = "identity", position = position_dodge2(padding=0.1,width=0.8), width=0.6) +
    scale_x_discrete(labels=c("0.025"="0.025%", "0.05"="0.05%", "0.1"="0.1%")) + 
    facet_wrap(.~ variable, scales="free_y", ncol=1) + 
    coord_cartesian(ylim=c(0.9*sort(df_temp$value)[7], NA)) + 
    scale_fill_manual(values = color_set) +
    theme(axis.line = element_line(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "none"
    )
}
# deal with roc auc cut y-axis at 0.5
df_temp <- df_long[df_long$variable==metric_name[length(metric_name)],]
p_list[[length(metric_name)]] <- ggplot(df_temp, aes(x = top_percent, y = value, fill = SNP_set)) +
  geom_bar(stat = "identity", position = position_dodge2(padding=0.1,width=0.8), width=0.6) +
  scale_x_discrete(labels=c("0.025"="0.025%", "0.05"="0.05%", "0.1"="0.1%")) + 
  facet_wrap(.~ variable, scales="free_y", ncol=1) + 
  coord_cartesian(ylim=c(0.5, NA)) + 
  labs(x="Top Percentage") + 
  scale_fill_manual(values = color_set,
                    breaks = c("P-value","Recall","Accuracy","Precision","Balanced Accuracy","ROC AUC","F1 Score","Adaptive LASSO")) + 
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom"
  )
p <- p_list[[1]]
for (i in 2:length(p_list)){
  p <- p + p_list[[i]]
}
p <- p + plot_layout(ncol = 1)

















