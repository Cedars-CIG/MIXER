library(dplyr)
library(ggplot2)

args <- commandArgs(TRUE)

top_vec <- c("0.025", "0.05", "0.1")

df_list <- vector("list", 3)
for (i in 1:3){
     df_list[[i]] <- read.table(paste0(args[1], top_vec[i], "_OR_data_test"), header=T)
     df_list[[i]]$dataset <- "UKBB Testing"
     df_list[[i]]$top <- paste0("Top ", top_vec[i], "%")
}
df_UKBB <- bind_rows(df_list)

for (i in 1:3){
     df_list[[i]] <- read.table(paste0(args[1], top_vec[i], "_OR_data_ADSP"), header=T)
     df_list[[i]]$dataset <- "ADSP"
     df_list[[i]]$top <- paste0("Top ", top_vec[i], "%")
}
df_ADSP <- bind_rows(df_list)

df <- rbind(df_UKBB, df_ADSP)

df <- df %>% mutate(SNP_set = recode(SNP_set, 
                                     "accuracy" = "Accuracy",
                                     "balanced_accuracy" = "Balanced Accuracy",
                                     "f1_score" = "F1 Score",
                                     "precision" = "Precision",
                                     "pval" = "P-value",
                                     "recall" = "Recall",
                                     "roc_auc" = "AUC",
                                     "MIXER" = "MIXER"))


SNP_set <- c("Accuracy","AUC","Balanced Accuracy","F1 Score","P-value","Precision","Recall","MIXER")
df$SNP_set <- factor(df$SNP_set, levels=SNP_set)

df$dataset <- factor(df$dataset, levels=c("UKBB Testing", "ADSP"))

color_set <- c("P-value"="#f57c6e","Recall"="#fbe0be","Accuracy"="#e3ede0",
               "Precision"="#d3eee2","Balanced Accuracy"="#d0e0ef","AUC"="#f1d0c6",
               "F1 Score"="#e4e3bf","MIXER"="#6bca73")
p <- ggplot(data=df,aes(x=percent,y=OR,color=SNP_set)) + 
        geom_point(size = 3) +
        geom_line() + 
        geom_errorbar(aes(ymin = lw, ymax = up), width = 0.4) +  # Add error bars
        ggh4x::facet_grid2(dataset ~ top, scales="free", independent = "all") + 
        labs(x = "Percentile Threshold of the High-Risk Group", 
             y = "Odds Ratio",
             color = "Selection Method") +
        #scale_x_continuous(limits=c(0.04, 0.96), breaks=seq(0.05, 0.95, by=0.1)) +
        scale_x_continuous(limits=c(0.04, 0.21), breaks=seq(0.05, 0.20, by=0.05)) +
        scale_color_manual("SNP Set", values=color_set) + 
        guides(color = guide_legend(nrow=1)) + 
        theme(axis.line = element_line(),
              panel.border = element_blank(), 
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              legend.position = "bottom",
              legend.box = "horizontal")

ggsave(args[2], p, device="jpeg", width=10, height=6)