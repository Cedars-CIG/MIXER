library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(patchwork)

args <- commandArgs(TRUE)


##############################
## get the top percent plot
##############################
metric_order <- c("P-value", "Accuracy", "AUC", "Balanced Accuracy", "F1 Score", "Precision", "Recall")

metric <- c("accuracy", "balanced_accuracy", "f1_score", "precision", "pval", "recall", "roc_auc")
metric_plot <- c("Accuracy", "Balanced Accuracy", "F1 Score", "Precision", "P-value", "Recall", "AUC")


df_list <- vector("list", length(metric))

for (i in 1:length(metric)){
    df <- as.data.frame(fread(paste0(args[1], metric[i], "_dist.csv"), header=T))
    df$metric <- metric_plot[i]
    print(head(df))
    df_list[[i]] <- df
}

df_plot <- bind_rows(df_list)
df_plot$plot_label <- "Selected SNP P-values"
df_plot$metric <- factor(df_plot$metric, levels = metric_order)

print(head(df_plot))
print(tail(df_plot))

p1 <- ggplot(df_plot, aes(x=pval)) + 
    geom_histogram(aes(y=after_stat(count/sum(count))), bins=100, fill = "#FFD877") + 
    ggh4x::facet_grid2(metric ~ plot_label, scales="free", independent = c("all")) + 
    scale_x_continuous(limits=c(0,1)) +
    labs(x = "P-value") + 
    theme(axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.title.y = element_blank(),
            strip.text.y = element_blank(),
            strip.text.x = element_text(size=15),
            ggh4x.axis.bottom = element_line(),  # ensures x-axis line shows even if shared
        )


##############################
## get the MAF plot
##############################

outlier_IQR <- function(df){

    # Calculate Q1, Q3 and IQR
    Q1 <- quantile(df[,2], 0.25)
    Q3 <- quantile(df[,2], 0.75)
    IQR <- Q3 - Q1

    # Define lower and upper bounds for non-outliers
    lower_bound <- Q1 - 1.5 * IQR
    upper_bound <- Q3 + 1.5 * IQR

    # Filter the dataframe to remove outliers
    df_clean <- df[df[,2] >= lower_bound & df[,2] <= upper_bound, ]
    return(df_clean)

}


df_list <- vector("list", length(metric))

for (i in 1:length(metric)){
    df <- as.data.frame(fread(paste0(args[1], metric[i], "_maf.csv"), header=T))
    df_list[[i]] <- outlier_IQR(df)
}

df_plot <- bind_rows(df_list)
df_plot$plot_label <- "MAF Distribution"
df_plot$metric <- factor(df_plot$metric, levels = metric_order)

p2 <- ggplot(df_plot, aes(x=MAF)) + 
    geom_histogram(aes(y=after_stat(count/sum(count))), bins = 50, fill = "#4FD3D9") +  
    ggh4x::facet_grid2(metric ~ plot_label, scales="free", independent = c("all")) + 
    scale_x_continuous(limits=c(0, 0.5)) + 
    labs(x = "Minor Allele Frequency") + 
    theme(axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.title.y = element_blank(),
            strip.text.y = element_blank(),
            strip.text.x = element_text(size=15),
            ggh4x.axis.bottom = element_line(),  # ensures x-axis line shows even if shared
        )


##############################
## get the betas plot
##############################


df_list <- vector("list", length(metric))

for (i in 1:length(metric)){
    df <- as.data.frame(fread(paste0(args[1], metric[i], "_gwas.csv"), header=T))
    df_list[[i]] <- outlier_IQR(df)
}

df_plot <- bind_rows(df_list)
df_plot$plot_label <- "Effect Size Distribution"
df_plot$metric <- factor(df_plot$metric, levels = metric_order)

p3 <- ggplot(df_plot, aes(x=beta)) + 
    geom_histogram(aes(y=after_stat(count/sum(count))), bins = 50, fill = "#FF9A80") + 
    ggh4x::facet_grid2(metric ~ plot_label, scales="free", independent = c("all")) + 
    scale_x_continuous(limits=c(min(df_plot$beta), max(df_plot$beta))) + 
    labs(x = "Effect Size") + 
    theme(axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.title.y = element_blank(),
            strip.text.x = element_text(size=15),
            ggh4x.axis.bottom = element_line(),  # ensures x-axis line shows even if shared
        )

p <- p1 + p2 + p3 + plot_layout(ncol=3)
ggsave(paste0(args[1], "top_SNP_dist.jpeg"), p, device="jpeg", width=12, height=11)
