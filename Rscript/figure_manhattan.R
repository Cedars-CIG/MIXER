library(dplyr)
library(readr)
library(data.table)
library(reshape2)
library(ggplot2)
library(latex2exp)
library(scales)
library(ggh4x)

args <- commandArgs(TRUE)


if (args[2] == "_finngen_data.csv"){
  pheno_list <- c("AD", "breast_cancer", "hypertension", "t2d", "bmi")
  plot_title_list <- c("AD", "Breast Cancer", "Hypertension", "T2D", "BMI")
} else {
  pheno_list <- c("AD", "breast_cancer", "hypertension", "t2d", "bmi", "cholesterol", "hdl", "ldl")
  plot_title_list <- c("AD", "Breast Cancer", "Hypertension", "T2D", "BMI", "Cholesterol", "HDL", "LDL")
}

df_bp_info <- read.csv("/common/lir2lab/Wanling/FDR_AB_PRS/Analysis_figure/manhattan/bp_info.csv", header=T)

df_long <- read.csv(paste0(args[1],pheno_list[1],args[2]),header=T)
df_long$pheno <- plot_title_list[1]

for (i in 2:length(pheno_list)) {
  temp <- read.csv(paste0(args[1],pheno_list[i],args[2]),header=T)
  temp$pheno <- plot_title_list[i]
  df_long <- rbind(df_long, temp)
}
df_long$pheno <- factor(df_long$pheno, levels=plot_title_list)

# Define shapes and colors
shapes <- c("Additive" = 15, "Dominant" = 16, "Recessive" = 17, "Co-dominant" = 18)
#colors <- c("Additive" = "#CC79A7", "Dominant" = "#0072B2", "Recessive" = "#009E73", "Co-dominant" = "#E69F00", "Sub-optimal" = "grey")
#colors <- c("Additive" = "#7FB4A6", "Dominant" = "#E3D5AE", "Recessive" = "#B5823E", "Co-dominant" = "#82531E", "Sub-optimal" = "grey")
colors <- c("Additive" = "#7FB4A6", "Dominant" = "#E3D5AE", "Recessive" = "#82531E", "Co-dominant" = "#B5823E", "Sub-optimal" = "grey")


# Custom transformation function
custom_trans <- function(y0) {
  trans_new(
    name = "custom",
    transform = function(y) ifelse(y <= y0, y, y0 + (y - y0) / 50),
    inverse = function(y) ifelse(y <= y0, y, y0 + (y - y0) * 50)
  )
}


# break point
y0 <- 20

BP_start <- df_bp_info$BP_start[-1]
df_long$optimal_type <- factor(df_long$optimal_type, levels =c("Additive", "Dominant", "Recessive", "Co-dominant", "Sub-optimal"))
## reduce grey points

df_add <- df_long[df_long$type == "Additive",]

df_non_add <- df_long[df_long$type != "Additive",]
df_non_add <- df_non_add[df_non_add$optimal_type != "Sub-optimal",]

df_adjust <- rbind(df_add, df_non_add)
df_adjust$pheno <- factor(df_adjust$pheno, levels=plot_title_list)
df_adjust$optimal_type <- factor(df_adjust$optimal_type, levels =c("Additive", "Dominant", "Recessive", "Co-dominant", "Sub-optimal"))

# Plot the data
p <- ggplot(df_adjust, aes(x = BPcum, y = pval, shape = type, color = optimal_type, alpha = is_max)) +
  geom_point(size = 3) +
  facet_wrap(~pheno, nrow = 8, scales = "free", strip.position="right") + 
  scale_shape_manual(values = shapes) +
  scale_color_manual(labels = c("Additive", "Dominant", "Recessive", "Co-dominant",  "Sub-optimal Additive"),
                     values = colors) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.8)) +  # Set alpha values based on is_max
  labs(x = "Chromosome", y = TeX("$-log_{10}(p-value)$"), color = "") +
  guides(shape = "none",  # Remove shape legend
         color = guide_legend(override.aes = list(shape = c(15, 16, 17, 18, 15))),
         alpha = "none") +  # Remove alpha from legend
  scale_x_continuous(label = df_bp_info$CHR, breaks= df_bp_info$BP_mid) +
  scale_y_continuous(trans = custom_trans(y0),
                     breaks = function(y) {
                      if (max(y) > 100) return(c(seq(0, y0-5, by = 5), seq(y0+80, max(y), by = 200))) else {return(c(seq(0, y0-5, by = 5), 70))}
                     }) + 
  geom_vline(xintercept=c(df_bp_info$BP_start, df_bp_info$BP_end[22]), linetype="dashed", colour="grey") + 
  annotate("text", y = c(y0-0.2, y0+20), x = -Inf, label = "~", size=unit(6,"pt")) + coord_cartesian(clip = "off") + 
  #geom_rect(aes(xmin=-Inf, xmax=+Inf, ymin=y0-0.4, ymax=y0+2, fill="white", color="white")) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "top",
        strip.text = element_text(angle=0, hjust=0.5),
        panel.spacing = unit(1, "lines"),
        axis.text.y = element_text()
        ) + 
  guides(y = guide_axis_truncated(
          trunc_lower = c(-Inf, y0+20),
          trunc_upper = c(y0-0.2, Inf)
        ))

if (args[2] == "_finngen_data.csv"){
  ggsave(
        args[3],
        plot = p,
        width = 8.27,
        height = 7.31,
        dpi = 240,
        units = "in"
    )
} else {
  ggsave(
        args[3],
        plot = p,
        width = 8.27,
        height = 11.69,
        dpi = 240,
        units = "in"
    )
}