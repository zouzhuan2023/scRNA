###Load the package###
library(ggplot2)
###Set the working path to###
setwd("/path")
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  cache = TRUE,
  out.width = "100%"
)
coords <- read.csv("coords.csv")
Fib_ratio <- read.csv("Fibroblast_data_with_rownames.csv")
counts <- read.csv("CDKN1A_matrix_export.csv")
combine <- merge(coords,y= c(Fib_ratio,counts))
combine <- combine[,-6]
combine <- combine[,-1]
df <- combine
median_CDKN1A <- median(df$CDKN1A.denoised.)
df$color_label <- with(df, ifelse(CDKN1A.denoised. > median_CDKN1A & Fibroblast > 0.1, "CDKN1A High & Fibroblast High",
                                  ifelse(CDKN1A.denoised. < median_CDKN1A & Fibroblast > 0.1, "CDKN1A Low & Fibroblast High",
                                         ifelse(CDKN1A.denoised. > median_CDKN1A & Fibroblast < 0.1, "CDKN1A High & Fibroblast Low",
                                                "CDKN1A Low & Fibroblast Low"))))

print(table(df$color_label))
###Define color mapping###
color_mapping <- c("CDKN1A High & Fibroblast High" = "red",
                   "CDKN1A Low & Fibroblast High" = "lightcoral",
                   "CDKN1A High & Fibroblast Low" = "blue",
                   "CDKN1A Low & Fibroblast Low" = "lightblue")
p <- ggplot(df, aes(x = x, y = y)) +
  geom_point(aes(color = color_label), size = 2) +
  scale_color_manual(values = color_mapping) +
  coord_fixed() +  
  theme_minimal() +
  theme(panel.grid = element_blank(),  
        axis.line = element_blank(),  
        axis.text = element_blank(),  
        axis.title = element_blank()) +  
  labs(title = "Scatter Plot with Conditional Coloring",
       color = "Condition")  
print(p)
ggsave("Fibroblast_0.1_CDKN1A.pdf",p,width = 10,height = 6)
