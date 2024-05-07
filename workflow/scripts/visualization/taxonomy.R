library(ggplot2)
library(RColorBrewer)
library(viridis)
#library(plotly)
library(patchwork)

#@ write outputs to file as log! 
log_file <- snakemake@log[['log_file']]
log <- file(log_file, open = "wt")
sink(log)
sink(log, type="message")

classification_df <-function(input_file, cut_off){
  
  lines <- readLines(input_file)
  virus_names <- character()
  classification_scores <- numeric()
  
  for (line in lines) {
    parts <- strsplit(line, " ")[[1]]
    virus_name <- paste(parts[-length(parts)], collapse = " ")
    classification_score <- as.numeric(parts[length(parts)])
    
    
    if (!is.na(classification_score) && classification_score != 0 && classification_score >= cut_off) {
      virus_names <- c(virus_names, virus_name)
      classification_scores <- c(classification_scores, classification_score)
    }
    result_df <- data.frame(virus_name = virus_names, classification_score = classification_scores)
  }
  
  return(result_df)
}


plot_classification <- function(df, n, sample_name){
  
  cut_off <- head(df[order(-df$classification_score), ], n)
  cut_off$virus_name <- gsub("\\[|\\]", "", cut_off$virus_name)
  
  custom_colors <- viridis_pal(option = "cividis")(nrow(cut_off))
  
  # Create the plot
  p <- ggplot(cut_off, aes(x = classification_score * 100, y = reorder(virus_name, classification_score))) +
    geom_bar(stat = "identity", fill = custom_colors) +
    labs(title = sample_name,
         x = "Classification Score",
         y = "Classification Content") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8, family = "Times", face = "bold"),
          panel.grid = element_blank())
  
  
  
  return(p)
}

taxonomy_path <- snakemake@input[['input_taxonomy']]
sample_name <-  snakemake@params[['sample_name']]
output_path <- snakemake@output[['output_vis']]
depth <- snakemake@params[['depth']]

cap = paste0("Normalized Assignment of ", sample_name)
df = classification_df(taxonomy_path, 0.000000001)
p = plot_classification(df, depth, cap)

ggsave(output_path, plot = p, width = 7, height = 9, units = "in", dpi = 300)
