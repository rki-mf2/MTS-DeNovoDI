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

hist_df <-function(input_file){
  
  lines <- readLines(input_file)
  virus_names <- character()
  lengths <- numeric()
  
  for (line in lines) {
    parts <- strsplit(line, " ")[[1]]
    virus_name <- paste(parts[-length(parts)], collapse = " ")
    length <- as.numeric(parts[length(parts)]) 
  
  virus_names <- c(virus_names, virus_name)
  lengths <- c(lengths, length)
    
    result_df <- data.frame(virus_name = virus_names, length = lengths)
  }
  
  return(result_df)
}


plot_hist <- function(df, n, db_name){
  
  cut_off <- head(df[order(-df$length), ], n)
  cut_off$virus_name <- gsub("\\[|\\]", "", cut_off$virus_name)
  
  custom_colors <- viridis_pal(option = "plasma")(nrow(cut_off))
  
  p <- ggplot(cut_off, aes(x = length, y = reorder(virus_name, length))) +
    geom_bar(stat = "identity", fill = custom_colors) +
    labs(title = db_name,
         x = "Length (AA)",
         y = "Sequence Record") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8, family = "Times", face = "bold"),
          panel.grid = element_blank())
  
  
  
  return(p)
}


input_path <- snakemake@input[['input_hist']]
df = hist_df(input_path)
depth <- snakemake@params[['depth']]
p = plot_hist(df, depth, "Database AA Histogram")
output_path <- snakemake@output[['output_hist']]
ggsave(output_path, plot = p, width = 7, height = 9, units = "in", dpi = 300)
