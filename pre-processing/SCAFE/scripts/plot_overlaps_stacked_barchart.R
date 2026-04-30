library(optparse)
library(tidyverse)
library(ggplot2)

option_list = list(make_option(c("--output"), type = "character", default = NULL, help = "Output file", metavar = "character"),
                   make_option(c("--input"), type = "character", default = NULL, help = "CRE_set", metavar = "character")
)

opt = parse_args(OptionParser(option_list = option_list))
output <- opt$output
input_path <- opt$input

df <- read_csv(input_path)
names(df) <- c("name", "None", "One", "Two", "Three")

df_long <- df %>% 
  pivot_longer(cols = None:Three, names_to = "Overlap", values_to = "Count") %>%
  mutate(Overlap = factor(Overlap, levels = c("None", "One", "Two", "Three"))) %>%  # Explicitly set factor levels
  group_by(name) %>%  # Group by the CRE set name
  mutate(Proportion = Count / sum(Count))  # Calculate proportions

pdf(output)
ggplot(df_long, aes(x = name, y = Proportion, fill = Overlap)) +
  geom_bar(stat = "identity") +
  labs(x = "CRE Set Name", y = "Proportion of Overlaps", title = "Normalized Stacked Bar Chart of CRE Set Overlaps") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1")
dev.off()
