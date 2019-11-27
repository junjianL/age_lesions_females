
library(ggplot2)
reads_c <- readr::read_tsv("data/bismark_alignment.tsv")
head(reads_c)
colnames(reads_c)

ggplot(reads_c) + 
  geom_col(aes(Category, `Aligned Uniquely`)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  #scale_y_continuous(limits = c(0,60000000))
ggsave("figures/barplot_readsAlignedUniquely.png")
