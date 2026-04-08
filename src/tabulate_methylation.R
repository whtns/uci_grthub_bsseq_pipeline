#!/usr/bin/env Rscript
library(tidyverse)
library(fs)

library(plyranges)
# library(methrix) 
# library(BSgenome.Hsapiens.UCSC.hg38)

# List your files
bed_files <- dir_ls("output/methyldackel", glob = "*_CpG.bedGraph")  |> 
set_names(path_file)  |> 
set_names(str_remove, ".sorted_CpG.bedGraph")


# name: YourSeq
# location: chr11:34513780-34514215 (+)
# score: 1000.0
# score: matches = 436
# mismatches = 0
# repeat matches = 0
# # inserts in query = 0
# # inserts in target = 0
# --------------
# chr11:34513780-34514215

bed_counts  <- 
bed_files  |> 
    imap(~read_bed_graph(.x) %>% 
    as_tibble() %>%
            select(chrom = 1, start = 2, end = 3, {{.y}} := 6))

methylation_table <- 
    bed_counts  |> 
    purrr::reduce(full_join, by = c("chrom", "start", "end")) %>%
    mutate(across(starts_with("percent"), ~replace_na(.x, 0))) %>%
    arrange(chrom, start)  |> 
    dplyr::select(-`mR480-L1-PrNotRecog`)  |> 
    # dplyr::mutate(chrom = 11, 
    #        start = 34513780 + start - 2,
    #        end = 34513780 + end - 2)  |> 
           identity()

write_csv(methylation_table, "results/methylation_table_CpG.csv")

# 4. Plot with ggplot2
methylation_table  |> 
tidyr::pivot_longer(cols = -c(chrom, start, end), 
                    names_to = "sample", 
                    values_to = "percent")  |>
                    ggplot(aes(x = start, y = sample, color = percent)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Methylation Profiles", y = "Methylation %", x = "Position")

ggsave("results/methylation_profiles_CpG.png", width = 10, height = 6)
