# Generate table of sequences used
# and generate GISAID acknowledgement table

library(tidyverse)
library(Biostrings)
library(writexl)

## sequences used in phylogenetic tree analysis
seqs_tree <- readDNAStringSet("../../1.Analysing_recombinant/results/seq_tree_fullgenome.fasta")
seq_name_helper <- names(seqs_tree)[!grepl("^WHP",names(seqs_tree))]
writeLines(seq_name_helper, "../results/seqname_tree_helper.txt")


## gisaid acknowledge table
df_gisaid_meta <- read_tsv("../results/gisaid_auspice_input_hcov-19_2022_03_23_07/1648021914078.metadata.tsv")
df_gisaid_meta <- df_gisaid_meta %>% arrange(originating_lab, submitting_lab, authors, gisaid_epi_isl)
df_gisaid_meta <- df_gisaid_meta %>% mutate(paste=paste0(originating_lab, submitting_lab, authors))
df_gisaid_meta_collapse <- df_gisaid_meta %>% group_by(paste) %>% summarise(`Accession ID`=paste(gisaid_epi_isl, collapse=", "), `Originating lab`=originating_lab[1],`Submitting lab`=submitting_lab[1], Authors=authors[1]) %>% ungroup()

write_xlsx(df_gisaid_meta_collapse %>% select(-paste), "../results/GISAID_table.xlsx")
