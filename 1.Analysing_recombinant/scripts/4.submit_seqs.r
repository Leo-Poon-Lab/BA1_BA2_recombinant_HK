library(Biostrings)

seqs_used <- readDNAStringSet("../results/seq_tree_fullgenome.fasta")
list_seqs <- names(seqs_used)[grep("^WHP", names(seqs_used))]
list_seqs <- c(list_seqs, "WHP5870")
list_seqs[list_seqs=="WHP6494-ARCTIC"] <- "WHP6494-ARCTIC-S20-iseq"
writeLines(list_seqs, "../results/list_of_samples_to_submit.txt")

file_consensus <- list.files("../../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/ivar_consensus/")
file_consensus_full <- list.files("../../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/ivar_consensus/", full.names=T)

samples_consensus <- gsub("_consensus.fa","",file_consensus)

seqs_tosubmit <- lapply(file_consensus_full[samples_consensus %in% list_seqs],readDNAStringSet)
seqs_tosubmit <- do.call(c, seqs_tosubmit)
names(seqs_tosubmit) <- gsub("_consensus_threshold_.+_quality_30", "", names(seqs_tosubmit))
names(seqs_tosubmit) <- gsub("Consensus_", "", names(seqs_tosubmit))
writeXStringSet(seqs_tosubmit, "../results/seqs_tosubmit.fasta")
