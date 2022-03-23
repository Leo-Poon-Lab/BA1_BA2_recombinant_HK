library(Biostrings)

seqs_consensus <- readDNAStringSet("../../0.Search_for_coinfection_recombinant/results/seqs_consensus_aln.fasta")
seqs_consensus <- seqs_consensus[c(4,3)]
names(seqs_consensus)[1] <- paste0("Reference:", names(seqs_consensus)[1])

seq_ref <- readDNAStringSet("/Volumes/GoogleDrive/My Drive/work/2020/2020-09-01_COVID_NGS_pipeline/NGS_data_input/reference.fasta")
seqs_parental <- readDNAStringSet("../results/seqs_closest_subset.fasta")
# names(seqs_parental) <- paste0("Parental_", LETTERS[seq_along(seqs_parental)], "(", names(seqs_parental), ")")
names(seqs_parental) <- gsub("_EPI", "_(EPI", names(seqs_parental))
names(seqs_parental) <- paste0(names(seqs_parental), ")")

seqs_consensus_new <- c(seq_ref, seqs_consensus, seqs_parental)
writeXStringSet(seqs_consensus_new, "../results/seqs_consensus_aln.fasta")
system("snipit ../results/seqs_consensus_aln.fasta -f pdf -d ../results/ ")
