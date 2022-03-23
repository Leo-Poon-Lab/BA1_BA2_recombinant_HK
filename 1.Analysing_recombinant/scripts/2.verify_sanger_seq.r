library(Biostrings)

seq_files <- list.files("../data/PCR_product_seqs/", ".seq$", full.names=T)
file_name <- list.files("../data/PCR_product_seqs/", ".seq$")
seqs <- lapply(seq_along(seq_files), function(i){
	# x <- seq_files[1]
	x <- seq_files[i]
	tmp <- readDNAStringSet(x)
	names(tmp) <- file_name[i] 
	if(grepl("R-", x) | grepl("R_", x)){
		tmp <- reverseComplement(tmp)
		names(tmp) <- paste0("RC_", names(tmp))
		}
	return(tmp)	
})
seqs <- do.call(c, seqs)

seq_name <- gsub("\\.seq", "", file_name)

names(seqs) <- seq_name
seqs <- seqs[order(names(seqs))]
writeXStringSet(seqs, "../results/sanger_seqs.fasta")

# align
system("mafft --thread -1 --localpair --maxiterate 1000 --adjustdirection --keeplength --add ../results/sanger_seqs.fasta ../results/representative.fasta > ../results/sanger_pcr_aligned.fasta")
