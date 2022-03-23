# Aim: to build a tree showing the placement of partial genomes of the recombinant.
# 1. prepare background sequences for building tree.
# 2. add the closest matches.
# 3. build tree
# 4. Use the sequences for dating analysis

library(tidyverse)
library(Biostrings)
library(readxl)
library(writexl)
library(ggtree)
library(ggrepel)
library(shadowtext)
source("https://raw.githubusercontent.com/Koohoko/Save-ggplot-to-pptx/main/scripts/save_pptx.r")

# 1. prepare background sequences for building tree.
# system("./extract_seq_from_meta.sh")
df_meta_bg <- read_tsv("../data/nextstrain_ncov_open_reference_metadata.tsv")
seq_bg <- readDNAStringSet("../results/Genbank_lineage_representative.fasta")
df_meta_bg <- df_meta_bg %>% filter(!duplicated(pango_lineage))
sort(df_meta_bg$pango_lineage)
df_meta_bg <- df_meta_bg %>% filter(pango_lineage %in% c("B", "B.1.1.7", "B.1.351", "P.1", "B.1.617.2", "BA.3"))
seq_bg <- seq_bg[names(seq_bg) %in% df_meta_bg$strain]
write_xlsx(df_meta_bg, "../results/df_meta_bg.xlsx")
writeXStringSet(seq_bg, "../results/seqs_bg.fasta")

# 2. add the closest matches.
df_closest <- read_excel("../results/df_comp_closest.xlsx")
df_closest <- df_closest %>% filter(grepl("p[12]_GISAID", type))
table(df_closest$type)
## pick random four of the BA.2 cloest sequences
set.seed(2022)
(idx <- sample(4:9978, 3))
df_closest_subset <- df_closest %>% arrange(type, num_diff) %>% .[c(1:3, 10000+idx),]
df_closest_subset$mutations
write_xlsx(df_closest, "../results/df_meta_closest.xlsx")
write_xlsx(df_closest_subset, "../results/df_meta_closest_subset.xlsx")

files_seqs_closest <- list.files("../results/", "GISAID_closest_n10000_alnseq.fasta", full.names=T)
seqs_all <- lapply(files_seqs_closest, readDNAStringSet)
seqs_all <- do.call(c, seqs_all)
seqs_closest <- seqs_all[names(seqs_all) %in% df_closest$sample]
seqs_closest <- seqs_closest[!duplicated(names(seqs_closest))]
writeXStringSet(seqs_closest, "../results/seqs_closest.fasta")
seqs_closest_subset <- seqs_closest[names(seqs_closest) %in% df_closest_subset$sample]

# names(seqs_closest_subset) <- paste0("Parental_", LETTERS[seq_along(seqs_closest_subset)], "_(", c(rep("BA.1_",3), rep("BA.2_",3)), names(seqs_closest_subset), ")")
names(seqs_closest_subset) <- paste0(c(rep("BA.1_",3), rep("BA.2_",3)), names(seqs_closest_subset))
writeXStringSet(seqs_closest_subset, "../results/seqs_closest_subset.fasta")

# 3. add other imported cases
seqs_aln_all <- readDNAStringSet("/Volumes/GoogleDrive/My Drive/work/2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/ivar_consensus/con_all_combined_aln.fasta")
df_meta_imported <- readxl::read_excel("../../0.Search_for_coinfection_recombinant/results/df_metadata_hk_import.xlsx")
df_meta_imported <- df_meta_imported %>% filter(lineage!="None")
df_problem <- read_tsv("../../0.Search_for_coinfection_recombinant/results/problematic_samples.tsv")

df_meta_imported$Sample <- gsub("_\\d+", "", df_meta_imported$Sample)
df_meta_imported <- df_meta_imported %>% filter(!Sample %in% df_problem$sample)
seqs_aln_imported <- seqs_aln_all[names(seqs_aln_all) %in% df_meta_imported$Sample]

# 4. build tree
seq_recombinant <- readBStringSet("../results/representative_all.fasta")
seq_recombinant_partial <- seq_recombinant[grepl("part", names(seq_recombinant))]
seq_tree <- c(seq_bg, seqs_closest_subset, seq_recombinant, seqs_aln_imported)
# seq_tree <- c(seq_bg, seqs_aln_imported, seq_recombinant)

tmp <- rep(FALSE, width(seq_tree)[1])
tmp[1:265] <- TRUE # mask UTR
tmp[29645:29903] <- TRUE # mask UTR
at <- matrix(rep(tmp, length(seq_tree)),
            nrow=length(seq_tree), ncol=width(seq_tree)[1], byrow=TRUE)
letter_subject <- DNAString(paste(rep.int("N", width(seq_tree)[1]), collapse=""))
letter <- as(Views(letter_subject, start=1, end=rowSums(at)), "XStringSet")
seq_tree <- replaceLetterAt(seq_tree, at, letter)

seq_few <- seq_tree[!names(seq_tree) %in% names(seqs_aln_imported)]
writeXStringSet(seq_few[!grepl("WHP6494-ARCTIC", names(seq_few))], "../results/seq_few_partial.fasta")
system("~/softwares/iqtree2/bin/iqtree2 --redo -s ../results/seq_few_partial.fasta --alrt 1000 -B 1000 -T 8 -o 'Wuhan-Hu-1/2019'")

options(scipen=999) 
plot_tree <- function(tree_file) {
	tree <- read.tree(tree_file)
	p <- ggtree(tree, size=0.3)
	p$data <- left_join(p$data, (bind_rows(df_meta_bg %>% mutate(label=strain)) %>% select(label, pango_lineage)), "label")
	p$data <- left_join(p$data, (bind_rows(df_meta_imported %>% mutate(label=Sample)) %>% select(label, lineage)), "label")
	# p$data$label <- gsub("__","_(",p$data$label)
	p$data$label <- gsub("_EPI",": EPI",p$data$label)
	p$data$label[grepl("partial1", p$data$label)] <- "Recombinant (5' end)"
	p$data$label[grepl("partial2", p$data$label)] <- "Recombinant (3' end)"
	p$data$label[grepl("WHP6494-ARCTIC", p$data$label)] <- "Recombinant_full_genome"

	p$data$label_new <- paste0(p$data$lineage, " (", p$data$label, ")")
	p$data$label_new2 <- paste0(p$data$pango_lineage, " (", p$data$label, ")")
	
	# bootstrap
	p$data$shalrt <- NA 
	p$data$ufboost <- NA 
	p$data$shalrt[!p$data$isTip] <- sapply(p$data$label[!p$data$isTip], function(x) {
		strsplit(x, "\\/")[[1]][1]
	})
	p$data$shalrt <- as.numeric(p$data$shalrt)
	p$data$ufboost[!p$data$isTip] <- sapply(p$data$label[!p$data$isTip], function(x) {
		strsplit(x, "\\/")[[1]][2]
	})
	p$data$ufboost <- as.numeric(p$data$ufboost)
	p$data$node_support <- (p$data$ufboost>=95) & (p$data$shalrt>=80)
	
	label_size=5
	p_out <- p + 
		geom_text_repel(aes(x=x, y=y, label=label), data=p$data %>% filter(node_support), color="Dark red",size=4, bg.color="white", nudge_x=-0.00023, nudge_y=0.5, min.segment.length=2, segment.size=0.2, segment.color="red")+
		geom_nodepoint(aes(x=x, y=y), data=p$data %>% filter(node_support), size=3, fill="Dark red", color="Dark red", alpha=0.8)+
		geom_tiplab(data=. %>% filter(grepl("Recombinant", label)), color="Dark red", size=label_size, offset=0.00006)+
		geom_tiplab(data=. %>% filter(grepl("^BA", label)), size=label_size, offset=0.00006)+
		geom_tiplab(aes(label=label_new), data=. %>% filter(!is.na(lineage)), color="Dark blue", size=label_size, offset=0.00006)+
		geom_tiplab(aes(label=label_new2), data=. %>% filter(!is.na(pango_lineage)), size=label_size, offset=0.00006)+
		geom_treescale(x=0.003, y=3, width=0.0005, fontsize=3)+
		xlim(0,0.0038)+
		NULL
	return(p_out)
}

p_out <- plot_tree("../results/seq_few_partial.fasta.treefile")
ggsave("../results/tree_few_partial.pdf", width=7.5, height=7.5*sqrt(2), plot=p_out)
save_pptx("../results/tree_few_partial.pptx", width=7.5, height=7.5*sqrt(2), plot=p_out)

# seq more
writeXStringSet(seq_tree[!grepl("WHP6494-ARCTIC", names(seq_tree))], "../results/seq_tree_partial.fasta")
writeXStringSet(seq_tree[!grepl("partial", names(seq_tree))], "../results/seq_tree_fullgenome.fasta")
system("~/softwares/iqtree2/bin/iqtree2 --redo -s ../results/seq_tree_partial.fasta --alrt 1000 -B 1000 -T 8 -o 'Wuhan-Hu-1/2019'")
system("~/softwares/iqtree2/bin/iqtree2 --redo -s ../results/seq_tree_fullgenome.fasta --alrt 1000 -B 1000 -T 8 -o 'Wuhan-Hu-1/2019'")

options(scipen=999) 
plot_tree <- function(tree_file) {
	tree <- read.tree(tree_file)
	p <- ggtree(tree, size=0.3)
	p$data <- left_join(p$data, (bind_rows(df_meta_bg %>% mutate(label=strain)) %>% select(label, pango_lineage)), "label")
	p$data <- left_join(p$data, (bind_rows(df_meta_imported %>% mutate(label=Sample)) %>% select(label, lineage)), "label")
	# p$data$label <- gsub("__","_(",p$data$label)
	p$data$label <- gsub("_EPI",": EPI",p$data$label)
	p$data$label[grepl("partial1", p$data$label)] <- "Recombinant (5' end)"
	p$data$label[grepl("partial2", p$data$label)] <- "Recombinant (3' end)"
	p$data$label[grepl("WHP6494-ARCTIC", p$data$label)] <- "Recombinant_full_genome"

	p$data$label_new <- paste0(p$data$lineage, " (", p$data$label, ")")
	p$data$label_new2 <- paste0(p$data$pango_lineage, " (", p$data$label, ")")
	
	# bootstrap
	p$data$shalrt <- NA 
	p$data$ufboost <- NA 
	p$data$shalrt[!p$data$isTip] <- sapply(p$data$label[!p$data$isTip], function(x) {
		strsplit(x, "\\/")[[1]][1]
	})
	p$data$shalrt <- as.numeric(p$data$shalrt)
	p$data$ufboost[!p$data$isTip] <- sapply(p$data$label[!p$data$isTip], function(x) {
		strsplit(x, "\\/")[[1]][2]
	})
	p$data$ufboost <- as.numeric(p$data$ufboost)
	p$data$node_support <- (p$data$ufboost>=95) & (p$data$shalrt>=80)
	
	label_size=1.5
	annt_text_size=1.5
	offset_label=0.00002
	p_out <- p + 
		geom_text_repel(aes(x=x, y=y, label=label), data=p$data %>% filter(node_support), color="Dark red",size=annt_text_size, bg.color="white", nudge_x=-0.00005, nudge_y=0.1, min.segment.length=0.01, segment.size=0.2, segment.color="red")+
		geom_nodepoint(aes(x=x, y=y), data=p$data %>% filter(node_support), size=1, fill="Dark red", color="Dark red", alpha=0.8)+
		geom_tiplab(data=. %>% filter(grepl("Recombinant", label)), color="Dark red", size=label_size, offset=offset_label)+
		geom_tiplab(data=. %>% filter(grepl("^BA", label)), size=label_size, offset=offset_label)+
		geom_tiplab(aes(label=label_new), data=. %>% filter(!is.na(lineage)), color="Dark blue", size=label_size, offset=offset_label)+
		geom_tiplab(aes(label=label_new2), data=. %>% filter(!is.na(pango_lineage)), size=label_size, offset=offset_label)+
		geom_treescale(x=0.002, y=3, width=0.0005, fontsize=3)+
		xlim(0,0.0027)+
		NULL
	return(p_out)
}

p_out <- plot_tree("../results/seq_tree_fullgenome.fasta.treefile")
ggsave("../results/tree_fullgenome.pdf", width=7.5, height=7.5*sqrt(2), plot=p_out)
save_pptx("../results/tree_fullgenome.pptx", width=7.5, height=7.5*sqrt(2), plot=p_out)

p_out <- plot_tree("../results/seq_tree_partial.fasta.treefile")
ggsave("../results/tree_partial.pdf", width=7.5, height=7.5*sqrt(2), plot=p_out)
save_pptx("../results/tree_partial.pptx", width=7.5, height=7.5*sqrt(2), plot=p_out)
