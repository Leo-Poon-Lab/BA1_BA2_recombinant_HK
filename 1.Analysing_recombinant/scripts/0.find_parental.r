# ## Analyzing steps
# 1. Visualize the mutations, MAF, Recombination. 
# 2. Blast the partial sequence to find parental sequences.
# 3. Run 3SEQ to get statistics supporting recombination.

library(Biostrings)
library(tidyverse)
library(naturalsort)
library(ggsci)
library(patchwork)
library(writexl)
source("https://raw.githubusercontent.com/Koohoko/Save-ggplot-to-pptx/main/scripts/save_pptx.r")
# 1. Visualize the mutations, MAF, Recombination. 
### generate the table for comparision
samples <- c("WHP5870", "WHP6494-AR")

files_bamstat <- list.files("../../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/pysamstats/", full.names=T)
check <- sapply(samples, function(x) {
	tmp <- grep(x, files_bamstat)
	if(length(tmp)==0){return(NA)}else{return(tmp)}
})
df_bamstat <- lapply(files_bamstat[check[!is.na(check)]], read_tsv)
df_bamstat <- bind_rows(df_bamstat)
df_bamstat$sample <- rep(samples[!is.na(check)], each=29903)

source("./helper/af_by_bamstat.r")
df_ba12 <- read_csv("../data/Omicron_BA.1_BA.2_mutations.csv")
# df_ba12 <- df_ba12 %>% filter(is.na(`B.1.1.529`)) %>% select(-`B.1.1.529`)
split_tmp <- strsplit(df_ba12$`mut (nuc)`, "[,] ")
df_vcf_int <- tibble(mutation_aa=rep(df_ba12$`mutation (aa)`, sapply(split_tmp, length)),`mut (nuc)`=rep(df_ba12$`mut (nuc)`, sapply(split_tmp, length)), mutation_nt=unlist(split_tmp))
df_vcf_int <- left_join(df_vcf_int, df_ba12 %>% mutate(mutation_aa=`mutation (aa)`) %>% select(mutation_aa, gene, `BA.2`:notes))

df_vcf_int$POS <- gsub("[_ ].+$", "", df_vcf_int$`mutation_nt`)
df_vcf_int <- df_vcf_int %>% arrange(POS)
df_vcf_int$mutation_nt_sim <- gsub(" \\(.+", "", df_vcf_int$`mutation_nt`)
df_vcf_int$POS <- gsub("\\D", "", df_vcf_int$POS)
df_vcf_int$POS <- as.numeric(df_vcf_int$POS)
df_vcf_int$REF <- sapply(strsplit(df_vcf_int$`mutation_nt_sim`, "\\d"), function(x){x[1]})
df_vcf_int$ALT <- sapply(strsplit(df_vcf_int$`mutation_nt_sim`, "\\d"), function(x){x[length(x)]})
df_vcf_int$ALT[grepl("_", df_vcf_int$`mutation_nt_sim`)|grepl("del", tolower(df_vcf_int$`mutation_nt_sim`))] <- "deletions"
df_vcf_int$ALT[grepl("ins", df_vcf_int$`mutation_nt_sim`)] <- "insertions"
# df_vcf_int <- df_vcf_int[-c(19,20), ]

df_af <- lapply(naturalsort(samples), function(x) {
	tmp <- af_by_bamstat(df_bamstat, x, df_vcf_int$POS, df_vcf_int$REF, df_vcf_int$ALT)
	tmp$Alt_freq <- tmp$ALT_depth/tmp$Total_depth
	tmp$MAF <- cal_maf(df_bamstat, x, df_vcf_int$POS)

	df_input_check_valid <- df_vcf_int %>% select(`mut (nuc)`, `BA.2`:`B.1.1.529`)
	df_tmp <- filter_valid(df=df_input_check_valid, af_tmp=tmp, variant="all", min_depth=10, min_AF=0.1)

	df_tmp$sample <- x
	return(df_tmp)
})
df_af <- bind_rows(df_af)
df_af$Alt_freq <- df_af$ALT_depth/df_af$Total_depth
df_af <- left_join(df_af, df_vcf_int %>% select(mutation_aa:gene, POS:ALT), "mut (nuc)")
sort(df_af$MAF)
writexl::write_xlsx(df_af, "../results/df_af.xlsx")

## visualize

df_plot <- df_af
# df_plot <- df_plot %>% filter(`mut (nuc)` != "C22674T") # remove confusing snp
### cal maf
sort(df_plot$MAF[df_plot$Alt_freq>0.9])
mean(df_plot$MAF[df_plot$Alt_freq>0.9], na.rm=T)*100
median(df_plot$MAF[df_plot$Alt_freq>0.9], na.rm=T)*100

df_plot <- df_plot %>% filter(is.na(B.1.1.529))
df_plot <- df_plot %>% arrange(POS)
df_plot$mutation_aa <- factor(df_plot$mutation_aa, levels=unique(df_plot$mutation_aa))
df_plot$`mut (nuc)` <- factor(df_plot$`mut (nuc)`, levels=unique(df_plot$`mut (nuc)`))
df_plot$mut_freq <- sapply(df_plot$Alt_freq, function(x){
	if(is.na(x)){return("No coverage")}
	if(x>0.99){
		return("> 0.99")
	} else if(x>0.9){
		return("> 0.90")
	} else if(x<0.02){
		return("< 0.02")
	} else {
		return("< 0.9")
	}
})

df_plot <- df_plot %>% pivot_longer(cols=c(`BA.1`, `BA.2`, `mut_freq`))
df_plot <- df_plot %>% filter(!is.na(value))
df_plot$value[df_plot$name %in% c("BA.1", "BA.2")] <- df_plot$name[df_plot$name %in% c("BA.1", "BA.2")]
df_plot$name[df_plot$name =="mut_freq"] <- df_plot$sample[df_plot$name =="mut_freq"]
df_plot$name[df_plot$name =="WHP5870"] <- "Case-patient B"
df_plot$name[df_plot$name =="WHP6494-AR"] <- "Case-patient A"

df_plot$name[df_plot$name %in% c("BA.1", "BA.2")] <- "BA.1/BA.2\ndefining mutations"
types <- c("BA.1", "BA.2", "> 0.99", "< 0.02", "> 0.98")
df_plot$value <- factor(df_plot$value, levels = types)
df_plot$gene_facet <- factor(df_plot$gene, unique(c("nuc", df_plot$gene)))
levels(df_plot$gene_facet)[1] <- "Not in ORF"

colors <- c("Dark red", "Dark blue", "Dark green", "#bfd9d8", "White")
names(colors) <- types

p1 <- df_plot %>% ggplot() +
	geom_tile(aes(x=name, y=`mut (nuc)`, fill=value), color="black")+
	scale_fill_manual(name="Variant defining mutations &\nSample mutation frequency", values=colors, na.value = "white")+
	theme_classic()+
	xlab("")+
	ylab("Mutations (Nt)")+
	# facet_grid(rows=vars(gene_facet), scales="free", space="free")+
	theme(legend.position='top', axis.text.y=element_text(colour=ifelse(levels(df_plot$`mut (nuc)`) %in% c('A20055G', 'C21618T'), 'red', 'black')))+
	NULL

colors_gene <- pal_jco()(length(levels(df_plot$gene_facet)))
names(colors_gene) <- levels(df_plot$gene_facet)
colors_gene[colors_gene=="#868686FF"] <- colors_gene[names(colors_gene)=="Not in ORF"]
colors_gene[names(colors_gene)=="Not in ORF"] <- "grey40"

p2 <- df_plot %>% filter(name=="Case-patient A") %>% mutate(name="Gene") %>% ggplot() +
	geom_tile(aes(x=name, y=mutation_aa, fill=gene_facet), color="black")+
	scale_fill_manual(name="Gene", values=colors_gene, na.value = "white")+
	theme_classic()+
	xlab("")+
	ylab("Mutations (AA)")+
	scale_y_discrete(position = "right")+
	theme(
		panel.background = element_rect(fill='transparent'),
		plot.background = element_rect(fill='transparent', color=NA),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.text.y=element_text(colour=ifelse(levels(df_plot$`mut (nuc)`) %in% c('A20055G', 'C21618T'), 'red', 'black'))
	)+
	NULL

p1_2 <- p1+p2+plot_layout(widths = c(3, 0.3), guides="keep")
ggsave("../results/fig1.pdf", width = 8, height=8, plot=p1_2)
save_pptx("../results/fig1.pptx", width = 8, height=8, plot=p1_2)

# 2. Blast the full sequence to find similar sequences.
# 2. Blast the partial sequence to find parental sequences.
## breaking point 20055 to 21618 
source("./helper/helper_comp_seq.r")

seq_ori <- readDNAStringSet("../data/comp.fasta")
df_comp_2patient <- comp_seqs(seq_ori[grep("WHP6494-ARCTIC", names(seq_ori))], seq_ori[grep("WHP5870", names(seq_ori))])
write_xlsx(df_comp_2patient, "../results/df_comp_2patient.xlsx")
comp_seqs(seq_ori[grep("WHP6494-ARCTIC", names(seq_ori))], seq_ori[grep("WHP5870", names(seq_ori))]) %>% .$mutations
comp_seqs(seq_ori[grep("WHP6494-ARCTIC", names(seq_ori))], seq_ori[grep("WHP6494-S24", names(seq_ori))]) %>% .$mutations
comp_seqs(seq_ori[grep("WHP6494-ARCTIC", names(seq_ori))], seq_ori[grep("261", names(seq_ori))]) %>% .$mutations

## substitute the ambiguous nts according to other sequences
subseq(seq_ori[grep("WHP6494-ARCTIC", names(seq_ori))], 15521, 15521) <- "T"
subseq(seq_ori[grep("WHP6494-ARCTIC", names(seq_ori))], 22752, 22752) <- "C"

seq_ori <- seq_ori[c(grep("WHP6494-ARCTIC", names(seq_ori)), grep("WHP6494-S24", names(seq_ori)))]
names(seq_ori) <- c("WHP6494-ARCTIC", "WHP6494-S24")
writeXStringSet(seq_ori, "../results/representative_ori.fasta") # in the final version, we use WHP6494 as the representative

seq_int <- seq_ori[1]
## fill gaps with another sequence
check_gap_us <- which(strsplit(as.character(seq_int), "")[[1]] %in% c("N", "-"))
check_gap_chp <- which(strsplit(as.character(seq_ori[2]), "")[[1]] %in% c("N", "-"))
check_gap_tofill <- check_gap_us[!check_gap_us %in% check_gap_chp]
sapply(check_gap_tofill, function(x){
	subseq(seq_int, x, x) <<- subseq(seq_ori[2], x, x)
	return(NA)
})

tmp <- rep(FALSE, width(seq_int)[1])
tmp[1:265] <- TRUE # mask UTR
tmp[29645:29903] <- TRUE # mask UTR
at <- matrix(rep(tmp, length(seq_int)),
            nrow=length(seq_int), ncol=width(seq_int)[1], byrow=TRUE)
letter_subject <- DNAString(paste(rep.int("N", width(seq_int)[1]), collapse=""))
letter <- as(Views(letter_subject, start=1, end=rowSums(at)), "XStringSet")
seq_int <- replaceLetterAt(seq_int, at, letter)

tmp <- rep(FALSE, width(seq_int)[1])
tmp[20055:width(seq_int)] <- TRUE
at <- matrix(rep(tmp, length(seq_int)),
            nrow=length(seq_int), ncol=width(seq_int)[1], byrow=TRUE)
letter_subject <- DNAString(paste(rep.int("N", width(seq_int)[1]), collapse=""))
letter <- as(Views(letter_subject, start=1, end=rowSums(at)), "XStringSet")
seq_partial1 <- replaceLetterAt(seq_int, at, letter)
names(seq_partial1) <- "representative_partial1"

tmp <- rep(FALSE, width(seq_int)[1])
tmp[1:21618] <- TRUE
at <- matrix(rep(tmp, length(seq_int)),
            nrow=length(seq_int), ncol=width(seq_int)[1], byrow=TRUE)
letter_subject <- DNAString(paste(rep.int("N", width(seq_int)[1]), collapse=""))
letter <- as(Views(letter_subject, start=1, end=rowSums(at)), "XStringSet")
seq_partial2 <- replaceLetterAt(seq_int, at, letter)
names(seq_partial2) <- "representative_partial2"

writeXStringSet(seq_int, "../results/representative.fasta")
writeXStringSet(seq_partial1, "../results/representative_p1.fasta")
writeXStringSet(seq_partial2, "../results/representative_p2.fasta")
writeXStringSet(c(seq_int, seq_partial1, seq_partial2), "../results/representative_all.fasta")

## get mut sites of representative
seq_ref <- readDNAStringSet("/Volumes/GoogleDrive/My Drive/work/2020/2020-09-01_COVID_NGS_pipeline/NGS_data_input/reference.fasta")
df_comp <- comp_seqs(seq_ref, seq_int, ignore_nt=c("N", "?"))
mut_sites <- df_comp$mutations
mut_sites <- gsub("\\D", "", strsplit(mut_sites, "|", fixed=T)[[1]])
writeLines(mut_sites, "../results/representative_mutsites.txt")
write_xlsx(df_comp, "../results/df_mutations_representative.xlsx")

## minimap2
system("chmod 755 ./minimap2.sh")
system("./minimap2.sh 2> minimap.err")

## mask UTR/recomb with Ns and find closest
system("chmod 755 ./find_closest.sh")
system("./find_closest.sh 2> find_closest.err")

files_metadata <- list.files("../results/", "n10000_metadata", full.names=T)
df_comp_all <- mclapply(files_metadata, function(x){
	file_prefix <- gsub(pattern="_G.+", replacement="", x)
	file_refseq <- paste0(file_prefix, "_masked.fasta")
	file_maskseq <- gsub("metadata.tsv", "maskseq.fasta", x)
	file_alnseq <- gsub("metadata.tsv", "alnseq.fasta", x)
	seq_ref <- readBStringSet(file_refseq)
	seq_mask <- readBStringSet(file_maskseq)
	# seq_aln <- readDNAStringSet(file_alnseq)
	df_meta <- read_tsv(x)
	if(grepl("GISAID", x)){
		df_meta_sim <- df_meta %>% select(`Accession ID`, `Pango lineage`, `Collection date`, `Location`)
		df_meta_sim$Location <- sapply(strsplit(df_meta_sim$Location, " / "), function(x){x[2]})
	}else{
		df_meta_sim <- df_meta %>% select(strain, `pango_lineage`, date, country)
	}
	names(df_meta_sim) <- c("sample", "pango_lineage", "date", "Country")
	df_meta_sim$date <- lubridate::ymd(df_meta_sim$date)
	df_comp <- comp_seqs(seq_ref, seq_mask)
	df_comp$reference <- names(seq_ref)
	df_comp$type <- gsub(".+\\/", "", file_maskseq) %>% gsub(pattern="_mask.+", replacement="", x)
	return(left_join(df_comp, df_meta_sim, "sample"))
}, mc.cores=8)

df_comp_all <- bind_rows(df_comp_all)

ggplot(df_comp_all)+
	geom_histogram(aes(x=num_diff_excluding_primer_binding))+
	facet_wrap(vars(type), ncol=1)

write_xlsx(df_comp_all, "../results/df_comp_closest.xlsx")

df_comp_all_sum <- df_comp_all %>% filter((grepl("p1", type) & num_diff<=1) | (grepl("p2", type) & num_diff==0)) %>% group_by(type, num_diff, mutations, Country) %>% summarise(N=n(), date_range=paste0(min(date, na.rm=T), " to ", max(date, na.rm=T))) %>% arrange(type, num_diff, desc(N)) %>% ungroup()
df_comp_all_sum$Genome_region <- ifelse(grepl("p1", df_comp_all_sum$type), "BA.1", "BA.2")
df_comp_all_sum$Database <- ifelse(grepl("GISAID", df_comp_all_sum$type), "GISAID", "Genbank")
df_comp_all_sum <- df_comp_all_sum %>% select(-type) %>% select(Genome_region, Database, everything())
write_xlsx(df_comp_all_sum, "../results/df_comp_all_sum.xlsx")

df_comp_all_sum %>% filter(Genome_region=="BA.2") %>% group_by(Database) %>% summarise(sum=sum(N))
df_comp_all_sum %>% filter(Genome_region=="BA.2") %>% summarise(sum=sum(N))
df_comp_all_sum %>% filter(Genome_region=="BA.2") %>% .$Country %>% unique() %>% length()

df_comp_all %>% filter(grepl("p1", type)) %>% .$num_diff %>% table()
df_comp_all %>% filter(grepl("p1", type)) %>% .$mutations %>% table()
df_comp_all %>% filter(grepl("p2", type)) %>% .$num_diff %>% table()
df_comp_all %>% filter(grepl("p2", type)) %>% group_by(type) %>% filter(num_diff==0) %>% summarise(n())
df_comp_all %>% filter(grepl("p1", type) & num_diff_excluding_primer_binding==1) %>% .$sample
df_comp_all %>% filter(grepl("p1", type) & num_diff_excluding_primer_binding==1) %>% .$mutations


# 3. Run 3SEQ to get statistics supporting recombination.
