# 1. plot heat map of mutation frequency to confirm connection
# 2. visualize consensus sequence to find potential recombinant
# 3. manual inspection of the raw reads in igv

library(tidyverse)
library(readxl)
library(ggsci)
library(patchwork)
library(Biostrings)
library(ggnewscale)
source("https://raw.githubusercontent.com/Koohoko/Save-ggplot-to-pptx/main/scripts/save_pptx.r")

# 1. plot heat map of mutation frequency to confirm connection
df_mut_freq <- read_csv("../results/mut_freq_of_potential_coinfection.csv")
df_vcf_int <- read_excel("../data/Omicron_BA.1_BA.2_mutations_cleaned.xlsx")
df_plot <- left_join(df_mut_freq, df_vcf_int %>% select(`mut (nuc)`, mutation_aa, gene) %>% unique())
df_plot <- df_plot %>% filter(`mut (nuc)` != "C22674T") # remove confusing snp

df_plot <- df_plot %>% filter(!grepl("WHP6494-S24", sample)) #duplicated sample
samples_coinfection <- unique(df_plot$sample)
df_plot <- df_plot %>% filter(!grepl("WHP6494", sample)) # only plot co-infection
df_plot <- df_plot %>% filter(!grepl("WHP5870", sample)) # only plot co-infection

df_plot$sample <- gsub("-.+", "", df_plot$sample)

df_plot <- df_plot %>% filter(is.na(B.1.1.529))
df_plot <- df_plot %>% select(-`B.1.1.529`)

# df_plot$mutation_aa <- factor(df_plot$mutation_aa, levels=unique(df_plot %>% arrange(`BA.1`) %>% .$mutation_aa))
# df_plot$`mut (nuc)` <- factor(df_plot$`mut (nuc)`, levels=unique(df_plot %>% arrange(`BA.1`) %>% .$`mut (nuc)`))
df_plot$mutation_aa <- factor(df_plot$mutation_aa, levels=unique(df_plot$mutation_aa))
df_plot$`mut (nuc)` <- factor(df_plot$`mut (nuc)`, levels=unique(df_plot$`mut (nuc)`))

func_Y_to_1 <- function(x) {
	x <- ifelse(x=="Y", 1, 0)
	x[is.na(x)] <- 0
	return(x)
}
df_plot <- df_plot %>% mutate_at(vars(`BA.2`:`BA.1`), func_Y_to_1)

df_plot <- df_plot %>% pivot_longer(cols=c(`BA.1`, `BA.2`, `Alt_freq`))
df_plot <- df_plot %>% filter(!is.na(value))
df_plot$name[df_plot$name =="Alt_freq"] <- df_plot$sample[df_plot$name =="Alt_freq"]

df_plot$gene_facet <- factor(df_plot$gene, unique(c("nuc", df_plot$gene)))
levels(df_plot$gene_facet)[1] <- "Not in ORF"

df_plot$value_break <- cut(df_plot$value, c(0.02, seq(0,100,20)/100), include.lowest = T)

p1 <- df_plot %>% ggplot() +
	geom_tile(aes(x=name, y=`mut (nuc)`, fill=value_break), color="black")+
	scale_fill_viridis_d(name="Sample mutation frequency", na.value = "white")+
	theme_classic()+
	xlab("")+
	# scale_x_discrete(guide = guide_axis(n.dodge=2))+
	ylab("Mutations (Nt)")+
	# facet_grid(rows=vars(gene_facet), scales="free", space="free")+
	theme(legend.position='top')+
	NULL
ggsave("../results/plot_mut_freq_coinfection_check_v1.pdf", width = 16, height=8, plot=p1)

## after checking, exclude the samples with possible lab contamination
df_plot <- df_plot %>% filter(!name %in% c("WHP5209", "WHP5407"))
samples_coinfection <- samples_coinfection[!samples_coinfection %in% c("WHP5209", "WHP5407")]

muts_BA1 <- df_plot %>% filter(name=="BA.1" & value==1) %>% .$`mut (nuc)` %>% unique() %>% as.character()
muts_BA2 <- df_plot %>% filter(name=="BA.2" & value==1) %>% .$`mut (nuc)` %>% unique() %>% as.character()

df_plot$name[df_plot$name %in% c("BA.1", "BA.2")] <- "BA.1/BA.2\ndefining mutations"
df_plot$`mut (nuc)`

p1 <- df_plot %>% ggplot() +
	geom_tile(aes(x=name, y=`mut (nuc)`, fill=value), color="black", data=. %>% filter(`mut (nuc)`%in%muts_BA1))+
	scale_y_discrete(drop = FALSE)+
	scale_fill_gradient2(name="Mutation frequency (BA.1)", low="#FFFFFF00", mid="#F4001B50", high="#F4001B", midpoint=0.3)+ # BA.1
	new_scale_fill()+
	geom_tile(aes(x=name, y=`mut (nuc)`, fill=value), color="black", data=. %>% filter(`mut (nuc)`%in%muts_BA2))+
	scale_fill_gradient2(name="Mutation frequency (BA.2)", low="#FFFFFF00", mid="#88CD4650", high="#88CD46", midpoint=0.3)+ # BA.2
	# scale_fill_viridis_d(name="Sample mutation frequency", na.value = "white")+
	geom_tile(aes(x=name, y=`mut (nuc)`), fill="white", color="black", data=. %>% filter(grepl("WHP", name) & round(value, 2)==0))+
	geom_text(aes(x=name, y=`mut (nuc)`, label=round(value, 2)), size=3, color="black", data=. %>% filter(grepl("WHP", name)))+
	theme_classic()+
	xlab("")+
	# scale_x_discrete(guide = guide_axis(n.dodge=2))+
	ylab("Mutations (Nt)")+
	# facet_grid(rows=vars(gene_facet), scales="free", space="free")+
	theme(legend.position='top', axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
	NULL

colors_gene <- pal_jco()(length(levels(df_plot$gene_facet)))
names(colors_gene) <- levels(df_plot$gene_facet)
colors_gene[colors_gene=="#868686FF"] <- colors_gene[names(colors_gene)=="Not in ORF"]
colors_gene[names(colors_gene)=="Not in ORF"] <- "grey40"

p2 <- df_plot %>% filter(name=="BA.1") %>% mutate(name="Gene") %>% ggplot() +
	geom_tile(aes(x=name, y=mutation_aa, fill=gene_facet), color="black")+
	scale_fill_manual(name="Gene", values=colors_gene, na.value = "white")+
	theme_classic()+
	xlab("")+
	ylab("Mutations (AA)")+
	scale_y_discrete(position = "right")+
	theme(
		# axis.text.y=element_text(colour=ifelse(levels(df_plot$`mut (nuc)`) %in% c('A20055G', 'C21618T'), 'red', 'black'),
		panel.background = element_rect(fill='transparent'),
		plot.background = element_rect(fill='transparent', color=NA),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank()
	)+
	NULL

p1_2 <- p1+p2+plot_layout(widths = c(5, 1), guides="keep")
ggsave("../results/plot_mut_freq_coinfection.pdf",  width = 6.5, height=7.5, plot=p1_2)
save_pptx("../results/plot_mut_freq_coinfection.pptx", width = 6.5, height=7.5, plot=p1_2)


# 2. visualize consensus sequence to find potential recombinant
## prepare alignment for sequences
seqs_aln_all <- readDNAStringSet("/Volumes/GoogleDrive/My Drive/work/2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/ivar_consensus/con_all_combined_aln.fasta")
samples_recombinant <- samples_coinfection[grepl("WHP5870",samples_coinfection) | grepl("WHP6494",samples_coinfection)]
seqs_aln_recombinant <- seqs_aln_all[names(seqs_aln_all) %in% samples_recombinant]

names(seqs_aln_recombinant) <- gsub("-.+", "", names(seqs_aln_recombinant))
names(seqs_aln_recombinant)[names(seqs_aln_recombinant)=="WHP5870"] <- "WHP5870_(case-patient_B)"
names(seqs_aln_recombinant)[names(seqs_aln_recombinant)=="WHP6494"] <- "WHP6494_(case-patient_A)"

seq_ref <- readDNAStringSet("/Volumes/GoogleDrive/My Drive/work/2020/2020-09-01_COVID_NGS_pipeline/NGS_data_input/reference.fasta")
seqs_bg <- readDNAStringSet("../data/BA1_2_nextclade.aligned.fasta")

seqs_consensus <- c(seqs_bg, seqs_aln_recombinant)
tmp <- rep(FALSE, width(seqs_consensus)[1])
tmp[1:265] <- TRUE # mask UTR
tmp[29645:29903] <- TRUE # mask UTR
at <- matrix(rep(tmp, length(seqs_consensus)),
            nrow=length(seqs_consensus), ncol=width(seqs_consensus)[1], byrow=TRUE)
letter_subject <- DNAString(paste(rep.int("N", width(seqs_consensus)[1]), collapse=""))
letter <- as(Views(letter_subject, start=1, end=rowSums(at)), "XStringSet")
seqs_consensus <- replaceLetterAt(seqs_consensus, at, letter)

writeXStringSet(seqs_consensus, "../results/seqs_consensus_aln.fasta")
system("snipit ../results/seqs_consensus_aln.fasta -d ../results/ ")


# plot for controls
df_mut_freq <- read_csv("../results/mut_freq_of_controls.csv")
df_vcf_int <- read_excel("../data/Omicron_BA.1_BA.2_mutations_cleaned.xlsx")
df_plot <- left_join(df_mut_freq, df_vcf_int %>% select(`mut (nuc)`, mutation_aa, gene) %>% unique())
df_plot <- df_plot %>% filter(`mut (nuc)` != "C22674T") # remove confusing snp

df_plot <- df_plot %>% filter(is.na(B.1.1.529))
df_plot <- df_plot %>% select(-`B.1.1.529`)

df_plot$mutation_aa <- factor(df_plot$mutation_aa, levels=unique(df_plot$mutation_aa))
df_plot$`mut (nuc)` <- factor(df_plot$`mut (nuc)`, levels=unique(df_plot$`mut (nuc)`))

func_Y_to_1 <- function(x) {
	x <- ifelse(x=="Y", 1, 0)
	x[is.na(x)] <- 0
	return(x)
}
df_plot <- df_plot %>% mutate_at(vars(`BA.2`:`BA.1`), func_Y_to_1)

df_plot <- df_plot %>% pivot_longer(cols=c(`BA.1`, `BA.2`, `Alt_freq`))
df_plot <- df_plot %>% filter(!is.na(value))
df_plot$name[df_plot$name =="Alt_freq"] <- df_plot$sample[df_plot$name =="Alt_freq"]

df_plot$gene_facet <- factor(df_plot$gene, unique(c("nuc", df_plot$gene)))
levels(df_plot$gene_facet)[1] <- "Not in ORF"

## after checking, exclude the samples with possible lab contamination
muts_BA1 <- df_plot %>% filter(name=="BA.1" & value==1) %>% .$`mut (nuc)` %>% unique() %>% as.character()
muts_BA2 <- df_plot %>% filter(name=="BA.2" & value==1) %>% .$`mut (nuc)` %>% unique() %>% as.character()

df_plot$name[df_plot$name %in% c("BA.1", "BA.2")] <- "BA.1/BA.2\ndefining mutations"

p1 <- df_plot %>% ggplot() +
	geom_tile(aes(x=name, y=`mut (nuc)`, fill=value), color="black", data=. %>% filter(`mut (nuc)`%in%muts_BA1))+
	scale_y_discrete(drop = FALSE)+
	scale_fill_gradient2(name="Mutation frequency (BA.1)", low="#FFFFFF00", mid="#F4001B50", high="#F4001B", midpoint=0.3)+ # BA.1
	new_scale_fill()+
	geom_tile(aes(x=name, y=`mut (nuc)`, fill=value), color="black", data=. %>% filter(`mut (nuc)`%in%muts_BA2))+
	scale_fill_gradient2(name="Mutation frequency (BA.2)", low="#FFFFFF00", mid="#88CD4650", high="#88CD46", midpoint=0.3)+ # BA.2
	# scale_fill_viridis_d(name="Sample mutation frequency", na.value = "white")+
	geom_tile(aes(x=name, y=`mut (nuc)`), fill="white", color="black", data=. %>% filter(grepl("WHP", name) & round(value, 2)==0))+
	geom_text(aes(x=name, y=`mut (nuc)`, label=round(value, 2)), size=3, color="black", data=. %>% filter(grepl("WHP", name)))+
	theme_classic()+
	xlab("")+
	# scale_x_discrete(guide = guide_axis(n.dodge=2))+
	ylab("Mutations (Nt)")+
	# facet_grid(rows=vars(gene_facet), scales="free", space="free")+
	theme(legend.position='top', axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
	NULL

colors_gene <- pal_jco()(length(levels(df_plot$gene_facet)))
names(colors_gene) <- levels(df_plot$gene_facet)
colors_gene[colors_gene=="#868686FF"] <- colors_gene[names(colors_gene)=="Not in ORF"]
colors_gene[names(colors_gene)=="Not in ORF"] <- "grey40"

p2 <- df_plot %>% filter(name=="BA.1/BA.2\ndefining mutations") %>% mutate(name="Gene") %>% ggplot() +
	geom_tile(aes(x=name, y=mutation_aa, fill=gene_facet), color="black")+
	scale_fill_manual(name="Gene", values=colors_gene, na.value = "white")+
	theme_classic()+
	xlab("")+
	ylab("Mutations (AA)")+
	scale_y_discrete(position = "right")+
	theme(
		# axis.text.y=element_text(colour=ifelse(levels(df_plot$`mut (nuc)`) %in% c('A20055G', 'C21618T'), 'red', 'black'),
		panel.background = element_rect(fill='transparent'),
		plot.background = element_rect(fill='transparent', color=NA),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank()
	)+
	NULL

p1_2 <- p1+p2+plot_layout(widths = c(50, 0.1), guides="keep")
ggsave("../results/plot_mut_freq_cotrols.pdf",  width = 30, height=7.5, plot=p1_2)
# save_pptx("../results/plot_mut_freq_coinfection.pptx", width = 6.5, height=7.5, plot=p1_2)
