library(tidyverse)
library(parallel)
library(readxl)
library(writexl)
library(lubridate)
library(ggsci)
source("https://raw.githubusercontent.com/Koohoko/Save-ggplot-to-pptx/main/scripts/save_pptx.r")

# examine the samples collected after 2021-12-01
df_metadata_hk <- read_excel("../../../../2021/2021-06-24_merge_metadata/results/cleaned_metadata.xlsx", guess_max=10000)
df_metadata_hk %>% filter(grepl("^import", tolower(Classification))) %>% filter(grepl("BA.[12]", lineage)) %>% .$`Report date` %>% unique() %>% sort()
df_metadata_hk %>% filter(grepl("^import", tolower(Classification))) %>% filter(`Report date`>="2021-12-01" & `Report date`<="2022-02-04") %>% nrow()

df_metadata_hk_import_all <- df_metadata_hk %>% filter(grepl("^import", tolower(Classification))) %>% filter(`Report date`>="2021-11-15" & `Report date`<="2022-02-04") %>% filter(sequenced_by_us)
# df_metadata_hk_import %>% group_by(lineage) %>% summarise(mean_cov=mean(coverage)) %>% arrange(mean_cov)
sort(unique(df_metadata_hk_import_all$`Country of importation`))
df_metadata_hk_import_all$`Country of importation` <- gsub("/.+", "", df_metadata_hk_import_all$`Country of importation`)
df_metadata_hk_import_all$`Country of importation` <- stringr::str_to_title(df_metadata_hk_import_all$`Country of importation`)
country_table <- sort(table(df_metadata_hk_import_all$`Country of importation`), decreasing=T)
df_country_summary <- df_metadata_hk_import_all %>% group_by(`Country of importation`) %>% summarise(N=n(), Date_range=paste0(min(`Report date`), " to ", max(`Report date`))) %>% arrange(desc(N))
df_country_summary$Date_range <- sapply(df_country_summary$Date_range, function(x) {
	# print(x)
	tmp <- strsplit(x, " to ", fixed=T)[[1]]
	if(tmp[1]==tmp[2]){
		return(tmp[1])
	}else{
		return(x)
	}
})
write_csv(df_country_summary, "../results/country_of_imported_cases_all.csv")

df_metadata_hk_import_all <- df_metadata_hk %>% filter(grepl("^import", tolower(Classification))) %>% filter(`Report date`>="2021-11-15" & `Report date`<="2022-02-04") %>% filter(sequenced_by_us)
df_metadata_hk_import <- df_metadata_hk_import_all %>% filter(coverage>=0.85)

df_metadata_hk_import$lineage_sim <- df_metadata_hk_import$lineage
sort(unique(df_metadata_hk_import$lineage_sim))
df_metadata_hk_import$lineage_sim[grepl("^AY", df_metadata_hk_import$lineage_sim)] <- "Delta"
df_metadata_hk_import$lineage_sim[grepl("^B.1.617.2", df_metadata_hk_import$lineage_sim)] <- "Delta"
df_metadata_hk_import$lineage_sim[grepl("^BA.1$", df_metadata_hk_import$lineage_sim)] <- "Omicron BA.1"
df_metadata_hk_import$lineage_sim[grepl("^BA.1.1$", df_metadata_hk_import$lineage_sim)] <- "Omicron BA.1.1"
df_metadata_hk_import$lineage_sim[grepl("^BA.2", df_metadata_hk_import$lineage_sim)] <- "Omicron BA.2"

(lineage_summary <- sort(table(df_metadata_hk_import$lineage_sim), decreasing=T))
df_lineage_summary <- df_metadata_hk_import %>% group_by(lineage_sim) %>% summarise(N=n(), Date_range=paste0(min(`Report date`), " to ", max(`Report date`))) %>% arrange(desc(N))
df_lineage_summary$Date_range <- sapply(df_lineage_summary$Date_range, function(x) {
	# print(x)
	tmp <- strsplit(x, " to ", fixed=T)[[1]]
	if(tmp[1]==tmp[2]){
		return(tmp[1])
	}else{
		return(x)
	}
})
df_lineage_summary$lineage_sim[df_lineage_summary$lineage_sim=="None"] <- "None (Recombinant)"
names(df_lineage_summary)[1] <- "lineage"
write_csv(df_lineage_summary, "../results/lineage_of_imported_cases.csv")

lineage_table <- sort(table(df_metadata_hk_import$lineage_sim), decreasing=T)
df_metadata_hk_import$lineage_sim <- factor(df_metadata_hk_import$lineage_sim, levels=names(lineage_table), labels=paste0(names(lineage_table), " (N=", lineage_table, ")"))
levels(df_metadata_hk_import$lineage_sim) <- gsub("None", "None (Recombinant)", levels(df_metadata_hk_import$lineage_sim))

sort(unique(df_metadata_hk_import$`Country of importation`))
df_metadata_hk_import$`Country of importation` <- gsub("/.+", "", df_metadata_hk_import$`Country of importation`)
df_metadata_hk_import$`Country of importation` <- stringr::str_to_title(df_metadata_hk_import$`Country of importation`)
country_table <- sort(table(df_metadata_hk_import$`Country of importation`), decreasing=T)
df_country_summary <- df_metadata_hk_import %>% group_by(`Country of importation`) %>% summarise(N=n(), Date_range=paste0(min(`Report date`), " to ", max(`Report date`))) %>% arrange(desc(N))
df_country_summary$Date_range <- sapply(df_country_summary$Date_range, function(x) {
	# print(x)
	tmp <- strsplit(x, " to ", fixed=T)[[1]]
	if(tmp[1]==tmp[2]){
		return(tmp[1])
	}else{
		return(x)
	}
})
write_csv(df_country_summary, "../results/country_of_imported_cases.csv")

country_keep <- country_table[country_table>=10]
country_others <- country_table[country_table<10]
df_metadata_hk_import$`Country of importation`[df_metadata_hk_import$`Country of importation` %in% names(country_others)] <- "Others"

df_metadata_hk_import$`Country of importation` <- factor(df_metadata_hk_import$`Country of importation`, levels=c(names(country_keep), "Others"), labels=c(paste0(names(country_keep), " (N=", country_keep, ")"), "Others (N<10)"))

p1 <- df_metadata_hk_import %>% ggplot() +
	geom_histogram(aes(x=ymd(`Report date`), fill=`Country of importation`), stat="count", color="black", size=0.2)+
	scale_x_date(date_breaks = "1 week", date_labels = "%b-%d", expand = c(0, 0), limits = c(ymd("2021-11-10"), NA), guide = guide_axis(n.dodge = 2))+
	scale_fill_jama(name="Country")+
	facet_wrap(vars(lineage_sim), ncol=1)+
	xlab("Report date")+
	ylab("Number of cases")+
	theme(legend.position="bottom")+
	NULL

ggsave("../results/time_series_imported_cases.pdf", width=7.5/sqrt(2)*1.2, height=7.5*1.2, plot=p1)
save_pptx("../results/time_series_imported_cases.pptx", width=7.5/sqrt(2), height=7.5, plot=p1)
write_xlsx(df_metadata_hk_import, "../results/df_metadata_hk_import.xlsx")

df_metadata_hk_BA12 <- df_metadata_hk_import %>% filter(grepl("BA.[12]", lineage)) 
samples_tocheck <- gsub("_\\d+", "", df_metadata_hk_BA12$Sample)
samples_tocheck <- gsub("-.+", "", samples_tocheck)
samples_tocheck <- unique(samples_tocheck)

files_bamstat <- list.files("../../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/pysamstats/", full.names=T)
files_bamstat_sim <- list.files("../../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/pysamstats/")
files_bamstat_sim <- gsub("-.+", "", files_bamstat_sim)
files_bamstat_sim <- gsub(".tsv", "", files_bamstat_sim)
samples_tocheck[!samples_tocheck %in% files_bamstat_sim]
files_bamstat_tocheck <- files_bamstat[files_bamstat_sim %in% samples_tocheck]

df_ba12 <- read_csv("../data/Omicron_BA.1_BA.2_mutations.csv")
split_tmp <- strsplit(df_ba12$`mut (nuc)`, "[,] ")
df_vcf_int <- tibble(mutation_aa=rep(df_ba12$`mutation (aa)`, sapply(split_tmp, length)),`mut (nuc)`=rep(df_ba12$`mut (nuc)`, sapply(split_tmp, length)), mutation_nt=unlist(split_tmp))
df_vcf_int <- left_join(df_vcf_int, df_ba12 %>% mutate(mutation_aa=`mutation (aa)`) %>% select(mutation_aa, gene, `BA.2`:notes))

df_vcf_int$POS <- gsub("[_ ].+$", "", df_vcf_int$`mutation_nt`)
df_vcf_int$mutation_nt_sim <- gsub(" \\(.+", "", df_vcf_int$`mutation_nt`)
df_vcf_int$POS <- gsub("\\D", "", df_vcf_int$POS)
df_vcf_int$POS <- as.numeric(df_vcf_int$POS)
df_vcf_int$REF <- sapply(strsplit(df_vcf_int$`mutation_nt_sim`, "\\d"), function(x){x[1]})
df_vcf_int$ALT <- sapply(strsplit(df_vcf_int$`mutation_nt_sim`, "\\d"), function(x){x[length(x)]})
df_vcf_int$ALT[grepl("_", df_vcf_int$`mutation_nt_sim`)|grepl("del", tolower(df_vcf_int$`mutation_nt_sim`))] <- "deletions"
df_vcf_int$ALT[grepl("ins", df_vcf_int$`mutation_nt_sim`)] <- "insertions"

df_vcf_int <- df_vcf_int[-c(19,20), ] # df_vcf_int$notes
writexl::write_xlsx(df_vcf_int, "../data/Omicron_BA.1_BA.2_mutations_cleaned.xlsx")

source("./helper/af_by_bamstat.r")
df_putative_coinfection <- mclapply(files_bamstat_tocheck, function(x) {
	print(x)
	sample_x <- gsub(".+\\/", "", x)
	sample_x <- gsub(".tsv", "", sample_x)
	df_bamstat <- read_tsv(x, show_col_types = FALSE)
	tmp <- af_by_bamstat(df_bamstat, df_vcf_int$POS, df_vcf_int$REF, df_vcf_int$ALT)
	tmp$Alt_freq <- tmp$ALT_depth/tmp$Total_depth
	tmp$MAF <- cal_maf(df_bamstat, df_vcf_int$POS)

	df_input_check_valid <- df_vcf_int %>% select(`mut (nuc)`, `BA.2`:`B.1.1.529`)
	df_tmp <- filter_valid(df=df_input_check_valid, af_tmp=tmp, variant="all", min_depth=10, min_AF=0.1)
	n_ba1 <- check_num_valid(df_tmp=df_tmp, variant="BA.1", min_depth=10, min_AF=0.1)
	n_ba2 <- check_num_valid(df_tmp=df_tmp, variant="BA.2", min_depth=10, min_AF=0.1)
	n_b11529 <- check_num_valid(df_tmp=df_tmp, variant="B.1.1.529", min_depth=10, min_AF=0.1)
	
	if(n_ba1>=3 & n_ba2>=3){ # check coninfection
		return(tibble(sample_x, n_ba1, n_ba2, n_b11529, df_freq=list(df_tmp)))
	} else {
		return(NA)
	}
}, mc.cores=4)

df_putative_coinfection <- bind_rows(df_putative_coinfection[!sapply(df_putative_coinfection, function(x) {length(x)==1})])
names(df_putative_coinfection)[1:4] <- c("sample", "n_BA_1", 'n_BA_2', "n_B_1_1_529")
df_putative_coinfection <- df_putative_coinfection %>% filter(n_B_1_1_529>=30) # at least 30 Omicron defining mutations
df_putative_coinfection[,1:4] %>% write_csv("../results/summary_potential_coinfection.csv")

df_mut_freq <- do.call(c, df_putative_coinfection[,5])
df_mut_freq <- bind_rows(df_mut_freq)
df_mut_freq$sample <- rep(df_putative_coinfection$sample, sapply(df_putative_coinfection[[5]], nrow))
write_csv(df_mut_freq, "../results/mut_freq_of_potential_coinfection.csv")

# control samples
samples_controls <- readLines("../data/novaseq_batch_samples.txt")
# samples_controls <- samples_controls[!samples_controls %in% df_putative_coinfection$sample]
# samples_bam <- (gsub(".+\\/", "", files_bamstat_tocheck) %>% gsub(pattern=".tsv", replacement=""))
files_bamstat_control <- files_bamstat_tocheck[samples_bam %in% samples_controls]

df_controls <- mclapply(files_bamstat_control, function(x) {
	print(x)
	sample_x <- gsub(".+\\/", "", x)
	sample_x <- gsub(".tsv", "", sample_x)
	df_bamstat <- read_tsv(x, show_col_types = FALSE)
	tmp <- af_by_bamstat(df_bamstat, df_vcf_int$POS, df_vcf_int$REF, df_vcf_int$ALT)
	tmp$Alt_freq <- tmp$ALT_depth/tmp$Total_depth
	tmp$MAF <- cal_maf(df_bamstat, df_vcf_int$POS)

	df_input_check_valid <- df_vcf_int %>% select(`mut (nuc)`, `BA.2`:`B.1.1.529`)
	df_tmp <- filter_valid(df=df_input_check_valid, af_tmp=tmp, variant="all", min_depth=10, min_AF=0.1)
	n_ba1 <- check_num_valid(df_tmp=df_tmp, variant="BA.1", min_depth=10, min_AF=0.1)
	n_ba2 <- check_num_valid(df_tmp=df_tmp, variant="BA.2", min_depth=10, min_AF=0.1)
	n_b11529 <- check_num_valid(df_tmp=df_tmp, variant="B.1.1.529", min_depth=10, min_AF=0.1)
	
	return(tibble(sample_x, n_ba1, n_ba2, n_b11529, df_freq=list(df_tmp)))
}, mc.cores=4)

df_controls <- bind_rows(df_controls[!sapply(df_controls, function(x) {length(x)==1})])
df_mut_freq_controls <- do.call(c, df_controls[,5])
df_mut_freq_controls <- bind_rows(df_mut_freq_controls)
df_mut_freq_controls$sample <- rep(df_controls$sample_x, sapply(df_controls[[5]], nrow))
write_csv(df_mut_freq_controls, "../results/mut_freq_of_controls.csv")
