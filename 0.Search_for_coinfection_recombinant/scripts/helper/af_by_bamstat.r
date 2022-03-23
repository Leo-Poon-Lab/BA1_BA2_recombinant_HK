
af_by_bamstat <- function(df_bam, POS, REF, ALT) {
	df_tmp <- df_bam
	stopifnot(length(POS)==length(ALT))
	stopifnot(length(POS)==length(REF))
	tmp <- lapply(seq_along(POS), function(i){
		pos_i <- POS[i]
		ref_i <- REF[i]
		alt_i <- ALT[i]
		if(alt_i == "deletions"){
			pos_i <- pos_i + 1
			mut <- alt_i
		} else if(alt_i == "insertions"){
			mut <- alt_i
		} else{
			if(nchar(ref_i)>nchar(alt_i)){ # deletion
				pos_i <- pos_i + 1
				mut <- "deletions"
			} else if(nchar(ref_i)<nchar(alt_i)){ # insertion
				mut <- "insertions"
			} else { # sub
				mut <- toupper(alt_i)
			}
		}
		depth_i <- df_tmp$reads_all[pos_i]
		depth_alt_i <- df_tmp[[mut]][pos_i]
		return(c(depth_alt_i, depth_i))
	})
	tmp <- matrix(unlist(tmp), ncol=2, byrow=T)
	colnames(tmp) <- c("ALT_depth", "Total_depth")
	return(as_tibble(tmp))
}


cal_maf <- function(df_bam, POS){
	df_tmp <- df_bam
	idx <- which(names(df_tmp) %in% c("A", "C", "T", "G", "N", "insertions", "deletions"))
	df_tmp <- df_tmp[,idx]
	df_tmp <- df_tmp[POS,]
	maf <- apply(df_tmp,1,function(x){
		x[order(x, decreasing=T)[2]]/sum(x)
	})
	return(maf)
}

filter_valid <- function(df, af_tmp, variant="all", min_depth, min_AF) {
	df_tmp <- bind_cols(df, af_tmp)
	id_col_df <- seq_len(ncol(df))
	if(variant!="all"){
		df_tmp <- df_tmp[which(df[[variant]]=="Y"),]
	} 
	mutations_to_merge <- unique(df_tmp$`mut (nuc)`[grepl(",", df_tmp$`mut (nuc)`)])
	if(length(mutations_to_merge)>0){
		sapply(mutations_to_merge, function(y){
			idx <- which(df_tmp$`mut (nuc)`==y)
			check <- all(df_tmp$Alt_freq[idx]>=min_AF)
			if(is.na(check)){
				df_tmp[idx[1],-id_col_df] <<- NA
				df_tmp <<- df_tmp[-idx[2:length(idx)],]
				return(NA)
			}
			if(check){
				df_tmp[idx[1],-id_col_df] <<- sapply(apply(df_tmp[idx,-id_col_df],2,mean),list)
				df_tmp <<- df_tmp[-idx[2:length(idx)],]
			} else{
				df_tmp[idx[1],-id_col_df] <<- sapply(apply(df_tmp[idx,-id_col_df],2,mean),list)
				df_tmp[idx[1],"Alt_freq"] <<- 0
				df_tmp[idx[1],"MAF"] <<- 0
				df_tmp <<- df_tmp[-idx[2:length(idx)],]
			}
		})
	}
	return(df_tmp)
}

check_num_valid <- function(df, af_tmp, variant="all", min_depth, min_AF, df_tmp=NULL){
	if(is.null(df_tmp)){
		df_tmp <- filter_valid(df, af_tmp, variant, min_depth, min_AF)
	} else{
		if(variant!="all"){
			df_tmp <- df_tmp[which(df_tmp[[variant]]=="Y"),]
		}
	}
	n <- sum((df_tmp$Total_depth>=min_depth) & (df_tmp$Alt_freq>=min_AF), na.rm=T)
	return(n)
}