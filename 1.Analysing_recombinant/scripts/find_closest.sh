#!/bin/bash

# pre-requisite: seqtk, gofasta, figleaf

eval "$(conda shell.bash hook)"
conda activate base 

# for GISAID seqs
for sample_t in "representative" "representative_p1" "representative_p2"
do
	{
		# build bed/region file for N region
		sed '$ s/$/A/' ../results/$sample_t".fasta"  | perl -ne 'chomp;if( />(.*)/){$head = $1; $i=0; next};@a=split("",$_); foreach(@a){$i++; if($_ eq "N" && $s ==0 ){$z=$i-1; print "$head\t$z"; $s =1}elsif($s==1 && $_ ne "N"){$j=$i-1;print "\t$j\n";$s=0}}' - | awk 'BEGIN {OFS="\t"} {print $2, $3}' > ../results/$sample_t"_mask.range"

		# mask reference seq
		figleaf -fi ../results/$sample_t".fasta" -r ../results/$sample_t"_mask.range" -fo ../results/$sample_t"_masked.fasta" --hard_mask_letter "?"

		# mask database seqs
		figleaf -fi /Volumes/SSD_480G/Downloads/GISAID_BA1_2_bg.aln.fasta.gz -r ../results/$sample_t"_mask.range" -fo /Volumes/SSD_480G/Downloads/GISAID_BA1_2_bg.aln.$sample_t"_mask.fasta.gz" --hard_mask_letter "?"

		# filter the masked seqs, keeping those with no "N"s on mut_sites of representative, and N proportion <=1%
		# python3 ./helper/sequence_cleaner.py /Volumes/SSD_480G/Downloads/test.fasta 2 "../results/representative_mutsites.txt" | grep ">" # test the function
		python3 ./helper/sequence_cleaner.py /Volumes/SSD_480G/Downloads/GISAID_BA1_2_bg.aln.$sample_t"_mask.fasta.gz" 1 "../results/representative_mutsites.txt" > /Volumes/SSD_480G/Downloads/GISAID_BA1_2_bg.aln.$sample_t"_mask_filter.fasta"

		# find closest against the filtered database seqs
		gofasta closest -t 8 -n 10000 --query ../results/$sample_t"_masked.fasta" --target /Volumes/SSD_480G/Downloads/GISAID_BA1_2_bg.aln.$sample_t"_mask_filter.fasta" -o ../results/$sample_t"_GISAID_closest_n10000.csv"

		# filter closest fasta
		sed -n '2p' ../results/$sample_t"_GISAID_closest_n10000.csv" | sed 's/[;,]/\n/g' | seqtk subseq /Volumes/SSD_480G/Downloads/GISAID_BA1_2_bg.aln.$sample_t"_mask_filter.fasta" - > ../results/$sample_t"_GISAID_closest_n10000_maskseq.fasta"
		sed -n '2p' ../results/$sample_t"_GISAID_closest_n10000.csv" | sed 's/[;,]/\n/g' | seqtk subseq /Volumes/SSD_480G/Downloads/GISAID_BA1_2_bg.aln.fasta.gz - > ../results/$sample_t"_GISAID_closest_n10000_alnseq.fasta"

		# filter closest metadata
		grep -e '^>' ../results/$sample_t"_GISAID_closest_n10000_alnseq.fasta" | sed 's/>//g' > ../results/header_$sample_t".gisaid.tmp"
		gzcat /Volumes/SSD_480G/Downloads/GISAID_BA1_2_metadata.tsv.gz | awk -F "\t" 'FNR==NR{a[$1];next} {if($1 in a || FNR==1) print $0}' ../results/header_$sample_t".gisaid.tmp" - > ../results/$sample_t"_GISAID_closest_n10000_metadata.tsv"

		# # manual check if the metadata is corresponding to the fasta sequences 
		# grep -e '^>' ../results/$sample_t"_GISAID_closest_n10000_alnseq.fasta" | awk 'BEGIN{FS="/"}{print $2}' | sort | uniq -c | sort
		# head -1 ../results/$sample_t"_GISAID_closest_n10000_metadata.tsv" | sed 's/\t/\n/g' | cat -n
		# awk -F "\t" '{print $10}' ../results/$sample_t"_GISAID_closest_n10000_metadata.tsv" | awk 'BEGIN{FS=" / "}{print $2}' | sort | uniq -c | sort
		# wc -l ../results/$sample_t"_GISAID_closest_n10000_metadata.tsv"
		# grep -e '^>' ../results/$sample_t"_GISAID_closest_n10000_alnseq.fasta" | wc -l

		rm /Volumes/SSD_480G/Downloads/GISAID_BA1_2_bg.aln.$sample_t"_mask_filter.fasta"
		rm ../results/header_$sample_t".gisaid.tmp"

	} &
done
wait 

# for Genbank seqs
for sample_t in "representative" "representative_p1" "representative_p2"
do
	{
		# # build bed/region file for N region
		# sed '$ s/$/A/' ../results/$sample_t".fasta"  | perl -ne 'chomp;if( />(.*)/){$head = $1; $i=0; next};@a=split("",$_); foreach(@a){$i++; if($_ eq "N" && $s ==0 ){$z=$i-1; print "$head\t$z"; $s =1}elsif($s==1 && $_ ne "N"){$j=$i-1;print "\t$j\n";$s=0}}' - | awk 'BEGIN {OFS="\t"} {print $2, $3}' > ../results/$sample_t"_mask.range"

		# # mask reference seq
		# figleaf -fi ../results/$sample_t".fasta" -r ../results/$sample_t"_mask.range" -fo ../results/$sample_t"_masked.fasta" --task exclude

		# mask database seqs
		figleaf -fi /Volumes/SSD_480G/Downloads/Genbank_BA1_2.fasta.gz -r ../results/$sample_t"_mask.range" -fo /Volumes/SSD_480G/Downloads/Genbank_BA1_2.$sample_t"_mask.fasta.gz" --hard_mask_letter "?"

		# filter the masked seqs, keeping those with no "N"s on mut_sites of representative, and N proportion <=1%
		python3 ./helper/sequence_cleaner.py /Volumes/SSD_480G/Downloads/Genbank_BA1_2.$sample_t"_mask.fasta.gz" 1 "../results/representative_mutsites.txt" > /Volumes/SSD_480G/Downloads/Genbank_BA1_2.$sample_t"_mask_filter.fasta"

		# find closest against the filtered database seqs
		gofasta closest -t 8 -n 10000 --query ../results/$sample_t"_masked.fasta" --target /Volumes/SSD_480G/Downloads/Genbank_BA1_2.$sample_t"_mask_filter.fasta" -o ../results/$sample_t"_Genbank_closest_n10000.csv"
		
		# filter closest fasta
		sed -n '2p' ../results/$sample_t"_Genbank_closest_n10000.csv" | sed 's/[;,]/\n/g' | seqtk subseq /Volumes/SSD_480G/Downloads/Genbank_BA1_2.$sample_t"_mask_filter.fasta" - > ../results/$sample_t"_Genbank_closest_n10000_maskseq.fasta"
		sed -n '2p' ../results/$sample_t"_Genbank_closest_n10000.csv" | sed 's/[;,]/\n/g' | seqtk subseq /Volumes/SSD_480G/Downloads/Genbank_BA1_2.fasta.gz - > ../results/$sample_t"_Genbank_closest_n10000_alnseq.fasta"

		# filter closest metadata
		grep -e '^>' ../results/$sample_t"_Genbank_closest_n10000_alnseq.fasta" | sed 's/>//g' > ../results/header_$sample_t".genbank.tmp"
		gzcat /Volumes/SSD_480G/Downloads/Genbank_BA1_2_metadata.tsv.gz | awk -F "\t" 'FNR==NR{a[$1];next} {if($1 in a || FNR==1) print $0}' ../results/header_$sample_t".genbank.tmp" - > ../results/$sample_t"_Genbank_closest_n10000_metadata.tsv"
		
		# manual check if the metadata is corresponding to the fasta sequences 
		# grep -e '^>' ../results/$sample_t"_Genbank_closest_n10000_alnseq.fasta" | sort | uniq -c | sort
		# awk -F "\t" '{print $1}' ../results/$sample_t"_Genbank_closest_n10000_metadata.tsv" | sort | uniq -c | sort
		# wc -l ../results/$sample_t"_Genbank_closest_n10000_metadata.tsv"
		# grep -e '^>' ../results/$sample_t"_Genbank_closest_n10000_alnseq.fasta" | wc -l

		rm /Volumes/SSD_480G/Downloads/Genbank_BA1_2.$sample_t"_mask_filter.fasta"
		rm ../results/header_$sample_t".genbank.tmp"


	} &
done
wait


conda deactivate
