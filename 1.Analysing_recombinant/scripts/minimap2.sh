#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate base 

# prepare GISAID data
for lineage_t in "BA_1" "BA_2"
do
	{
		tar -xvf /Volumes/SSD_480G/Downloads/lineage_${lineage_t}_tsv_2022_03_07.tar.xz -C /Volumes/SSD_480G/Downloads/ '*.tsv'
		tar -xvf /Volumes/SSD_480G/Downloads/lineage_${lineage_t}_fasta_2022_03_07.tar.xz -C /Volumes/SSD_480G/Downloads/ '*.fasta'
		
		seqtk seq -l0 /Volumes/SSD_480G/Downloads/lineage_${lineage_t}.fasta | sed -n 'n;p' > /Volumes/SSD_480G/Downloads/lineage_${lineage_t}.seq
		rm /Volumes/SSD_480G/Downloads/lineage_${lineage_t}.fasta
		# wc -l /Volumes/SSD_480G/Downloads/lineage_${lineage_t}.seq
		# wc -l /Volumes/SSD_480G/Downloads/lineage_${lineage_t}.tsv
		
		tail -n +2 /Volumes/SSD_480G/Downloads/lineage_${lineage_t}.tsv | awk 'BEGIN{FS="\t";RS="\n";OFS=""}{print ">",$1}' | paste -d "\n" - /Volumes/SSD_480G/Downloads/lineage_${lineage_t}.seq  > /Volumes/SSD_480G/Downloads/GISAID_${lineage_t}_seqs.fasta 

		rm /Volumes/SSD_480G/Downloads/lineage_${lineage_t}.seq

	} &
done
wait

cat /Volumes/SSD_480G/Downloads/GISAID_BA_1_seqs.fasta /Volumes/SSD_480G/Downloads/GISAID_BA_2_seqs.fasta | pigz --fast -p 8 > /Volumes/SSD_480G/Downloads/GISAID_BA1_2_seqs.fasta.gz

tail -n +2 /Volumes/SSD_480G/Downloads/lineage_BA_2.tsv | cat /Volumes/SSD_480G/Downloads/lineage_BA_1.tsv - | pigz --fast -p 8  > /Volumes/SSD_480G/Downloads/GISAID_BA1_2_metadata.tsv.gz

# gzcat /Volumes/SSD_480G/Downloads/GISAID_BA1_2_seqs.fasta.gz | grep "^>" | wc -l # 1222642
# gzcat /Volumes/SSD_480G/Downloads/GISAID_BA1_2_metadata.tsv.gz | wc -l # 1222643
# gzcat /Volumes/SSD_480G/Downloads/GISAID_BA1_2_seqs.fasta.gz  | sed -n 'p;n' | head
# gzcat /Volumes/SSD_480G/Downloads/GISAID_BA1_2_metadata.tsv.gz | tail -n +2 | awk 'BEGIN{FS="\t";RS="\n";OFS=""}{print ">",$1}' | head

rm /Volumes/SSD_480G/Downloads/GISAID_BA_1_seqs.fasta
rm /Volumes/SSD_480G/Downloads/GISAID_BA_2_seqs.fasta
rm /Volumes/SSD_480G/Downloads/lineage_BA_1.tsv
rm /Volumes/SSD_480G/Downloads/lineage_BA_2.tsv

## minimap2
minimap2 -t 8 -a -x asm20 --score-N=0 ../data/reference.fasta /Volumes/SSD_480G/Downloads/GISAID_BA1_2_seqs.fasta.gz > /Volumes/SSD_480G/Downloads/GISAID_BA1_2_bg.sam 

## sam to multialign
gofasta sam toMultiAlign -t 8 -s /Volumes/SSD_480G/Downloads/GISAID_BA1_2_bg.sam  -o /Volumes/SSD_480G/Downloads/GISAID_BA1_2_bg.aln.fasta 2> /dev/null
rm /Volumes/SSD_480G/Downloads/GISAID_BA1_2_bg.sam
# cat /Volumes/SSD_480G/Downloads/GISAID_BA1_2_bg.aln.fasta | grep "^>" | wc -l #1222642

pigz --fast -p8 /Volumes/SSD_480G/Downloads/GISAID_BA1_2_bg.aln.fasta > /Volumes/SSD_480G/Downloads/GISAID_BA1_2_bg.aln.fasta.gz


# genbank data (already aligned) 
## this Genbank data is downloaded from https://data.nextstrain.org/files/ncov/open/aligned.fasta.xz on 2022-03-09
# sed -n '1,10p' ../results/representative_p2_closest_n1000_aln.fasta | awk '{print length}'
# xzcat -d /Volumes/SSD_480G/Downloads/aligned.fasta.xz | head | sed -n '1,10p' | awk '{print length}'
# gzcat /Volumes/SSD_480G/Downloads/metadata.tsv.gz | head -n1 | sed 's/\t/\n/g' | grep -n lineage $pango_lineage at the 20th column
# gzcat /Volumes/SSD_480G/Downloads/metadata.tsv.gz | awk -F "\t" '{print $20}' | sort | uniq -c | tee /Volumes/SSD_480G/Downloads/Genbank.lineage 
gzcat /Volumes/SSD_480G/Downloads/metadata.tsv.gz | awk -F "\t" 'NR==1{print $0} NR!=1{if($20~/^BA\.[12]/)print $0}' | pigz -p 8 > /Volumes/SSD_480G/Downloads/Genbank_BA1_2_metadata.tsv.gz 
gzcat /Volumes/SSD_480G/Downloads/Genbank_BA1_2_metadata.tsv.gz | awk -F "\t" 'NR!=1{print $1,$20}' > ../results/Genbank_BA1_2.list
cat ../results/Genbank_BA1_2.list | awk -F " " '{print $2}' | sort | uniq -c 
xzcat /Volumes/SSD_480G/Downloads/aligned.fasta.xz | seqtk subseq - ../results/Genbank_BA1_2.list | pigz -p 8 > /Volumes/SSD_480G/Downloads/Genbank_BA1_2.fasta.gz

# gzcat /Volumes/SSD_480G/Downloads/Genbank_BA1_2.fasta.gz | grep "^>" | wc -l # 767379
# gzcat /Volumes/SSD_480G/Downloads/Genbank_BA1_2_metadata.tsv.gz | wc -l # 767399

conda deactivate




