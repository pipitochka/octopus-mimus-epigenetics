#!/bin/bash
cd ~/genome_project/downstream_analysis

echo "1. Сливаем перекрывающиеся регионы Zhunt в уникальные куски..."
sort -k1,1 -k2,2n mimus_zhunt.bed | bedtools merge -i - > mimus_zhunt_merged.bed
sort -k1,1 -k2,2n rubescens_zhunt.bed | bedtools merge -i - > rubescens_zhunt_merged.bed

echo "2. Извлекаем размеры хромосом..."
awk '/^>/{if (seqlen){print chrom"\t"seqlen}; chrom=substr($1,2); seqlen=0; next} {seqlen+=length($0)} END{print chrom"\t"seqlen}' mimus_clean.fna > mimus.chrom.sizes
awk '/^>/{if (seqlen){print chrom"\t"seqlen}; chrom=substr($1,2); seqlen=0; next} {seqlen+=length($0)} END{print chrom"\t"seqlen}' rubescens_clean.fna > rubescens.chrom.sizes

echo "3. Готовим регионы генов Mimus..."
awk '$3=="gene" {print $1"\t"($4-1)"\t"$5"\t"$10"\t.\t"$7}' mimus_clean.gtf | tr -d '";' | sort -k1,1 -k2,2n > mimus_genes.bed
awk '$3=="CDS" || $3=="exon" {print $1"\t"($4-1)"\t"$5"\t"$10"\t.\t"$7}' mimus_clean.gtf | tr -d '";' | sort -k1,1 -k2,2n | bedtools merge -i - > mimus_exons.bed
bedtools subtract -a mimus_genes.bed -b mimus_exons.bed > mimus_introns.bed
bedtools flank -i mimus_genes.bed -g mimus.chrom.sizes -l 1000 -r 0 -s > mimus_promoters.bed
bedtools flank -i mimus_genes.bed -g mimus.chrom.sizes -l 0 -r 200 -s > mimus_downstream.bed
cat mimus_genes.bed mimus_promoters.bed mimus_downstream.bed | sort -k1,1 -k2,2n | bedtools merge -i - > mimus_all.bed
bedtools complement -i mimus_all.bed -g mimus.chrom.sizes > mimus_intergenic.bed

echo "4. Готовим регионы генов Rubescens..."
awk '$3=="gene" {print $1"\t"($4-1)"\t"$5"\t"$10"\t.\t"$7}' rubescens_clean.gtf | tr -d '";' | sort -k1,1 -k2,2n > rubescens_genes.bed
awk '$3=="CDS" || $3=="exon" {print $1"\t"($4-1)"\t"$5"\t"$10"\t.\t"$7}' rubescens_clean.gtf | tr -d '";' | sort -k1,1 -k2,2n | bedtools merge -i - > rubescens_exons.bed
bedtools subtract -a rubescens_genes.bed -b rubescens_exons.bed > rubescens_introns.bed
bedtools flank -i rubescens_genes.bed -g rubescens.chrom.sizes -l 1000 -r 0 -s > rubescens_promoters.bed
bedtools flank -i rubescens_genes.bed -g rubescens.chrom.sizes -l 0 -r 200 -s > rubescens_downstream.bed
cat rubescens_genes.bed rubescens_promoters.bed rubescens_downstream.bed | sort -k1,1 -k2,2n | bedtools merge -i - > rubescens_all.bed
bedtools complement -i rubescens_all.bed -g rubescens.chrom.sizes > rubescens_intergenic.bed

echo " "
echo "============================================="
echo "ТАКСИ ДЛЯ ПРЕЗЕНТАЦИИ (Octopus Mimus):"
echo "============================================="
echo "--- G-КВАДРУПЛЕКСЫ ---"
echo "Exons: $(bedtools intersect -u -a mimus_g4.bed -b mimus_exons.bed | wc -l)"
echo "Introns: $(bedtools intersect -u -a mimus_g4.bed -b mimus_introns.bed | wc -l)"
echo "Promoters: $(bedtools intersect -u -a mimus_g4.bed -b mimus_promoters.bed | wc -l)"
echo "Downstream: $(bedtools intersect -u -a mimus_g4.bed -b mimus_downstream.bed | wc -l)"
echo "Intergenic: $(bedtools intersect -u -a mimus_g4.bed -b mimus_intergenic.bed | wc -l)"

echo " "
echo "--- Z-DNA (Алгоритм Zhunt) ---"
echo "Exons: $(bedtools intersect -u -a mimus_zhunt_merged.bed -b mimus_exons.bed | wc -l)"
echo "Introns: $(bedtools intersect -u -a mimus_zhunt_merged.bed -b mimus_introns.bed | wc -l)"
echo "Promoters: $(bedtools intersect -u -a mimus_zhunt_merged.bed -b mimus_promoters.bed | wc -l)"
echo "Downstream: $(bedtools intersect -u -a mimus_zhunt_merged.bed -b mimus_downstream.bed | wc -l)"
echo "Intergenic: $(bedtools intersect -u -a mimus_zhunt_merged.bed -b mimus_intergenic.bed | wc -l)"

echo " "
echo "============================================="
echo "ТАКСИ ДЛЯ ПРЕЗЕНТАЦИИ (Octopus Rubescens):"
echo "============================================="
echo "--- G-КВАДРУПЛЕКСЫ ---"
echo "Exons: $(bedtools intersect -u -a rubescens_g4.bed -b rubescens_exons.bed | wc -l)"
echo "Introns: $(bedtools intersect -u -a rubescens_g4.bed -b rubescens_introns.bed | wc -l)"
echo "Promoters: $(bedtools intersect -u -a rubescens_g4.bed -b rubescens_promoters.bed | wc -l)"
echo "Downstream: $(bedtools intersect -u -a rubescens_g4.bed -b rubescens_downstream.bed | wc -l)"
echo "Intergenic: $(bedtools intersect -u -a rubescens_g4.bed -b rubescens_intergenic.bed | wc -l)"

echo " "
echo "--- Z-DNA (Алгоритм Zhunt) ---"
echo "Exons: $(bedtools intersect -u -a rubescens_zhunt_merged.bed -b rubescens_exons.bed | wc -l)"
echo "Introns: $(bedtools intersect -u -a rubescens_zhunt_merged.bed -b rubescens_introns.bed | wc -l)"
echo "Promoters: $(bedtools intersect -u -a rubescens_zhunt_merged.bed -b rubescens_promoters.bed | wc -l)"
echo "Downstream: $(bedtools intersect -u -a rubescens_zhunt_merged.bed -b rubescens_downstream.bed | wc -l)"
echo "Intergenic: $(bedtools intersect -u -a rubescens_zhunt_merged.bed -b rubescens_intergenic.bed | wc -l)"
