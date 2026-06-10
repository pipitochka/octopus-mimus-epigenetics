#!/bin/bash
cd ~/genome_project/downstream_analysis

echo "1. Сливаем перекрывающиеся регионы ZDNABERT в уникальные куски..."
sort -k1,1 -k2,2n mimus_zdnabert.bed | bedtools merge -i - > mimus_zdnabert_merged.bed
sort -k1,1 -k2,2n rubescens_zdnabert.bed | bedtools merge -i - > rubescens_zdnabert_merged.bed

echo " "
echo "============================================="
echo "ФИНАЛЬНАЯ ТАБЛИЦА ZDNABERT (Octopus Mimus):"
echo "============================================="
echo "Exons: $(bedtools intersect -u -a mimus_zdnabert_merged.bed -b mimus_exons.bed | wc -l)"
echo "Introns: $(bedtools intersect -u -a mimus_zdnabert_merged.bed -b mimus_introns.bed | wc -l)"
echo "Promoters: $(bedtools intersect -u -a mimus_zdnabert_merged.bed -b mimus_promoters.bed | wc -l)"
echo "Downstream: $(bedtools intersect -u -a mimus_zdnabert_merged.bed -b mimus_downstream.bed | wc -l)"
echo "Intergenic: $(bedtools intersect -u -a mimus_zdnabert_merged.bed -b mimus_intergenic.bed | wc -l)"

echo " "
echo "============================================="
echo "ФИНАЛЬНАЯ ТАБЛИЦА ZDNABERT (Octopus Rubescens):"
echo "============================================="
echo "Exons: $(bedtools intersect -u -a rubescens_zdnabert_merged.bed -b rubescens_exons.bed | wc -l)"
echo "Introns: $(bedtools intersect -u -a rubescens_zdnabert_merged.bed -b rubescens_introns.bed | wc -l)"
echo "Promoters: $(bedtools intersect -u -a rubescens_zdnabert_merged.bed -b rubescens_promoters.bed | wc -l)"
echo "Downstream: $(bedtools intersect -u -a rubescens_zdnabert_merged.bed -b rubescens_downstream.bed | wc -l)"
echo "Intergenic: $(bedtools intersect -u -a rubescens_zdnabert_merged.bed -b rubescens_intergenic.bed | wc -l)"
echo "============================================="
