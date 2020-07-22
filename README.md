# ROC

Generates ROC curves from bedfiles (.bed)

Usage: ROCCurve.py

For example, 2 experimental groups called "group1" and "group2"

Required folder structure:

True_Positives, False_Positives, True_Negatives, False_Negatives
--> group1_1.bed, group1_2.bed, group1_(n).bed
--> group2_1.bed, group2_2.bed, group2_(n).bed

Test_data
--> group1.bed
--> group2.bed

Currently supported .bed file formats:
- ChIP-R
- IDR

Warnings:
1. Currently, the calculation at the base pair level is quite slow. Each of my experimental runs (.bed) are ~10mb in size and the runtime all up was ~60 minutes
2. Prefix names need to be completely unique (not contain eachother)  
