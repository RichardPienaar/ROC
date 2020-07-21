# ROC

Generates ROC curves from bedfiles (.bed)

Usage: ROCCurve.py

For example, 2 experimental groups called "group1" and "group2"

Required folder structure:

True_Positives, False_Positives, True_Negatives, False_Negatives
--> group1_1.bed, group1_2.bed, group1(n).bed
--> group2_1.bed, group2_2.bed, group2_(n).bed

Test_data
--> group1.bed
--> group2.bed

Currently supported .bed file formats:
- ChIP-R
- IDR
