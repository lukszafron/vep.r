# $${\color{lightgreen}Vep.r}$$

This R script was developed to extend the possibilities of the ENSEMBL Variant Effect Predictor app. In combination with the vep.comparison.r, it offers the ability to identify genetic variants with a high or moderate impact on the structure/function of the corresponding proteins. The results are supplemented with detailed statistical analyses followed by the visualization of the results, allowing for the determination of genes that are characterized by a higher frequency of genetic alterations in some groups of samples compared to others.

This program requires eleven arguments.

                                  The first one is a folder where VEP results for each sample are stored (for SNP and non-SNP results, it has to be named 'SNP' and 'NON_SNP', respectively),
                                  The second one is a destination folder,
                                  the third one is a semicolon-delimited CSV file with sample names and other necessary clinical variables,
                                  the fourth determines the sample ID variable,
                                  the fifth - grouping variable,
                                  the sixth - independent factor,
                                  the seventh - a path to the file with a gene list (one gene per line) to calculate the gene signature,
                                  the eighth is the number of threads to use,
                                  the ninth says if samples are paired,
                                  the tenth determines the column with pair indicators,
                                  while the eleventh is a logical indicator determining if the FDR correction should be applied in the TOST analysis of equivalence.
