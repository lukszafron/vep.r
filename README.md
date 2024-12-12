# $${\color{lightgreen}Vep.r}$$

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
