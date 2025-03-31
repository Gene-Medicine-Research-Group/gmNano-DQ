# gmNano-DQ
Python script to extract DNA Phred quality scores from a fastq file

The intened use is to analyse DNA qulaity scores from a pDNA Nanopore sequencing project.
gmNano-DQ takes the consensus DNA sequence fastq file generated and outputs Phred quality scores for each base
gmNano-DQ can also calculate the average Phred quality score for the consensus sequence (which most people think isn't very useful)
and the expected/most likely number of base errors - which is more useful.
