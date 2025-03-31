# gmNano-DQ.py
# v02 2025-03
# Steve Hyde

# gmNano-DQ takes as input a consensus Nanopore plasmid sequencing fastq file which contaons ONE consensus DNA sequence with its associated Phred quality score string
# gmNano-DQ can't process as fastq file with more than one embedded DNA sequence, it tells you so and then quits 
# When gmNano-DQ outputs to the terminal it creates four columns: base#, base, Phred quality ASCII code, numerical Phred quality score
# When gmNano-DQ outputs to a tab deliminated TXT file it creates two columns: base#, numerical Phred quality score
# additional options output headers that can be used when directly importing the output file into GraphPad Prism for visualisation.

# suggested command line usage:
# python gmNano-DQ.py -h                                    <<<< provides brief help on the arguments gmNano-DQ expects
# python gmNano-DQ.py -i filename.fastq                     <<<< outputs base#, base, Phred quality score to the terminal window
# python gmNano-DQ.py -i filename.fastq -o filename.txt     <<<< outputs base#, base, Phred quality score to filename.txt
# python gmNano-DQ.py -i filename.fastq -o filename.txt -op <<<< outputs base#, base, Phred quality score to filename.txt, includes a GraphPad Prism Project Information header 

# example
# python gmNano-DQ.py -i 613162701_1346-3_barcode31.final.fastq -o myPrism.txt -op -v

# useful information on fastq files and Phred quality scores:
# https://en.wikipedia.org/wiki/FASTQ_format
# https://en.wikipedia.org/wiki/Phred_quality_score

# average Phred quality score is calculated according to the approach used by Nanoplot (convert quality scores to error probability, average, convert back)
# https://gigabaseorgigabyte.wordpress.com/2017/06/12/oxford-nanopore-basecall-quality-scores/
# https://gigabaseorgigabyte.wordpress.com/2018/08/30/update-on-oxford-nanopore-basecall-quality-scores/
# to me this seems to provide a value that is too low... may need to refelct on this.. consider using a simple average of the quality score??

# Perhaps more helpful than average Phred quality score is the expecetd and most probable number of base errors in the provided fastq file
# This can be calculated as described here:
# https://www.drive5.com/usearch/manual/exp_errs.html
# https://pubmed.ncbi.nlm.nih.gov/26139637/

# useful information on GraphPad Prism headers
# https://www.graphpad.com/guides/prism/latest/user-guide/using_paste_special_notes.htm
