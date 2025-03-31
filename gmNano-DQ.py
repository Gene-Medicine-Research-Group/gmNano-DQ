#!/usr/bin/python

# gmNano-DQ.py
# v02 2025-03
# Steve Hyde

# gmNano-DQ takes as input a consensus Nanopore plasmid sequencing fastq file which contaons ONE consensus DNA sequence with its associated Phred quality score string
# gmNano-DQ can't process as fastq file with more than one embedded DNA sequence, it tells you so and then quits 
# gmNano-DQ outputs either to the terminal or to a tab deliminated TXT file three columns: base#, base, Phred quality score

# suggested command line usage:
# python gmNano-DQ.py -h                                    <<<< provides brief help on the arguments gmNano-DQ expects
# python gmNano-DQ.py -i filename.fastq                     <<<< outputs base#, base, Phred quality score to the terminal window
# python gmNano-DQ.py -i filename.fastq -o filename.txt     <<<< outputs base#, base, Phred quality score to filename.txt 

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

# ther provided output include a

#needed for math.log10() and math.floor() [equivalent of rounddown() in excel, fwiw math.ceil() is equivalent of roundup()]
import math 

# setup parsed commands
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input file name to process")
parser.add_argument("-o", "--output", help="Output analysis file name to write")
parser.add_argument("-op", "--output_prism_project_info", help="Output Prism Project Info", action="store_true")
parser.add_argument("-v", "--verbose", help="Provide greater explanation", action="store_true")
args = parser.parse_args()

if args.verbose:
    print ("Verbose option selected, will provide as much information as possible")

if args.input:
    myInputFileName = args.input
    if args.verbose:
        print ('gmNano-DQ processing '+ myInputFileName)
else:
    print ("gmNano-DQ requires a fastq file. None provided. try --input filename.fastq or the simpler -i filename.fastq")
    exit (1)

if args.output:
    myOutputFileName = args.output
    if args.verbose:
        print ("Will output results to "+myOutputFileName)
    myPrintFile = False
else:
    myPrintFile = True

if args.verbose and myPrintFile:
    print ("No output filename provided, will list output in terminal")

if args.output_prism_project_info:
    myOutputPrismProject = args.output_prism_project_info
    if args.verbose:
        print ("Will output Prism Project Information")
    myOutputPrismProject = True
else:
    myOutputPrismProject = False

# check the supplied fastq file. It should have extactly 4 lines
# line 1 sequence_ID
# line 2 DNA sequence
# line 3 "+" or "+_and_some_arbitrary_information_often_the_sequence_ID_again"
# line 4 Phred quality score/base

# check its the right length
myInputFile = open(myInputFileName, "r")
myTotalSequencesCount = int (sum(1 for _ in myInputFile) / 4)
myInputFile.close()
if myTotalSequencesCount != 1:
    print ("Input file not single DNA sequence fastq format")
    exit (1)    # exit with error

# open the input file at the top
myInputFile = open(myInputFileName, "r")

# read line 1 - the sequence_id
myLine = (myInputFile.readline()) # read the sequence_id line
mySequenceName = myLine.strip() # strip the newline from the sequence_id

# read line 2 - the DNA sequence 
myLine = (myInputFile.readline()) # read the associated DNA base sequence
mySequenceString = myLine.strip()
mySequenceLength = len(myLine.strip()) # strip the final newline and determine length of DNA sequence 

# read and check line 3
myLine = (myInputFile.readline()) # in fastq this line always starts with a "+" check/test - exit if not
myCheckfastqFormat = myLine.strip()
myCheckfastqFormatLength = len(myCheckfastqFormat)
myError = False

if myCheckfastqFormatLength < 1:
    myError = True

if not myError:
    myChar = myCheckfastqFormat[:1]
    if not myChar == "+":
        myError = True

if myError:
    print ("Input file not single DNA sequence fastq format")
    myInputFile.close()     # close the open file before exit
    exit (1)                # exit with error

# read line 4 - the Phred quality score 
myLine = (myInputFile.readline()) # read the Phred quality score for each base
mySequencePhredString = myLine.strip() #strip the newline from the Phred quality string
mySequencePhredLength = len(myLine.strip())

# finished reading the inpout file so close it
myInputFile.close()

# check DNA sequence and quality score are same length
if mySequenceLength != mySequencePhredLength:
    print ("Sequence and quality records of different length - unable to process this fastq file - consider stripping gaps")
    exit (1)

# if the user wants to know, supply basic DNA sequence information to terminal
if args.verbose:
    print ("Processing DNA Sequence ID: "+mySequenceName)
    print ("DNA sequence length: "+str(mySequenceLength)+" base pairs")
    
# determine the average Phred quality score
myRunningPhredProbability = 0 # variable to add up all the Phred score probabilities

for character in mySequencePhredString:
    myBasePhredScore = ord(character) - ord ('!') # convert character to Phred score where the ASCII char ! = 0
    myBasePhredProbability = pow(10,-myBasePhredScore/10) # convert the Phred score to an error probability
    myRunningPhredProbability = myRunningPhredProbability + myBasePhredProbability

myAveragePhredProbability = myRunningPhredProbability/mySequencePhredLength # calculate the average error probability
myAveragePhredScore = -10 * math.log10 (myAveragePhredProbability) # convert the average error probability to an average Phred score

if not (myPrintFile):
    myOutputFile = open(myOutputFileName, "a")

if myOutputPrismProject and (not (myPrintFile)):    # output Prism Project Info
    myOutputFile.write ("<Info>")
    myOutputFile.write("\n")

    myOutputFile.write ("DNA Sequence Length")       # provide average Phred quality score as base -3
    myOutputFile.write("\t")
    myOutputFile.write (str(mySequenceLength))
    myOutputFile.write("\n")
    
    myOutputFile.write ("Average Phred Quality Score")       # provide average Phred quality score as base -3
    myOutputFile.write("\t")
    myOutputFile.write (str("{:.9f}".format(myAveragePhredScore)))
    myOutputFile.write("\n")

    myOutputFile.write ("Expected Number Of Nanopore Sequencing Base Errors")       # provide expected number of base errors as base -2
    myOutputFile.write("\t")
    myOutputFile.write (str("{:.9f}".format(myRunningPhredProbability)))
    myOutputFile.write("\n")

    myOutputFile.write ("Most Probable Number Of Nanopore Sequencing Base Errors")       # provide most probable number of base errors as base -1
    myOutputFile.write("\t")
    myOutputFile.write (str(math.floor(myRunningPhredProbability)))
    myOutputFile.write("\n")

    myOutputFile.write ("</Info>")
    myOutputFile.write("\n")

    myOutputFile.write ("DNA Read Length (bp)")       # provide most probable number of base errors as base -1
    myOutputFile.write("\t")
    myOutputFile.write ("Base Phred Quality Score")
    myOutputFile.write("\n")


if myOutputPrismProject:        # output to terminal
    print ("Average Phred Quality Score: "+str("{:.9f}".format(myAveragePhredScore)))
    print ("Expected number of errors: "+str("{:.9f}".format(myRunningPhredProbability)))
    print ("Most probable number of errors: "+str(math.floor(myRunningPhredProbability)))

# process each base/quality score and output to terminal/file as appropriate
for myBaseCounter in range(1, mySequenceLength+1):
    myNextBase = mySequenceString [myBaseCounter-1:myBaseCounter]
    myNextPhredChar = mySequencePhredString [myBaseCounter-1:myBaseCounter]
    myNextBasePhredScore = ord(myNextPhredChar) - ord ('!') # convert character to Phred score where the ASCII char ! = 0

    if myPrintFile:     # output to terminal
        print(str(myBaseCounter)+"\t"+myNextBase+"\t"+myNextPhredChar+"\t"+str(myNextBasePhredScore))   
    else:               # output to file
        myOutputFile.write(str(myBaseCounter))
        myOutputFile.write("\t")
        myOutputFile.write(str(myNextBasePhredScore))
        myOutputFile.write("\n")


if args.output:         # if outputing to file, we are done now, so close the file 
    myOutputFile.close()

if args.verbose:        # if user wants to know we are all done, then tell them
    print ("gmNano-DQ complete")

exit(0) # we finished, exit without error