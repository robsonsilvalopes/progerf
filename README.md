# ProGeRF
Proteome and Genome Repeat Finder Tool

Proteome and Genome Repeat Finder Tool(ProGeRF) is an efficient, fast and accurate tool to search for repeat sequence. ProGeRF allows to perform searches to identifier Perfect and Imperfect or Simple Sequence Repeats (SSR's) or Short Tandem Repeats (STR's) from genome or proteome sequences.

Using the ProGeRF it's possible to:

Search perfect repeats in genome and proteome;
Search imperfect repeat in genome and proteome;
Define the imperfection limit for repeat unit of each size;
Define Repeats of a particular size or for sizes within a range of values;
Define the minimum number of repeat units of a tract of each size or the same value for all;
Define the percentage of overlap allowed between sequence found;

ProGeRF was developed using the dictionary approach to extract perfects repeat region in (multi)fasta file, similarity the TROLL what allowing a linear time complexity. However, a function that allows identifier repeat regions with imperfection and to define diferent overlapping degree was developed.

## Install
ProGeRF was designed to be web server tool, but it's possible to download and performed stand-alone by command line interface.

Clone this repository.

````
git clone https://github.com/robsonsilvalopes/progerf/
````

Type the follow comand.

````
#unzip progerf.zip
#cd progerf/
#./makefile
````
Primary command is to unzip the file. Second is to entry in the directory and, finally, third command is to compile algorithms in C language. Now, it's possible performed the ProGeRF.

## Command Line
````
USAGE:
   perl progerf.pl -q [FILE] -o [STRING] -i [INTERGER] -y [INTERGER] -r [INTERGER] -g [INTERGER] -v [INTERGER] -d [INTERGER]


REQUIRED ARGUMENTS:
   -q [FILE]
      Input file with DNA/Proteome sequence in (multi)FASTA format.
   -o [STRING]
      Output file to be create. default results;
   -i [INTERGER]
      Minimum length of motif. default 2;
   -y [INTEGER]
      Maximum length of motif. default 5;
   -r [INTERGER] OR -rl [INTERGER-...-INTERGER]
      Minimum repeated times of motif. default -r 5;   
   -g [INTERGER]
      Maximum allowed Gaps between motifs of a tandem repeat. default 0;
   -v [INTERGER]
      Maximum allowed overllap.  Values allowed  0 - 100. default 0;
   -d [INTERGER]
      Maximum allowed Degeneration motifs. default 0;
   -m [STRING]
      Run mode. n to nucleotide or p to protein. Default is n.
   

EXAMPLE:
   perl progerf.pl -q Linfantum_JPCM5.fasta -o output -i 2 -y 6 -r 5 -g 3 -v 1 -d 20 -m n
````


The command line interaction may be performed indicating the (multi)fasta file address containing DNA or Proteome sequence(s), the motif length range, the minimum repeated times for all motif length or the minimum repeated time for each one, the maximum gaps allowed between motif, the maximum degeneration percentage, the motif shifting percentage, the rum mode where is defined if will performed to DNA or Protein and the output file name. 

For example, the command sequence perl progerf.pl -q Linfantum_JPCM5.fasta -o output_file -i 2 -y 6 -r 5 -g 3 -v 0 -d 20 -m n will search SSRs at file Linfantum_JPCM5.fasta of motif with length range between 2 and 6, with gaps maximum of 3, motif overlapping of 0%, degeneration of 20%, rum mode nucleotide and the result will be save on Out_output_file.ssr

# OUTPUT

The output of the command line is a file with tabular data, where the columns represent the following information in this order:

- sequence id
- sequence length
- minimal repeat
- tract length
- start position
- end position
- total number of gaps
- repetitive motif
- statistic(only nucleotide)
- repeat sequence

# REFERENCES
Lopes Rda S, Moraes WJ, Rodrigues Tde S, Bartholomeu DC. ProGeRF: proteome and genome repeat finder utilizing a fast parallel hash function. Biomed Res Int. 2015;2015:394157. doi:10.1155/2015/394157

