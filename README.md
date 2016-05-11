# python implement fast BLAST 
Python implementation of Basic Local Alignment Search Tool (BLAST) , which is the core algorithm in sequence alignmenrt for genomes.

This is the course project for Bioinformatics(BI3204 2015.09-2016.01) at [SUSTC](http://www.sustc.edu.cn/).

**Table of Contents**

- [Introduction to BLAST](#Introduction-to-BLAST)
- [BLAST implementation in python: For human genome](#BLAST-implementation-in-python:-For-human-genome)
  - [construct library](#construct-library)
  - [alignment algrithms](#alignment-algrithms)

##Introduction to BLAST

![1](https://github.com/JiaShun-Xiao/python-implement-fast-BLAST-Basic-Local-Alignment-Search-Tool/blob/master/images/1.jpg)
![2](https://github.com/JiaShun-Xiao/python-implement-fast-BLAST-Basic-Local-Alignment-Search-Tool/blob/master/images/2.jpg =100x20)
![3](https://github.com/JiaShun-Xiao/python-implement-fast-BLAST-Basic-Local-Alignment-Search-Tool/blob/master/images/3.png)

>In bioinformatics, BLAST for Basic Local Alignment Search Tool is an algorithm for comparing primary biological sequence information, such as the amino-acid sequences of different proteins or the nucleotides of DNA sequences. A BLAST search enables a researcher to compare a query sequence with a library or database of sequences, and identify library sequences that resemble the query sequence above a certain threshold.

>BLAST.  In Wikipedia. Retrieved May 12, 2016, from https://en.wikipedia.org/wiki/BLAST


##BLAST implementation in python: For human genome
####construct library
![4](https://github.com/JiaShun-Xiao/python-implement-fast-BLAST-Basic-Local-Alignment-Search-Tool/blob/master/images/4.png)
>construct library for human genome. Break whole genome (respectively for chromosome) sequence into 11 bases length words overlappedly, every word as a name of a txt file which contain all lolation of this words in genome sequence. So was the query reads.

```bash
perl BWT-Aligner.pl lambda.tally reads_1.fq reads_1.sam
# Since the index is built, you can now align the raw reads onto the reference.
# This will output the alignment in SAM format, reads_1.sam.
```
####alignment algrithms
![5](https://github.com/JiaShun-Xiao/python-implement-fast-BLAST-Basic-Local-Alignment-Search-Tool/blob/master/images/5.png)

```bash
samtools sort -o reads_1.srt.bam -O bam -T temp reads_1.sam  # Convert it to sorted bam
 samtools index reads_1.srt.bam # Index it by samtools index
# View in IGV or other alignments viewer.
```
>File source specification:

Building BWT index for huge genome like human is a little more complicated because of the sorting process. This perl implementation uses 140bp header to sort in lambda virus reference. For human genome 1000bp has been used to sort but there are still some identical sequences. So I have sorted it step by step. I will push this part after I clean up my code.
