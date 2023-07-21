BetterMetaGenomics 
==================

The objective if this Python project is to conduct metagenomic analysis on a given sample.  The project receives multiple genome sequences and reads from the sample as input. Its outcome involves determining the likely source genome for each read.  Each sample contains a specific quantity of genomes.

The program is hardset to only sample 100,000 reads. So you must make sure your reads file contain at least 100,000 reads.

The program randomly samples 200 reads and pseudo aligns it to the all the genomes. This pseudo alignment ultilzes minimizer matching. Where the reads and genomes will be broken down into minimizers and only the minimizers are stored. Genomes that pass a set match threshold of at least 18 read minimizer matches will be selected for full alignment and the other genomes will be discarded. 

The input of this project will be a zip file that contains all the genome files and one read file. 

It is important that the genome files contain the word "genome" in the file name and the read file contains the words "reads.fasta" in this exact format because this is how the program distinguishes and extracts the genome files and the read files from the zipfile. 

The output is a list of reads and the genome that it matches to in the following format:
>read_0	Genome_Number_45
>read_1	Genome_Number_65
>read_2	Genome_Number_65
>read_3	Genome_Number_50
>read_4	Genome_Number_32
...


Deliverables:
-------------

bettermetagenomics.py -- code for metagenomic analysis

predictions.csv -- a list of reads and the genome that it matches to 

predictions.zip -- zipped csv of predictions.txt


Usage
-----
The program takes in one input, the reads fasta without the genome positions 

To run the program, navigate to the project directory and run:

> python3 bettermetagenomics.py samples.zip

The program takes the following argument:

* `--samples.zip`: A zip file that contains all the genome files and one read file. It is important that the genome files contain the word "genome" in the file name and the read file contains the words "reads.fasta" in this exact format because this is how the program distinguishes and extracts the genome files and the read files from the zipfile. 

It is also important to note that when unzipping the zip file the program will output the files from the zip file in your working directory. Please delete the files if you do not want them in the directory once the program finishes running. 

Examples
--------

Here is an example of how to run the program: 

>python3 bettermetagenomics.py project4-sample.zip

Performance
-----------

Runtime is based on laptop computer with the following specifications:
* Processor: 1.1 GHz Quad-Core Intel Core i5
* RAM: 8GB
* Operating System: MacOS Catalina 

For alignment of ~5000 genomes of ~200,000 nucleotides and ~100,000 reads the runtime was around 6+ hours. 

Future Improvements
-------------------
The program is quite slow and there is a lot to improve upon. Future implementations of this project might consider the use of bloom filters in addition to minimizers.