
# Commonly used file formats for next-generation sequencing (NGS) data

## FASTA

Among the most common and simplest file formats for representing nucleotide sequences is FASTA.  Essentially, each sequence is represented by a 'header' line that begins with a '>', followed by lines containing the actual nucleotide sequence. By convention, the first 'word' in the header line is a unique identifier, which is usually as accession number. Consider this example of a FASTA-formatted nucleotide sequence:

    >LC719646.1 Influenza A virus (A/swine/Tottori/B34/2020(H1N1)) segment 8 NS1, NEP genes for nonstructural protein 1, nuclear export protein, complete cds
    ATGGAATCCAACACCATGTCAAGCTTTCAGGTAGACTGTTTTCTTTGGCATATTCGCAAGCGATTTGCAG
    ACAATGGATTGGGTGATGCCCCATTCCTTGATCGGCTACGCCGAGATCAAAAGTCCTTAAAAGGAAGAGG
    CAACACCCTTGGCCTCGACATCAAAACAGCCACTCTTGTTGGGAAACAAATTGTGGAATGGATTTTGAAA
    GAGGAATCCAGCGAGACACTTAGAATGGCAATTGCATCTGTACCTACTTCGCGTTACATTTCTGACATGG
    CCCTCGAGGAAATGTCACGAGACTGGTTCATGCTTATGCCTAGGCAAAAGATAATAGGCCCTCTTTGCGT
    GCGATTGGACCAGGCGGTCATGGATAAGAACGTAGTACTGGAAGCAAACTTCAGTGTAATCTTCAACCGA
    TTAGAGACCTTGATACTACTAAGGGCTTTCACTGAGGAGGGAACAATAGTTGGAGAAATTTCACCATTAC
    CTTCTCTTCCAGGACATACTTATGAGGATGTCAAAAATGCAGTTGGGGTYCTCATCGGAGGACTTGAGTG
    GAATGGTAACACGGTTCGAGTCTCTGAAAATATACAGAGATTCGCTTGGAGAAGCTGTGATGAGAATGGG
    AGACCTTCACTACCTCCAGAGCAGAAATGAGAAGTGGCGGGAACAATTGGGACAGAAATTTGAGGAAATA
    AGGTGGTTAATTGAAGAAATACGACACAGATTGAAAGCGACAGAGAATAGTTTCGAACAAATAACATTTA
    TGCAAGCCTTACAACTACTGCTTGAAGTAGAGCAAGAGATAAGAGCTTTCTCGTTTCAGCTTATTTAA

- The first line begins with '>' indicating that it is the header line.
- This is immediately followed by 'LC719646.1', which is an accession number for [this sequence in the GenBank database](https://www.ncbi.nlm.nih.gov/nuccore/LC719646.1).
- Then follows the actual nucleotide sequence, split over several lines, beginning with 'ATGGAATCCAACA...' and ending with '...TTATTTAA'.

It is very common to combine multiple sequences into a single multi-FASTA file like this:

    >ON084923.1 Influenza A virus (A/ostriches (Struthio camelus)/Egypt/Mansoura1/2022(H5N8)) segment 4 hemagglutinin, HA2 region, (HA) gene, partial cds
    GTACCACCATAGCAATGAGCAGGGGAGTGGGTACGCTGCAGACAAAGAATCCACTCAAAAGGCAATAGAT
    GGAGTTACCAATAAGGTCAACTCAATCATTGACAAAATGAACACTCAATTTGAGGCAGTTGGAAGGGAGT
    TTAATAACTTAGAAAGGAGGATAGAGAATTTGA
    
    >MW170960.1 Influenza A virus (A/swine/Italy/410927/2018(H1N2)) segment 6 neuraminidase (NA) gene, partial cds
    CCTTATGCAGATTGCTATCCTGGTAACTACTGTTACATTTCACTTCAAGCAATATGAATACAATTTCTAC
    CCAAACAACCAAGTAATGCCATGTGAACCAACGATAATTGAAAGAAACATAACAGAAATAGTGTACCTGG
    CCAACACCAC
    
    >MW170083.1 Influenza A virus (A/swine/Italy/134212/2019(H1N2)) segment 6 neuraminidase (NA) gene, partial cds
    GTAGTAACTGCCTGAGTCCTAATAATGAAGAAGGGGGTCATGGGGTAAAAGGCTGGGCCTTTGATGATGG
    AAATGATGTTTGGATGGGAAGAACGATCAGCGAAAAGTTACGATTAGGTTATGAAACCTTCAAGGTCATC
    GACGGTTGGTCCAAGCC
    
    >MW169741.1 Influenza A virus (A/swine/Italy/8745/2019(H3N2)) segment 2 polymerase PB1 (PB1) gene, partial cds
    TCGTTCCATCCTCAATACTAGCCAAAGGGGAATTCTTGAGGATGAGCAAATGTATCAGAAGTGCTGCAAT
    TTATTTGAGAAATTCTTCCCTAGCAGTTCATACAGGAGGCCAGTGGGAATTTCAAGCATGGTGGAGGCCA
    TGGTATCTAGGGCCAGAATTGATGCACGGATTGATTTCGAGTCTGGAAGGATTAATAAAGAAGAATTTGC
    TGAGATCATGAAGATCTGTTCCACCATAGAAGAGTTCAGACGGCAAAAGTAG
    
    >OM149369.1 Influenza A virus (A/Hilly chicken/Bangladesh/Avian Influenza Virus/2019(H9)) segment 4 hemagglutinin (HA) gene, partial cds
    AATTTCTTAGCTAGCAAAATGGAAACAATAACACTGATGACTACACTACTATTAACAACAACGAGCCTTG
    CAGACAAAATCTGTATCGGCCACCAATCGACAAATTCTACAGAAACTGTAGACACACTAACAGAAACTAA
    CGTTCCTGTGACACATGCCAAAGAGTTGCTCCATACGGATCACAATGGAATGCTGTGTGCAACAAATCTA
    GGACATCCCCTCATCCTAGATAAATGTAACGTAGAAGGACTGATCTACGGCAACCCTTCTTGTGATCT


If you want a more detailed history of the FASTA file format, then you could take a look at the Wikipedia page here: https://en.wikipedia.org/wiki/FASTA_format.

## FASTQ
The widely used FASTA file format has the great advantage of simplicity. However, this simplicity can be restrictive if we want to include additional data/metadata in addition to the sequence.
Given the non-negligible error rates of NGS technologies, often we need to accompany our sequence data with quality scores that estimate our confidence in the accuracy of the sequence data. As we will see later, this allows us to perform quality control checks and filter-out poor-quality data before performing analyses.
FASTQ is a simple text-based format that allows us to include quality scores. A single sequence is represented by four lines of text:

    @ERR8261968.1 1 length=97
    ACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTA
    +ERR8261968.1 1 length=97
    CCCCCFDDFFFFGGGGGGGGGGHHHHHHHHHHHGGGGHHHHHHHHHHHHHHHGHHGHHIIHHGGGGGGHHHHHHHHHHHHHHHHHHHGGGHHHHHHH

- The first line is a 'header' containing a unique identifier for the sequence and, optionally, further description.
- The second line contains the actual nucleotide sequence.
- The third line is redundant  and can be safely ignored. Sometimes it simply repeats the first line. Sometimes it is blank or just contains a '+' character.
- The fourth line contains a string of characters that encode quality scores for each nucleotide in the sequence. Each single character encodes a score, typically   a number between 0 and 40; this score is encoded by a single character, as we saw during the introductory lecture.

| Character | ASCII | FASTQ quality score (ASCII – 33) 
| --|--|--
| ! | 33 | 0
| “ | 34 | 1
| # | 35 | 2
| $ | 36 | 3
| % | 37 | 4
| ... | ... | ...
| C | 67 | 34
| D | 68 | 35
| E | 69 | 36
| F | 70 | 37
| G | 71 | 38
| H | 72 | 39
|40 | 73 | 40

So, in the example above, we can see that most of the positions within the 97-nucleotide sequence have scores in the high 30s, which indicates a high degree of confidence in their accuracy.
- A score of 30 denotes a 1 in 1000 chance of an error, i.e. 99.9 %accuracy.
- A score of 40 denotes a 1 in 10,000 chance of an error, i.e. 99.99 %accuracy.

You can read more about the FastQ file format and quality scores here:
Cock, P. J., Fields, C. J., Goto, N., Heuer, M. L., & Rice, P. M. (2010). The Sanger FASTQ file format for sequences with quality scores, and the Solexa/Illumina FASTQ variants. *Nucleic Acids Research*, **38**, 1767–1771. https://doi.org/10.1093/nar/gkp1137.


## SAM and BAM

A SAM file (usually named *.sam) is used to represent aligned sequences. It is particularly useful for storing the results of aligning genomic or transcriptomic sequence reads aligned against a reference genome sequence. The BAM file format is a compressed form of SAM. This has the disadvantage that it is not readable by a human but has the advantage of being smaller than the corresponding SAM file and thus easier to share and copy between locations.

You can read about SAM and BAM formats here:
 - Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., & 1000 Genome Project Data Processing Subgroup (2009). The Sequence Alignment/Map format and SAMtools. *Bioinformatics*, **25**, 2078–2079. https://doi.org/10.1093/bioinformatics/btp352 and
-  [https://samtools.github.io/hts-specs/SAMv1.pdf](https://samtools.github.io/hts-specs/SAMv1.pdf).

We can view BAM files graphically using a specialised genome browser software such as:
- [IGV](https://igv.org/)
- [Tablet](https://ics.hutton.ac.uk/tablet/)
- [Artemis / BAMview](http://sanger-pathogens.github.io/Artemis/BamView/) 

# Public repositories of NGS data 
The Sequence Read Archive (SRA) contains a huge number of sequence reads generated by various NGS methods. We can browse this data on the web via the NCBI's web portal. We can also download NGS datasets in FastQ format and analyse them locally, for example in our virtual machine. Let's a take a look at an example dataset: [SRR19504912](https://www.ncbi.nlm.nih.gov/sra/?term=SRR19504912)

. Which virus does this sequencing dataset come from?

Now let's download the sequence data from this sequencing run from the SRA. Unfortunately it is not easy to download the data directly from the NCBI website; instead we have to use the *fasterq-dump* tool from the [NCBI's SRA Toolkit](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit). I recommend that we work in the directory /home/manager/ViralBioinfAsia2022/course_data/NGS_file_formats_and_data_QC/. So, first execute this command in the Terminal:

    cd /home/manager/ViralBioinfAsia2022/course_data/NGS_file_formats_and_data_QC

The SRA Toolkit software is already installed on the virtual machine and we can execute it in the Terminal simply like this:

    fasterq-dump SRR19504912 

You should then see some output something like this:

    spots read      : 306,691
    reads read      : 613,382
    reads written   : 613,382


And when you type `ls -lh` you will notice that new files have been created called *SRR19504912_1.fastq*  and *SRR19504912_2.fastq*. There are two files because this dataset consists of paired sequence reads. 

Take a look at the first few lines of each file by executing this command:

    head *.fastq

How many sequence reads are in each file? Use `wc -l *.fastq` to count the numbers of lines in each file. How many lines does a single sequence read take up in a FASTQ file?

# Quality control for FastQ-formatted data

## Visualising quality metrics using FastQC

Let's examine the quality of these sequence data using a really useful piece of software called [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

We can launch the graphical interface by simply executing `fastqc` on the Terminal command line. However, it is often more convenient to use the software in  command-line mode. Execute the following command in the Terminal:

    fastqc *.fastq

This will generate some messages something like this:

    Started analysis of SRR19504912_1.fastq
    Approx 5% complete for SRR19504912_1.fastq
    Approx 10% complete for SRR19504912_1.fastq
    ...
    Approx 90% complete for SRR19504912_1.fastq
    Approx 95% complete for SRR19504912_1.fastq
    Analysis complete for SRR19504912_1.fastq
    Started analysis of SRR19504912_2.fastq
    Approx 5% complete for SRR19504912_2.fastq
    Approx 10% complete for SRR19504912_2.fastq
    Approx 15% complete for SRR19504912_2.fastq
    ...
    Approx 90% complete for SRR19504912_2.fastq
    Approx 95% complete for SRR19504912_2.fastq
    Analysis complete for SRR19504912_2.fastq


Now, execute the command `ls- lh` and you should see some new files have appeared:

 

    total 400M
    -rw-rw-r-- 1 manager manager 1.2K Jul 25 16:20 readME.md
    -rw-rw-r-- 1 manager manager 199M Jul 31 13:26 SRR19504912_1.fastq
    -rw-rw-r-- 1 manager manager 666K Jul 31 13:31 SRR19504912_1_fastqc.html
    -rw-rw-r-- 1 manager manager 378K Jul 31 13:31 SRR19504912_1_fastqc.zip
    -rw-rw-r-- 1 manager manager 199M Jul 31 13:26 SRR19504912_2.fastq
    -rw-rw-r-- 1 manager manager 673K Jul 31 13:31 SRR19504912_2_fastqc.html
    -rw-rw-r-- 1 manager manager 391K Jul 31 13:31 SRR19504912_2_fastqc.zip


We are most interested in the HTML files, which contain the FastQC reports for our two fastq files. Let's open them with the following command:

    firefox *.html &

You should then see something like this:

![enter image description here](https://github.com/WCSCourses/ViralBioinfAsia2022/raw/main/course_data/NGS_file_formats_and_data_QC/images/Screenshot%202022-07-31%20at%2015.28.02.png)










## Trimming and filtering to remove poor-quality data







