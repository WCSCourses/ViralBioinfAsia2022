# [Viral Genomics and Bioinformatics Asia 2022 course](https://coursesandconferences.wellcomeconnectingscience.org/event/viral-genomics-and-bioinformatics-asia-20220822/)
* Monday 22nd - Friday 26th August 2022  
* [Wellcome Connecting Science](https://www.wellcomeconnectingscience.org) and [COG-Train](https://www.cogconsortium.uk/priority-areas/training/about-cog-train/)  
* [Course website](https://coursesandconferences.wellcomeconnectingscience.org/event/viral-genomics-and-bioinformatics-asia-20220822/)  
* [Course GitHub repository](https://github.com/WCSCourses/ViralBioinfAsia2022)

[Richard Orton](https://www.gla.ac.uk/schools/infectionimmunity/staff/richardorton/)   
[MRC-University of Glasgow Centre for Virus Research](https://www.gla.ac.uk/research/az/cvr/)  
464 Bearsden Road  
Glasgow  
G61 1QH  
UK  
E-mail: Richard.Orton@glasgow.ac.uk  

## Contents

**CONTENT IS STILL UNDER DEVELOPMENT - DO NOT USE BEFORE FRIDAY 26TH AUGUST 2022!**

* [1: SARS-CoV-2 Reference Alignment](#1-sars-cov-2-reference-alignment)
	+ [1.1: MinION Tutorial](#11-minion-tutorial)
		+ [1.1.1: Setup and Data](#111-setup-and-data)
		+ [1.1.2: ARTIC consensus sequence generation walkthrough](#112-artic-consensus-sequence-generation-walkthrough)
		+ [1.1.3: Exercise: Generating ARTIC conensus sequences yourself](#113-exercise-generating-artic-conensus-sequences-yourself)
		+ [1.1.4: Additional Resources](#114-additional-resources)
	+ [1.2: Illumina Tutorial](#12-illumina-tutorial)
		+ [1.2.1: Setup and data](#121-setup-and-data)
		+ [1.2.2: Illumina consensus generation walkrhough](#122-illumina-consensus-generation-walkthrough)
		+ [1.2.3: Exercise: Generating Illumina conensus sequences yourself](#123-exercise-generating-illumina-conensus-sequences-yourself)

* [2: SARS-CoV-2 Lineages and Mutations](#2-sars-cov-2-lineages-and-mutations)
* [3: SARS-CoV-2 Phylogenetics](#3-sars-cov-2-phylogenetics)
* [4: SARS-CoV-2 Group Practical](#4-sars-cov-2-group-practical)

## 1: SARS-CoV-2 Reference Alignment

In this module we will be performing reference alignment of SARS-CoV-2 (Severe acute respiratory syndrome coronavirus 2) samples to create consensus sequences for a number of samples. This will be done for both MinION and Illumina data sets using a different approach for each technology.

By this point in the course (day 5) you should be comfortable working with the terminal command line on the course Ubuntu linux virtual machine: navigating through folders using ```cd```, list folder contents using ```ls```, making directories using ```mkdir```, deleting files using ```rm```, entering bioinformatics commands, and using **TAB Completetion** which makes entering long filenames and paths easy.

Previosuly in the course you would of learnt about the FASTQ format, what SAM and BAM files, and been aligning primarily illumina reads to reference sequences to call a consensus sequence. This session will build upon this, tweaking the steps to adapt them for handling the large number of overlapping amplicons used in the ARTIC SARS-CoV-2 protocols, and using a consensus caller specifically designed for viral samples.

Commands that you need to enter into the terminal are presented in this tutorial in grey code boxes like this:

```
cd  dont_enter_this_command_its_just_an_example
```
**NB:** all commands presented here are single line commands, but if the command is long with many arguments or lengthy file paths it may get wrapped around onto multiple lines, wrapping does not signify new lines, only hit the enter button at the end of a command.

## 1.1: MinION Tutorial
This MinION tutorial is largely based on the [ARTIC](https://artic.network) network's [nCoV-2019 bioinformatics protocol](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html) which we will use to create consensus genome sequences for a number of MinION SARS-CoV-2 samples. 

### 1.1.1: Setup and Data

The ARTIC nCoV-2019 [Conda](https://docs.conda.io/en/latest/) enviroment has already been installed on the [VirtualBox](https://www.virtualbox.org) Ubuntu virtual machine (VM) that you are using for this course, but we need to 'activate' the Conda environment (which is called artic-ncov2019) each time we use it, in order to use the ARTIC pipeline:

```
conda activate artic-ncov2019
```

Next, lets change directory (cd) into the folder where the example MinION data is located for this practical:

```
cd ~/SARS-CoV-2/MinION/20201229_1542_X1_FAO14190_c9e59aa7_Batch124A/

```

**NB:** Using **TAB completion** makes navigating into long folder names easy!!

This run was performed on an [Oxford Nanopore Technologies](https://nanoporetech.com) (ONT) GridION machine, the run was started on 29th Decemeber 2020 at 15:42 (20201229_1542), using the first (X1) of the five flow cells slots on the GridION, the flowcell ID number was FAO14190, and the run was named 'Batch124A' locally within the [Medical Research Council-University of Glasgow Centre for Virus Research](https://www.gla.ac.uk/research/az/cvr/) (CVR) as part of a Covid-19 Genomics UK Consotium ([COG-UK](https://www.cogconsortium.uk)) sequencing run. The samples were sequenced used Version 2 (V2) of the ARTIC [nCoV-2019](https://github.com/artic-network/primer-schemes/tree/master/nCoV-2019) amplcion primers.

The run was live basecalled and demultiplexed using the ONT basecaller Guppy (this is available for download from the [ONT website](https://nanoporetech.com) after registration with ONT), making use of the GridION's onboard Graphical Processing Unit (GPU). If you list the contents of this directory you should see a number of files and folders:

```
ls
```

The key files and folders in this run folder are (amongst others):

* **fast5\_pass** - this folder contains all the raw FAST5 reads that PASSED the basic quality control filters of the guppy basecaller. 
* **fast5\_fail** - this folder contains all the raw FAST5 reads that FAILED the basic quality control filters of the guppy basecaller
* **fastq\_pass** - this folder contains all the FASTQ reads that were converted from the those within tge fast5\_pass fiolder
* **fastq_fail**- this folder contains all the FASTQ reads that were converted from the those within tge fast5\_fail fiolder
* **sequencing\_summary\_FAO14190\_ad60b376.txt** - this sequencing_summary file is produced by the basecaller and contains a summary of each read such as it's name, length, barcode and what FAST5 and FASTQ files it is loctaed it.

**NB:** To save disk space on the VM, the fast5\_fail and fastq\_fail folders do not contain any read data, as failed reads are not used in the pipeline.

As the data has already been demultiplexed, there is one folder for each barcode detected on the run:

```
ls fastq_pass
```

You should see the following barcode folders (06, 07 and 12) representing the 3 samples on the run that we will be analysing:

* **barcode06**
* **barcode07**
* **barcode12**
* **unclassified**

**NB**: unclassified is where reads whose barcode could not be determined are placed - to save disk space on the VM there are no reads in the unclassified folders, but they typically contain substantial amounts of reads.

Typically the FASTQ data for each sample is stored in multiple files of around 4000 reads each. For barcode06, you should see 25 different FASTQ files, numerically labelled at the end of their filename from 0 to 24:

```
ls fastq_pass/barcode06
```


### 1.1.2: ARTIC consensus sequence generation walkthrough

The first sample we will be working with is barcode06. The ARTIC bioinformatics protocol has two distinct steps:

1. **artic guppyplex** - combines all a samples FASTQ reads into a single file and size filters them (it can also perform a quality score check, which is not needed here as the reads are already split into pass and fail folders based on quality by guppy)
2. **artic minion** - aligns the reads to the [Wuhan-Hu-1](https://www.ncbi.nlm.nih.gov/nuccore/MN908947) reference sequence, trims the amplicon primer sequences from the aligned reads, downsamples amplicons to reduce the data, and creates a consensus seuqnece utilising [nanopolish](https://github.com/jts/nanopolish) for variant calling to correct for common MinION errors (such as those associated with homopolymer regions).

First we will create a folder to work in and store our output files:

```
mkdir ~/SARS-CoV-2/MinION_Results
```

Then we will move into the folder to work:

```
cd ~/SARS-CoV-2/MinION_Results
```

**artic guppyplex** - now we will run artic guppyplex on sample barcode06:

```
artic guppyplex --skip-quality-check --min-length 400 --max-length 700 --directory ~/SARS-CoV-2/MinION/20201229_1542_X1_FAO14190_c9e59aa7_Batch124A/fastq_pass/barcode06 --prefix cvr124a
```

Breaking this command down:

* **artic gupplyplex** = the name of the program/module to use (installed as part of conda environment)
* **--skip-quality-check** = don't filter reads based on quality (our reads are already filtered)
* **--min-length 400** = minimum read length to accept is 400 bases
* **--max-length 700** = maximum read length to accept is 700 bases
* **--directory** = PATH to input directory containing FASTQ reads to process
* **--prefix** = output name prefix to label output file (I choose cvr124a to signify batch 124a from the CVR)

This should create an output FASTQ file called **cvr124a_barcode06.fastq**:

```
ls
```

**artic minion** - next we will run artic minion using the FASTQ file created above:


```
artic minion --normalise 200 --threads 4 --scheme-directory ~/artic-ncov2019/primer_schemes --read-file cvr124a_barcode06.fastq --fast5-directory ~/SARS-CoV-2/MinION/20201229_1542_X1_FAO14190_c9e59aa7_Batch124A/fast5_pass/barcode06 --sequencing-summary ~/SARS-CoV-2/MinION/20201229_1542_X1_FAO14190_c9e59aa7_Batch124A/sequencing_summary_FAO14190_ad60b376.txt nCoV-2019/V2 barcode06
```

Breaking this comamnd down:

* **artic minion** = the name of the program/module to use (installed as part of conda environment)
* **--normalise 200** = normalise (downsample) each amplicon so there are only 200 reads in each direction (forward and reverse) - this enables the nanopolish variant and consensus calling to complete realtively quickly
* **--threads 4** = the number of computer threads to use (depends on how poerful your machine is, the more the better)
* **--scheme-directory** = path to the artic primer scheme directory (installed on the VM as part of the conda environment)
* **--read-file** = the name of the input FASTQ file to align
* **--fast5-directory** = the path to the corresponding FAST5 folder of the reads
* **--sequencing-summary** = the path to the sequencing_summary.txt file of the run 
* **nCoV-2019/V2** = the primer scheme to use for amplicon primer trimming - this folder is located in the scheme\_directory so on this VM this corresponds to the folder ~/artic-ncov2019/primer_schemes/nCoV-2019/V2
* **barcode06** = the output prefix name to label output files, this can be anything you want such as the sample name or barcode number or anything you want

Overall, this artic minion command uses the aligner [minimap2](https://github.com/lh3/minimap2) to align the reads to the [Wuhan-Hu-1](https://www.ncbi.nlm.nih.gov/nuccore/MN908947) reference sequence, [samtools](http://www.htslib.org) to convert and sort the SAM file into BAM format, custom [artic](https://github.com/artic-network/fieldbioinformatics) scripts for amplicon primer trimming and normalisation (downsampling), [nanopolish](https://github.com/jts/nanopolish) for variant calling, and custom [artic](https://github.com/artic-network/fieldbioinformatics) scripts for creating the consensus sequence using the reference and VCF files, and then masking low coverage regions with Ns. This will create the following key files (amongst many others), all starting with the prefix **barcode06** in this instance:

* **barcode06.sorted.bam** - BAM file containing all the reads aligned to the reference sequence (there is no amplicon primer trimming in this file)
* **barcode06.trimmed.rg.sorted.bam** - BAM file containing normalised (downsampled) reads with amplicon primers left on - this is the file used for variant calling
* **barcode06.primertrimmed.rg.sorted.bam** - BAM file containing normalised (downsampled) reads with amplicon primers trimmed off
* **barcode06.pass.vcf.gz** - detected variants that PASSed the filters in VCF format (gzipped)
* **barcode06.fail.vcf** - detected variants that FAILed the filters in VCF format
* **barcode06.consensus.fasta** - the consensus sequence of the sample

**NB:** this command will give a warning message at the end, this can be ignored as we don't need the muscle alignment file. The warning is due to an updated version of muscle within the artic environment that expects the command structured -align and -ouput rather thab -in and -out (also reported as a GitHub issue [here](https://github.com/artic-network/artic-ncov2019/issues/93)):

```
DO NOT ENTER THIS - IT IS A RECORD OF THE ERROR MESSAGE
Command failed:muscle -in barcode06.muscle.in.fasta -out barcode06.muscle.out.fasta
```

We can count the number of reads that have mapped to the reference sequence (discounting supplementary and additional alignments which are multi-mappings) using samtools view with the -c argument to count reads, and by utilising [SAM Flags](https://samtools.github.io/hts-specs/SAMv1.pdf) to only count read alignments that are mapped (F4), that are primary alignments (F256), and are not supplementary alignments (F2048): F4 + F256 + F2048 = F2308:

```
samtools view -c -F2308 barcode06.sorted.bam
```

If you compare the mapped read count to that from the normalised (downsampled) BAM file, you should see less mapped reads due to the normalisation step:

```
samtools view -c -F2308 barcode06.trimmed.rg.sorted.bam
```

We can view the FASTA consensus sequence via the command line (we will be analysing consensus sequences in more depth in later sessions to see what mutations they contain etc):

```
more barcode06.consensus.fasta
```

### 1.1.3: Exercise: Generating ARTIC conensus sequences yourself

Your task now is to adapt the above artic gupplyplex and minion commands to run on samples barcode07 and/or barcode12.

A reminder of the two commands used for barcode06 is:


```
artic guppyplex --skip-quality-check --min-length 400 --max-length 700 --directory ~/SARS-CoV-2/MinION/20201229_1542_X1_FAO14190_c9e59aa7_Batch124A/fastq_pass/barcode06 --prefix cvr124a
```

```
artic minion --normalise 200 --threads 4 --scheme-directory ~/artic-ncov2019/primer_schemes --read-file cvr124a_barcode06.fastq --fast5-directory ~/SARS-CoV-2/MinION/20201229_1542_X1_FAO14190_c9e59aa7_Batch124A/fast5_pass/barcode06 --sequencing-summary ~/SARS-CoV-2/MinION/20201229_1542_X1_FAO14190_c9e59aa7_Batch124A/sequencing_summary_FAO14190_ad60b376.txt nCoV-2019/V2 barcode06
```

Essentially all you will need to do is change every occurence of **barcode06** (once in guppyplex and three times in minion) to either **barcode07** or **barcode12**. 

**QUESTION** - what is number of mapped reads in each of the samples you have looked at?

As the commands are well structured and all that is needed to change is the input and output names within a run, it means that these commands are easily scriptable using bash (which you have learnt earilier on in this course) or pipeline tools such as [snakemake](https://snakemake.github.io) and [nextflow](https://www.nextflow.io).

At the end of this session we should decativate our conda environment:

```
conda deactivate
```

### 1.1.4: Additional Resources


We were initally planning on using the [ncov2019-artic-nf](https://github.com/connor-lab/ncov2019-artic-nf) nextflow pipeline for running the ARTIC network tools as it offers a number of extra functions and easier automation. However, we could not get it working on the VM correctly (DSL error resolved by replacing nextflow.preview.dsl = 2 with nextflow.enable.dsl = 2 within the scripts, but the following error could not be resolved by dropping down Nextflow versions (GitHub issue to be created):

```
Unqualified output value declaration has been deprecated - replace `tuple sampleName,..` with `tuple val(sampleName),..`
```

We did also try using the [nf-core/viral-recon](https://nf-co.re/viralrecon) viral netflow pipeline, which is highly recommended (not just for SARS-CoV-2), but performance was not great on our small VM.

There is a very nice [SARS-CoV-2 sequencing data analysis tutorial](https://github.com/cambiotraining/sars-cov-2-genomics) created by CamBioTraining, which includes [nf-core/viral-recon](https://nf-co.re/viralrecon) instructions.

This MinION tutorial is based on the original [ARTIC](https://artic.network) network's [nCoV-2019 bioinformatics protocol](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html).

## 1.2: Illumina Tutorial

In the above tutorial, we were working with MinION data. In this tutroial we will be using paired end Illumina SARS-CoV-2 data. Although this uses different computational tools, the two approaches are in essence very similar.

### 1.2.1: Setup and data

We do not need to use a conda enviroment as all the tools we need are already installed direclty on the VM: 

* [trim_galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) for read trimming (the reads are pre-trimmed using trim_galore)
* [bwa](https://github.com/lh3/bwa) for read alignment
* [samtools](http://www.htslib.org) for SAM/BAM conversion
* [ivar](https://andersen-lab.github.io/ivar/html/manualpage.html) for primer trimming and consensus calling
* [weeSAM](https://github.com/centre-for-virus-research/weeSAM) for coverage plots.

First, lets move into the Illumina data directory:

```
cd ~/SARS-CoV-2/Illumina/200703_M01569_0148_000000000-J53HN_Batch70/
```

The data in this folder is from a run on an Illumina MiSeq machine. The name of the folder implies it was run on the 3rd July 2020 (200703), the machine ID is M01569, the run ID is 0148_000000000-J53HN, and this was called Batch70 locally within the [Medical Research Council-University of Glasgow Centre for Virus Research](https://www.gla.ac.uk/research/az/cvr/) (CVR) as part of a Covid-19 Genomics UK Consotium ([COG-UK](https://www.cogconsortium.uk)) sequencing run. The samples were sequenced used Version 1 (V1) of the ARTIC [nCoV-2019](https://github.com/artic-network/primer-schemes/tree/master/nCoV-2019) amplcion primers.

There are four samples in this run called:

* CVR2058
* CVR2078
* CVR2092
* CVR2101

If you list the contents of the directory you should see paired end reads (R1.fastq and R2.fastq) for each of the four samples:

```
ls
```

As these are paired end reads, each sample's R1 and R2 read files should contain the same number of reads (and therefore the same number of lines). To check, lets count the number of lines in every FASTQ file in the directory (samples will have different read numbers, but the R1 and R2 files from a single sample should have the same number):

```
wc -l *.fastq
```

As each FASTQ read consists of four lines, the number of lines should be divided by 4 to get the number of reads, this can be done on a file by file basis on the command line using:

```
expr `(wc -l CVR2058_R1.fastq | cut -f1 -d " ")` / 4
```
**NB:** this first counts the number of lines (```wc -l```) in the file, ```wc``` outputs the number of lines followed by a space and the name of the file (e.g. 561120 CVR2058_R1.fastq), ```cut``` splits this output using the space delimiter (```-d " "```) and we select only the 1st column (```-f1```), this number (e.g. 561120) is then divided by 4 (```/ 4```) using ```expr```.

These samples have already been trimmed/filtered using [trim_galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/), so we do not need to do QC them. For information, this is command that was used for each sample:

```
DO NOT ENTER THE BELOW COMMAND- IT IS JUST FOR INFORMATION 
trim_galore -q 20 --length 50 --paired SAMPLE_R1.fastq SAMPLE_R2.fastq
```
### 1.2.2: Illumina consensus generation walkthrough

We will be using the [Wuhan-Hu-1](https://www.ncbi.nlm.nih.gov/nuccore/MN908947) SARS-CoV-2 isolate as the reference sequence, which has the GenBank accession number [MN908947](https://www.ncbi.nlm.nih.gov/nuccore/MN908947); this is also the sequence used as the basis for the SARS-CoV-2 RefSeq sequence [NC_045512](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512). We will use [bwa](https://github.com/lh3/bwa) to align our reads to the reference.

First we need to index the reference sequence that we will be aligning our reads to. Indexing enables the aligner to quickly look up the best places to start aligning a read by using small sequences call 'seeds' from within the read and looking up if and where those seeds occur in the reference genome using the index.

```
bwa index ~/SARS-CoV-2/MN908947.fasta 
```

This should of created a range of bwa index files (MN908947.fasta.amb/.ann/.bwt/.pac/.sa files), so list (```ls```) the contents of the directory to check:

```
ls ~/SARS-CoV-2/
```

You only need to index a genone once, if you are aligning many samples to the same genome sequence (as we will be on this course) you do not need to re-run the reference index step; but don't confuse this step with the BAM indexing step which does need to be done for each sample.

Now we will align the reads from sample CVR2058 to the reference sequence using bwa:


```
bwa mem -t4 ~/SARS-CoV-2/MN908947.fasta CVR2058_R1.fastq CVR2058_R2.fastq > CVR2058.sam

```
Breaking this command down:

* **bwa** = the name of the program
* **mem** = the name of the bwa algorithm to use (it is recommended for reads > 70bp)
* **-t4** = use 4 computational threads
* **~/SARS-CoV-2/MN908947.fasta** = the path to the reference file (and index)
* **CVR2058_R1.fastq** = FASTQ read pair 1
* **CVR2058_R2.fastq** = FASTQ read pair 2
* **> CVR2058.sam** = direct the output into the file CVR2058.sam rather than to the command line

We now have a SAM file, which we immediatley convert into a BAM file as they are compressed and faster to work with downstream. We can sort the file at the same time as converting it:


```
samtools sort -@4 CVR2058.sam -o CVR2058.bam
```
Breaking this command down:

* **samtools** = the name of the program
* **sort** = the name of the funciton within samtools to use
* **-@4** = use 4 threads
* **CVR2058.sam** = the input file
* **-o CVR2058.bam** =  the output file

We no longer need the SAM file so can delete (```rm```) it:

```
rm CVR2058.sam 

```
Now we need to index the BAM file, this makes downstream analyses faster and many downstream tools will only work if the index file has been created:

```
bwa index CVR2058.bam
```

If we list the contents of the directory we should see the index file with a .bai extension has been created:

```
ls
```


As we have used ARTIC amplicons, each read will typically start and end with a primer sequence. The primer sequence does not come from the sample's viral genome (and the primer may actually be slightly different to the virus) so all primer sequences should be removed from the read ends so they don't interfere with consensus and variant calling. 

To do this we use a tool called [ivar](https://andersen-lab.github.io/ivar/html/manualpage.html) which requires a [BED](https://software.broadinstitute.org/software/igv/BED) file containing the primer coordinates on the reference genome:

```
ivar trim -i CVR2058.bam -b ~/artic-ncov2019/primer_schemes/nCoV-2019/V1/nCoV-2019.bed -p CVR2058_trim.bam
```
**NB:** Use **TAB Completion** to help enter the primer bed file!

Breaking this command down:

* **ivar** = the name of the program
* **trim** = the name of the function within ivar to use (trimming primers)
* **-i CVR2058.bam** = the name of the input BAM fie (it is in the current directory)
* **-b ~/artic-ncov2019/primer_schemes/nCoV-2019/V1/nCoV-2019.bed** = the path to primer BED file
* **-p CVR2058_trim.bam** = the name of the output file to create (which will be a primer trimmed version of the input BAM file)

For those who are interested ivar works by soft clipping the primer sequences from the read alignments in the BAM file (rather than actually trimming the reads sequences) based on the alignment coordinates. Soft clipping is a form of trimming that is embedded within the CIGAR string of a read alignment. The CIGAR string is the 6th field of the SAM/BAM file, if you were to examine the BAM file manually you should see lots of 'S' characters in the CIGAR field: see the [CIGAR specificiation](https://samtools.github.io/hts-specs/SAMv1.pdf) for more information.

We now need to sort and then index this BAM file. Sort:

```
samtools sort -@4 CVR2058_trim.bam -o CVR2058_trim_sort.bam 
```

Rename the file back to CVR2058\_trim.bam as we don't need the unsorted BAM file:

```
mv CVR2058_trim_sort.bam CVR2058_trim.bam
```

Index:

```
bwa index CVR2058_trim.bam
```

We now have a BAM file (CVR2058\_trim.bam) that contains our samples's paired end reads (CVR2058\_R1.fastq CVR2058\_R2.fastq) aligned to the Wuhan-Hu-1 (MN908947) reference genome, with the amplicon primer sequences clipped off. So we are now ready to call a consensus sequence:

```
samtools mpileup -aa -A -d 0 -Q 0 CVR2058_trim.bam | ivar consensus -p CVR0258 -t 0.4
```

Breaking this command down, there are two parts:

1. samtools [mpileup](http://www.htslib.org/doc/samtools-mpileup.html) which essentially outputs the base and indel counts for each genome position
	* **-aa** = output data for absolutely all positions (even zero coverage ones)
	* **-A** = count orphan reads (reads whose pair did not map)
	* **-d 0** = override the maximum depth (default is 8000 which is typically too low for viruses)
	* **-Q 0** = minuimum base quality, 0 essentially means all the data
2. ivar [consensus](https://andersen-lab.github.io/ivar/html/manualpage.html) - this calls the consensus - the output of the samtools mpileup command is piped '|' directly into ivar
	* -p CVR2058 = prefix with which to name the output file
	* -t 0.4 = the minimum frequency threshold that a base must match to be used in calling the consensus base at a position. In this case, an ambiguity code will be used if more than one base is > 40% (0.4). See [ivar manual]

By default, ivar consensus uses a minimum depth (-m) of 10 and a minimum base quality (-q) of 20 to call the consensus; these defaults can be chnaged by using the appropriate arguments. If a genome position has a depth less than the minimum, an 'N' base will be used in the consensus sequence by default.

ivar will output some basic statistics to the screen such as:

```
DO NOT ENTER THIS - IT IS THE IVAR OUTPUT YOU SHOULD SEE
Reference length: 29903
Positions with 0 depth: 121
Positions with depth below 10: 121
```
and you should see our consensus sequence (CVR0258.fa) in the directory:

```
ls
```

which you can view via the command line (again, we will be doing variants in later sessions):

```
more CVR0258.fa 
```

### 1.2.3: Exercise: Generating Illumina conensus sequences yourself

There are three other samples in the Illumina data directory:

* CVR2078
* CVR2092
* CVR2101

You should now choose atleast one sample to create a consensus sequence for yourself by running through the above steps, but adapting them for the next sample (you simply need to change the input read names, and the output file names from CVR2058 to your next sample name). A reminder that the commands used were:

```
bwa mem -t4 ~/SARS-CoV-2/MN908947.fasta CVR2058_R1.fastq CVR2058_R2.fastq > CVR2058.sam
```

```
samtools sort -@4 CVR2058.sam -o CVR2058.bam
```

```
rm CVR2058.sam 
```

```
bwa index CVR2058.bam
```

```
ivar trim -i CVR2058.bam -b ~/artic-ncov2019/primer_schemes/nCoV-2019/V1/nCoV-2019.bed -p CVR2058_trim.bam
```

```
samtools sort -@4 CVR2058_trim.bam -o CVR2058_trim_sort.bam 
```

```
mv CVR2058_trim_sort.bam CVR2058_trim.bam
```

```
bwa index CVR2058_trim.bam
```

```
samtools mpileup -aa -A -d 0 -Q 0 CVR2058_trim.bam | ivar consensus -p CVR0258 -t 0.4
```

**QUESTION** - what is number of mapped reads in each of the samples you have looked at? Hint:

```
samtools view -c -F2308 input.bam
```

Overall, you should again see that we are simply running the same set of commands over and over again for different samples but just changing the input and output names. This is where the power of simple bash scripting and bioinformatics pipelines come into play, as relatively simple scripts can be written to automate this process.


## 2: SARS-CoV-2 Lineages and Mutations



* Pangolin
* SPEAR

```
cd ~/SARS-CoV-2/Variants/

```
```
grep ">" variant_seqs.fasta 
```

* Scotland/LSPA-3E62463/2022
* Scotland/LSPA-3E3B700/2022
* Scotland/QEUH-3DC59CA/2022
* Scotland/QEUH-3D90306/2022
* Scotland/QEUH-3D43C43/2022
* Scotland/CVR14531/2022
* Scotland/QEUH-36897FC/2022
* Scotland/QEUH-36491AF/2022
* Scotland/QEUH-2D86BD1/2021
* Scotland/QEUH-2D7F704/2021
* Scotland/QEUH-1BA3933/2021
* Scotland/QEUH-1B0246E/2021
* Scotland/QEUH-1725CCB/2021
* Scotland/QEUH-1585B0A/2021
* Scotland/QEUH-158D786/2021
* Scotland/CAMC-14DE972/2021
* Scotland/QEUH-147E6F5/2021
* Scotland/QEUH-13ADEF6/2021
* Scotland/QEUH-138F944/2021


## 3: SARS-CoV-2 Phylogenetics

* USHER
* CIVET

## 4: SARS-CoV-2 Group Practical

In this session, we will be working on some more Illumina paired end read data. The FASTQ data was downloaded from the [European Nucleotide Archive](https://www.ebi.ac.uk/ena/browser) (ENA), and there are 4 samples in total (the samples are not related to one another), with R1 and R2 FASTQ files for each:

* ERR9105817
* ERR9731990
* ERR9761275
* ERR9788433

Your task is to work as a group to analyse this samples. Initial read QC (with trim_galore) is not required. You should:

* align the reads to the Wuhan-Hu-1 reference sequence
* Report the number of mapped reads
* Trim the primers (primer versions will be added below)
* Call a consensus sequence
* Use pangolin to determine the linage
* Use SPEAR to call the mutations

One of the goals of this session is to try and write a bash script that does all of the above automatically for your. Remember, all the commands are the same each time, we are just varying the input file names and output file names.

Some initial bash scritps will be provided here to help get you started which you will need to exapnd to get working.

Completed example bash scripts will be uploaded here after the session.







