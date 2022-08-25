# [Viral Genomics and Bioinformatics Asia 2022 course](https://coursesandconferences.wellcomeconnectingscience.org/event/viral-genomics-and-bioinformatics-asia-20220822/)
* Monday 22nd - Friday 26th August 2022  
* [Wellcome Connecting Science](https://www.wellcomeconnectingscience.org) and [COG-Train](https://www.cogconsortium.uk/priority-areas/training/about-cog-train/)  
* [Course website](https://coursesandconferences.wellcomeconnectingscience.org/event/viral-genomics-and-bioinformatics-asia-20220822/)  
* [Course GitHub repository](https://github.com/WCSCourses/ViralBioinfAsia2022)

## Contact

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
		+ [1.1.3: Exercise: Generating ARTIC consensus sequences yourself](#113-exercise-generating-artic-consensus-sequences-yourself)
		+ [1.1.4: Additional Resources](#114-additional-resources)
	+ [1.2: Illumina Tutorial](#12-illumina-tutorial)
		+ [1.2.1: Setup and data](#121-setup-and-data)
		+ [1.2.2: Illumina consensus generation walkthrough](#122-illumina-consensus-generation-walkthrough)
		+ [1.2.3: Exercise: Generating Illumina consensus sequences yourself](#123-exercise-generating-illumina-consensus-sequences-yourself)
* [2: SARS-CoV-2 Lineages and Mutations](#2-sars-cov-2-lineages-and-mutations)
	+ [2.1: Pangolin lineages](#21-pangolin-lineages)
	+ [2.2: SPEAR Mutations](#22-spear-mutations)
* [3: SARS-CoV-2 Phylogenetics](#3-sars-cov-2-phylogenetics)
	+ [3.1: USHER](#31-usher)
	+ [3.2: CIVET](#32-civet)
* [4: SARS-CoV-2 Group Practical](#4-sars-cov-2-group-practical)
* [5: Warnings](#5-warnings)

## 1: SARS-CoV-2 Reference Alignment

In this module we will be performing reference alignment of SARS-CoV-2 (Severe acute respiratory syndrome coronavirus 2) samples to create consensus sequences for a number of samples. This will be done for both MinION and Illumina data sets using a different approach for each technology.

By this point in the course (day 5) you should be comfortable working with the terminal command line on the course Ubuntu linux virtual machine: navigating through folders using ```cd```, list folder contents using ```ls```, making directories using ```mkdir```, deleting files using ```rm```, entering bioinformatics commands, and using **TAB Completion** which makes entering long filenames and paths easy.

Previously in the course you would of learnt about the FASTQ format, what SAM and BAM files, and been aligning primarily illumina reads to reference sequences to call a consensus sequence. This session will build upon this, tweaking the steps to adapt them for handling the large number of overlapping amplicons used in the ARTIC SARS-CoV-2 protocols, and using a consensus caller specifically designed for viral samples.

Commands that you need to enter into the terminal are presented in this tutorial in grey code boxes like this:

```
cd  dont_enter_this_command_its_just_an_example
```
**NB:** commands are presented within code blocks - some are long and stretch off the page within a scrollpane.

## 1.1: MinION Tutorial
This MinION tutorial is largely based on the [ARTIC](https://artic.network) network's [nCoV-2019 bioinformatics protocol](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html) which we will use to create consensus genome sequences for a number of MinION SARS-CoV-2 samples. 

### 1.1.1: Setup and Data

The 'artic-ncov2019' [Conda](https://docs.conda.io/en/latest/) environment has already been installed on the [VirtualBox](https://www.virtualbox.org) Ubuntu virtual machine (VM) that you are using for this course, but we need to 'activate' the artic-ncov2019 Conda environment each time we want to use the ARTIC pipeline:

```
conda activate artic-ncov2019
```

Next, lets change directory (```cd```) into the folder where the example MinION data is located for this practical:

```
cd ~/SARS-CoV-2/MinION/20201229_1542_X1_FAO14190_c9e59aa7_Batch124A/

```

**NB:** Using **TAB completion** makes navigating into long folder names easy!!

This run was performed on an [Oxford Nanopore Technologies](https://nanoporetech.com) (ONT) GridION machine, the run was started on 29th December 2020 at 15:42 (20201229_1542), using the first (X1) of the five flow cells slots on the GridION, the flowcell ID number was FAO14190, and the run was named 'Batch124A' locally within the [Medical Research Council-University of Glasgow Centre for Virus Research](https://www.gla.ac.uk/research/az/cvr/) (CVR) as part of a Covid-19 Genomics UK Consortium ([COG-UK](https://www.cogconsortium.uk)) sequencing run. The samples were sequenced used Version 2 (V2) of the ARTIC [nCoV-2019](https://github.com/artic-network/primer-schemes/tree/master/nCoV-2019) amplicon primers.

The run was live basecalled and demultiplexed using the ONT basecaller Guppy (this is available for download from the [ONT website](https://nanoporetech.com) after registration with ONT), making use of the GridION's onboard Graphical Processing Unit (GPU). If you list the contents of this directory you should see a number of files and folders:

```
ls
```

The key files and folders in this run folder are (amongst others):

* **fast5\_pass** - this folder contains all the raw FAST5 reads that PASSED the basic quality control filters of the guppy basecaller. 
* **fast5\_fail** - this folder contains all the raw FAST5 reads that FAILED the basic quality control filters of the guppy basecaller
* **fastq\_pass** - this folder contains all the FASTQ reads that were converted from the those within tge fast5\_pass fiolder
* **fastq_fail**- this folder contains all the FASTQ reads that were converted from the those within tge fast5\_fail fiolder
* **sequencing\_summary\_FAO14190\_ad60b376.txt** - this sequencing_summary file is produced by the basecaller and contains a summary of each read such as it's name, length, barcode and what FAST5 and FASTQ files it is located in.

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
2. **artic minion** - aligns the reads to the [Wuhan-Hu-1](https://www.ncbi.nlm.nih.gov/nuccore/MN908947) reference sequence, trims the amplicon primer sequences from the aligned reads, downsamples amplicons to reduce the data, and creates a consensus sequence utilising [nanopolish](https://github.com/jts/nanopolish) for variant calling to correct for common MinION errors (such as those associated with homopolymer regions).

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

Breaking this command down:

* **artic minion** = the name of the program/module to use (installed as part of conda environment)
* **--normalise 200** = normalise (downsample) each amplicon so there are only 200 reads in each direction (forward and reverse) - this enables the nanopolish variant and consensus calling to complete relatively quickly
* **--threads 4** = the number of computer threads to use (depends on how powerful your machine is, the more the better)
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

**NB:** this command will give a warning message at the end, this can be ignored as we don't need the muscle alignment file. The warning is due to an updated version of muscle within the artic environment that expects the command structured -align and -output rather than -in and -out (also reported as a GitHub issue [here](https://github.com/artic-network/artic-ncov2019/issues/93)):

```
#DO NOT ENTER THIS - IT IS A RECORD OF THE ERROR MESSAGE
Command failed:muscle -in barcode06.muscle.in.fasta -out barcode06.muscle.out.fasta
```

We can count the number of reads that have mapped to the reference sequence (discounting supplementary and additional alignments which are multi-mappings) using samtools view with the -c argument to count reads, and by utilising [SAM Flags](https://samtools.github.io/hts-specs/SAMv1.pdf) to only count read alignments that are mapped (F4), that are primary alignments (F256), and are not supplementary alignments (F2048): F4 + F256 + F2048 = F2308:

```
samtools view -c -F2308 barcode06.sorted.bam
```

**NB:** Alternatively, you could use the command ```samtools flagstats barcode06.sorted.bam``` which you learnt previously.

If you compare the mapped read count to that from the normalised (downsampled) BAM file, you should see less mapped reads due to the normalisation step:

```
samtools view -c -F2308 barcode06.trimmed.rg.sorted.bam
```

We can view the FASTA consensus sequence via the command line (we will be analysing consensus sequences in more depth in later sessions to see what mutations they contain etc):

```
more barcode06.consensus.fasta
```

### 1.1.3: Exercise: Generating ARTIC consensus sequences yourself

Your task now is to adapt the above artic guppyplex and minion commands to run on samples barcode07 and/or barcode12.

A reminder of the two commands used for barcode06 is:


```
artic guppyplex --skip-quality-check --min-length 400 --max-length 700 --directory ~/SARS-CoV-2/MinION/20201229_1542_X1_FAO14190_c9e59aa7_Batch124A/fastq_pass/barcode06 --prefix cvr124a
```

```
artic minion --normalise 200 --threads 4 --scheme-directory ~/artic-ncov2019/primer_schemes --read-file cvr124a_barcode06.fastq --fast5-directory ~/SARS-CoV-2/MinION/20201229_1542_X1_FAO14190_c9e59aa7_Batch124A/fast5_pass/barcode06 --sequencing-summary ~/SARS-CoV-2/MinION/20201229_1542_X1_FAO14190_c9e59aa7_Batch124A/sequencing_summary_FAO14190_ad60b376.txt nCoV-2019/V2 barcode06
```

Essentially all you will need to do is change every occurrence of **barcode06** (once in guppyplex and three times in minion) to either **barcode07** or **barcode12**. 

**QUESTION** - what is number of mapped reads in each of the samples you have looked at?

As the commands are well structured and all that is needed to change is the input and output names within a run, it means that these commands are easily scriptable using bash (which you have learnt earlier on in this course) or pipeline tools such as [snakemake](https://snakemake.github.io) and [nextflow](https://www.nextflow.io).

At the end of this session we should deactivate our conda environment:

```
conda deactivate
```

### 1.1.4: Additional Resources


We were initially planning on using the [ncov2019-artic-nf](https://github.com/connor-lab/ncov2019-artic-nf) nextflow pipeline for running the ARTIC network tools as it offers a number of extra functions and easier automation. However, we could not get it working on the VM correctly (DSL error resolved by replacing nextflow.preview.dsl = 2 with nextflow.enable.dsl = 2 within the scripts, but the following error could not be resolved by dropping down Nextflow versions (GitHub issue to be created):

```
Unqualified output value declaration has been deprecated - replace `tuple sampleName,..` with `tuple val(sampleName),..`
```

We did also try using the [nf-core/viral-recon](https://nf-co.re/viralrecon) viral nextflow pipeline, which is highly recommended (not just for SARS-CoV-2), but performance was not great on our small VM.

There is a very nice [SARS-CoV-2 sequencing data analysis tutorial](https://github.com/cambiotraining/sars-cov-2-genomics) created by CamBioTraining, which includes [nf-core/viral-recon](https://nf-co.re/viralrecon) instructions.

This MinION tutorial is based on the original [ARTIC](https://artic.network) network's [nCoV-2019 bioinformatics protocol](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html).

## 1.2: Illumina Tutorial

In the above tutorial, we were working with MinION data. In this tutorial we will be using paired end Illumina SARS-CoV-2 data. Although this uses different computational tools, the two approaches are in essence very similar.

### 1.2.1: Setup and data

We do not need to use a conda environment as all the tools we need are already installed directly on the VM: 

* [trim_galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) for read trimming (the reads are pre-trimmed using trim_galore)
* [bwa](https://github.com/lh3/bwa) for read alignment
* [samtools](http://www.htslib.org) for SAM/BAM conversion
* [ivar](https://andersen-lab.github.io/ivar/html/manualpage.html) for primer trimming and consensus calling
* [weeSAM](https://github.com/centre-for-virus-research/weeSAM) for coverage plots.

First, lets move into the Illumina data directory:

```
cd ~/SARS-CoV-2/Illumina/200703_M01569_0148_000000000-J53HN_Batch70/
```

The data in this folder is from a run on an Illumina MiSeq machine. The name of the folder implies it was run on the 3rd July 2020 (200703), the machine ID is M01569, the run ID is 0148_000000000-J53HN, and this was called Batch70 locally within the [Medical Research Council-University of Glasgow Centre for Virus Research](https://www.gla.ac.uk/research/az/cvr/) (CVR) as part of a Covid-19 Genomics UK Consortium ([COG-UK](https://www.cogconsortium.uk)) sequencing run. The samples were sequenced using Version 1 (V1) of the ARTIC [nCoV-2019](https://github.com/artic-network/primer-schemes/tree/master/nCoV-2019) amplicon primers.

There are four samples in this run called:

* CVR2058
* CVR2078
* CVR2092
* CVR2101

If you list the contents of the directory you should see paired end reads (R1.fastq and R2.fastq) for each of the four samples:

```
ls
```

These samples have already been trimmed/filtered using [trim_galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/), so we do not need to QC them. For information, this is command that was used for each sample:

```
#DO NOT ENTER THE BELOW COMMAND- IT IS JUST FOR INFORMATION 
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

You only need to index a genome once, if you are aligning many samples to the same genome sequence (as we will be on this course) you do not need to re-run the reference index step; but don't confuse this step with the BAM indexing step which does need to be done for each sample.

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

We now have a SAM file, which we immediately convert into a BAM file as they are compressed and faster to work with downstream. We can sort the file at the same time as converting it:


```
samtools sort -@4 CVR2058.sam -o CVR2058.bam
```
Breaking this command down:

* **samtools** = the name of the program
* **sort** = the name of the function within samtools to use
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
**NB:** Use **TAB Completion** to help enter the -b primer bed file!

Breaking this command down:

* **ivar** = the name of the program
* **trim** = the name of the function within ivar to use (trimming primers)
* **-i CVR2058.bam** = the name of the input BAM fie (it is in the current directory)
* **-b ~/artic-ncov2019/primer_schemes/nCoV-2019/V1/nCoV-2019.bed** = the path to primer BED file
* **-p CVR2058_trim.bam** = the name of the output file to create (which will be a primer trimmed version of the input BAM file)

For those who are interested ivar works by soft clipping the primer sequences from the read alignments in the BAM file (rather than actually trimming the reads sequences) based on the alignment coordinates. Soft clipping is a form of trimming that is embedded within the CIGAR string of a read alignment. The CIGAR string is the 6th field of the SAM/BAM file, if you were to examine the BAM file manually you should see lots of 'S' characters in the CIGAR field: see the [CIGAR specification](https://samtools.github.io/hts-specs/SAMv1.pdf) for more information.

We now need to sort and then index this BAM file. Sort:

```
samtools sort -@4 CVR2058_trim.bam -o CVR2058_trim_sort.bam 
```

Rename the file back to CVR2058\_trim.bam as we don't need the unsorted BAM file:

```
mv CVR2058_trim_sort.bam CVR2058_trim.bam
```

Index the BAM:

```
samtools index CVR2058_trim.bam
```

We now have a sorted and indexed BAM file (CVR2058\_trim.bam) that contains our sample's paired end reads (CVR2058\_R1.fastq CVR2058\_R2.fastq) aligned to the Wuhan-Hu-1 (MN908947) reference genome, with the amplicon primer sequences clipped off. So we are now ready to call a consensus sequence:

```
samtools mpileup -aa -A -d 0 -Q 0 CVR2058_trim.bam | ivar consensus -p CVR0258 -t 0.4
```

Breaking this command down, there are two parts:

1. samtools [mpileup](http://www.htslib.org/doc/samtools-mpileup.html) which essentially outputs the base and indel counts for each genome position
	* **-aa** = output data for absolutely all positions (even zero coverage ones)
	* **-A** = count orphan reads (reads whose pair did not map)
	* **-d 0** = override the maximum depth (default is 8000 which is typically too low for viruses)
	* **-Q 0** = minimum base quality, 0 essentially means all the data
2. ivar [consensus](https://andersen-lab.github.io/ivar/html/manualpage.html) - this calls the consensus - the output of the samtools mpileup command is piped '|' directly into ivar
	* -p CVR2058 = prefix with which to name the output file
	* -t 0.4 = the minimum frequency threshold that a base must match to be used in calling the consensus base at a position. In this case, an ambiguity code will be used if more than one base is > 40% (0.4). See [ivar manual]

By default, ivar consensus uses a minimum depth (-m) of 10 and a minimum base quality (-q) of 20 to call the consensus; these defaults can be changed by using the appropriate arguments. If a genome position has a depth less than the minimum, an 'N' base will be used in the consensus sequence by default.

ivar will output some basic statistics to the screen such as:

```
#DO NOT ENTER THIS - IT IS THE IVAR OUTPUT YOU SHOULD SEE
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

### 1.2.3: Exercise: Generating Illumina consensus sequences yourself

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
samtools index CVR2058_trim.bam
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

In this session, we will be focussing on analysis steps you can do after you have created a SARS-CoV-2 consensus sequence:

* Determining the lineage of the viral genome sequence using [Pangolin](https://cov-lineages.org/index.html)
* Calling the mutations present within the viral genome sequence using [SPEAR](https://github.com/m-crown/SPEAR)

### 2.1: Pangolin lineages

A **lineage** can be defined as a group of closely related viruses with a common ancestor. SARS-CoV-2 has now evolved into many different lineages, some of which are known as a Variant of Concern (VOC) or Variant Under Investigation (VUI), and some have World Health Organisation (WHO) labels. 

In this tutorial we will be focussing on lineages defined by the [Pangolin](https://cov-lineages.org/index.html) tool. ([Pangolin](https://cov-lineages.org/index.html)) stands for **P**hylogenetic **A**ssignment of **N**amed **G**lobal **O**utbreak **LIN**eages.  Although [NextStrain](https://nextstrain.org/blog/2020-06-02-SARSCoV2-clade-naming) and [GISAID](https://gisaid.org/resources/statements-clarifications/clade-and-lineage-nomenclature-aids-in-genomic-epidemiology-of-active-hcov-19-viruses/) maintain their own clade and lineage nomenclatures we will simply be focussing on Pangolin in this tutorial.

Some of the most well known lineages are (sorted by lineage number):

| Lineage | Shortcut Lineage | WHO Label |
| --- | --- | ---|
| [B.1.1.7](https://cov-lineages.org/lineage.html?lineage=B.1.1.7) | | Alpha | 
| [B.1.1.28.1](https://cov-lineages.org/lineage.html?lineage=P.1) | P.1 | Gamma |
| [B.1.1.28.2](https://cov-lineages.org/lineage.html?lineage=P.2) | P.2 | Zeta |
| [B.1.1.28.1.3](https://cov-lineages.org/lineage.html?lineage=P.3) | P.3 | Theta |
| [B.1.1.529.1](https://cov-lineages.org/lineage.html?lineage=BA.1) | BA.1 | Omicron |
| [B.1.1.529.2](https://cov-lineages.org/lineage.html?lineage=BA.2) | BA.2 | Omicron |
| [B.1.1.529.3](https://cov-lineages.org/lineage.html?lineage=BA.3) | BA.3 | Omicron |
| [B.1.1.529.4](https://cov-lineages.org/lineage.html?lineage=BA.4) | BA.4 | Omicron |
| [B.1.1.529.5](https://cov-lineages.org/lineage.html?lineage=BA.5) | BA.5 | Omicron |
| [B.1.351](https://cov-lineages.org/lineage.html?lineage=B.1.351) | | Beta |
| [B.1.427](https://cov-lineages.org/lineage.html?lineage=B.1.427) | | Epsilon |
| [B.1.429](https://cov-lineages.org/lineage.html?lineage=B.1.429) | | Epsilon |
| [B.1.525](https://cov-lineages.org/lineage.html?lineage=B.1.525) | | Eta |
| [B.1.526](https://cov-lineages.org/lineage.html?lineage=B.1.526) | | Iota |
| [B.1.617.1](https://cov-lineages.org/lineage.html?lineage=B.1.617.1) | | Kappa |
| [B.1.617.2](https://cov-lineages.org/lineage.html?lineage=B.1.617.2) | | Delta |
| [B.1.621](https://cov-lineages.org/lineage.html?lineage=B.1.621) | | Mu |
| [C.37](https://cov-lineages.org/lineage.html?lineage=C.37) | | Lamda |
| [XD](https://cov-lineages.org/lineage.html?lineage=XD) | | Deltacron (not a WHO label) |


When you have obtained a consensus sequence for a sample, determining it's lineage is a typical analysis step. To do this we will use the [Pangolin](https://github.com/cov-lineages/pangolin) command line tool (although there is also a [Pangolin web tool](https://pangolin.cog-uk.io) available).

The pangolin conda environment has already been installed on the VM, but we do need to activate it when we want to use it:

```
conda activate pangolin
```

Now lets move into the folder where we will work:

```
cd ~/SARS-CoV-2/Variants/
```

If you list the contents of this directory you should see a file called ```variant_seqs.fasta``` which contains 19 sequences from Scotland from 2021-2022, you can see their names by typing:

```
grep ">" variant_seqs.fasta 
```

To run Pangolin to get a lineage assignment on these sequences we simply type:

```
pangolin variant_seqs.fasta
```

It is important to remember that a pangolin lineage assignment is a **best guess** at what the lineage of a sequence may be based on available data. Once finished this will create a file called ```lineage_report.csv```. This is a comma separated file which can simply view using ```more``` or ```less``` or ```cat``` etc. The file contains one line of output per input sequence, and each line consists of a number of fields (full descriptions can be found [here](https://cov-lineages.org/resources/pangolin/output.html)):

1. **taxon:** the name of sequence
2. **lineage:** the **most likely** linage assigned by pangolin
3. **conflict:** if the number (N) in this field is > 0, then it means the sequence could have fitted into N different lineage categories. 
4. **ambiguity_score:** SARS-CoV-2 amplicon generated genomes can frequently have failed amplicons creating tracts of Ns in the genome sequence. The ambiguity score is a function of the quantity of missing data at relevant positions in a sequence i.e. lineage defining mutation positions. It represents the proportion of relevant sites in a sequnece which were imputed to the reference values. A score of 1 indicates that no sites were imputed, while a score of 0 indicates that more sites were imputed than were not imputed. This score only includes sites which are used by the decision tree to classify a sequence.
5. **scorpio_call:** [Scorpio](https://github.com/cov-lineages/scorpio) is a tool that assigns a constellation to the sequence. A [constellation](https://github.com/cov-lineages/constellations) is a collection of mutations which are functionally meaningful, but which may have arisen independently a number of times. One of the uses of Scorpio constellations is to define [lineages of concern](https://cov-lineages.org/constellations.html), and if appropriate assign the input sequence to one e.g. 'Omicron (BA.5-like)'.
6. **scorpio_support:** see [Pangolin documentation](https://cov-lineages.org/resources/pangolin/output.html))
7. **scorpio_conflict:** see [Pangolin documentation](https://cov-lineages.org/resources/pangolin/output.html))
8. **scorpio_notes:** see [Pangolin documentation](https://cov-lineages.org/resources/pangolin/output.html))
9. **version:** see [Pangolin documentation](https://cov-lineages.org/resources/pangolin/output.html))
10. **pangolin_version:** see [Pangolin documentation](https://cov-lineages.org/resources/pangolin/output.html))
11. **scorpio_version:** see [Pangolin documentation](https://cov-lineages.org/resources/pangolin/output.html))
12. **constellation_version:** see [Pangolin documentation](https://cov-lineages.org/resources/pangolin/output.html))
13. **is_designated:** see [Pangolin documentation](https://cov-lineages.org/resources/pangolin/output.html))
14. **qc_status:** whether the sequence passed the min length and max N base thresholds
15. **qc_notes:** see [Pangolin documentation](https://cov-lineages.org/resources/pangolin/output.html))
16. **note:** see [Pangolin documentation](https://cov-lineages.org/resources/pangolin/output.html))

As there are many fields in the file, we will ```cut``` the first 5 columns out to ease reading:

```
cut -f1-5 -d ',' lineage_report.csv 
```

**Question:** what is the pangolin lineage assignment for each of the sequences in the ```variant_seqs.fasta``` file? The sequence names in the input file were:

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

As we have now finished with Pangolin we should deactivate the Conda environment:

```
conda deactivate
```

To finish, one of the best places to keep an eye on potentially new lineages is the Pangolin [cov-lienages GitHub issues](https://github.com/cov-lineages/pango-designation/issues) web page.

Other potentially useful websites are:

* [SARS-CoV-2 Variant definitions](https://github.com/phe-genomics/variant_definitions)
* [WHO SARS-CoV-2 Variants](https://www.who.int/en/activities/tracking-SARS-CoV-2-variants/)
* [European Centre for Disease Control - Variants of Concern](https://www.ecdc.europa.eu/en/covid-19/variants-concern)
* [cov-lineages lineage list](https://cov-lineages.org/lineage_list.html)
* [Pangolin Lineage aliases](https://github.com/cov-lineages/pango-designation/blob/master/pango_designation/alias_key.json)
* [Outbreak.info](https://outbreak.info)
* [Virological](http://virological.org/)
* [NextStrain-ncov](https://nextstrain.org/ncov)
* [COV-GLUE](http://cov-glue.cvr.gla.ac.uk/)
* [COG-UK Mutation Explorer](http://sars2.cvr.gla.ac.uk/cog-uk/)
* [ARTIC Network Field Bioinformatics](https://github.com/artic-network/fieldbioinformatics)


### 2.2: SPEAR mutations

[SPEAR](https://github.com/m-crown/SPEAR) stands for **S**ystematic **P**rot**E**in **A**nnotato**R**, it is a mutation calling and annotation tool for SARS-CoV-2 genomes that provides a comprehensive annotation of all the protein products, in particular, Spike (S) mutations are annotated with a range of scores that provide indications of their likely effects on ACE2 binding, and likely contribution to immune escape.

SPEAR has different modes and can run on VCF, alignment, or genome consensus sequence files.  We will be using the consensus option to run SPEAR on the SARS-CoV-2 consensus sequences used above. SPEAR has been installed on the VM, but we need to activate it's Conda environment when we want to use it:

```
conda activate spear
```

Unfortunately, SPEAR does not seem to work on a multi-sequence FASTA file so we will need to split the ```variant_seqs.fasta``` into individual sequence FASTA files. I wrote a simple BASH script called **seq_splitter.sh** to do this which you can download from the [course GitHub repository](https://github.com/WCSCourses/ViralBioinfAsia2022/blob/main/course_data/SARS-CoV-2_workflows/seq_splitter.sh) using ```wget```:

```
wget https://tinyurl.com/vgba2022/seq_splitter.sh
```
**NB:** I made a shortcut URL using [tinyurl](https://tinyurl.com), the full URL was https://raw.githubusercontent.com/WCSCourses/ViralBioinfAsia2022/main/course\_data/SARS-CoV-2\_workflows/seq\_splitter.sh

We now run the BASH script on the ```variants_seqs.fasta``` file:

```
bash seq_splitter.sh variants_seqs.fasta
```

If you list the contents of the directory you should now see lots of ```.fa``` files.

```
ls
```

We are not ready to use SPEAR. We can now run SPEAR on all the ```.fa``` files in the current directory using this command:

```
spear consensus --pangolin=none Scotland_CAMC-14DE972_2021.fa spear_camc
```

**NB:** Remember **TAB Completion** will make entering the sequence filename easy.

Breaking this command down:

* **spear:** the name of the program
* **consensus:** the name of the function within the program to use
* **--pangolin=none:** for this tutorial we are telling SPEAR not to run pangolin to lineage assign the sequences this is because it will check for pangolin updates which may take a while to download/upload on our VM. I would normally not include this option.
* **Scotland\_CAMC-14DE972\_2021.fa:** the name of the input FASTA sequence to analyse
* **spear_camc:** the name of the output file to store the results

SPEAR will output colourful summary tables to the terminal whilst it is running. When it is finished we can examine some of it's outputs in more detail. 

First the ```qc.csv``` file contains basic Quality Control information on the input sequence(s). 

```
more spear_camc/qc.csv
```
The file contains these fields

* **sample\_id:** the sequence name
* **global\_n:** the percentage of Ns in the whole sequence
* **s\_n:** the percentage of Ns in the Spike gene
* **s\_n\_contig:** the longest tract of Ns in the Spike gene
* **rbd\_n:** the longest tract of Ns in the Receptor Binding Domain (RBD) of Spike

It is important to check the quality of your sequences - a high proportion of Ns in a sequence will lead to uncertain lineage assignments and hinder downstream epidemiological investigations. In addition, a mutation can not be called if the codon contains Ns, leading to the old adage that absence of evidence is not evidence of absence.

SPEAR also output a VCF file for each sequence in the ```final_vcfs``` folder of the output folder, as well as a detailed annotation file: 

```
more spear_camc/spear_annotation_summary.tsv 
```

For a full definition of the fields see the [SPEAR documentation](https://github.com/m-crown/SPEAR) under the **Column headings for annotation files** section. If we ```cut``` some of the most important fields to ease visibility:

```
cut -f 1-5,7,9,11 spear_camc/spear_annotation_summary.tsv 
```
This cuts the following fields (columns):

* **1 sample_id:** the name of the sequence
* **2 POS:** the genome position of the variant
* **3 REF:** reference genome base(s)
* **4 ALT:** the alternative base(s) i.e. the mutation/indel
* **5 Gene_Name:** ORF1ab, S, ORF3a, E, M, ORF6, ORF7a, ORF7b, ORF8, N, variants between ORFs are flagged as intergenic.
* **7 consequence_type:** e.g. missense_variant, synonymous_variant
* **9 description:** Free text description of product, ORF1ab is further broken down here into nsp's	
* **11 residues:** Amino acid residue changes in the format of N501Y (REF-AA residue ALT-AA)

For our sample **Scotland\_CAMC\_14DE972\_2021** (which Pangolin assigned to the B.1.526 Iota lineage) we should see the following mutations in Spike:

* L5F
* T95I
* D253G
* E484K
* D614G
* A701V

An alternative is to open the ```spear_annotation_summary.tsv``` file in **LibreOffice-Calc** which is an Excel like program. Navigate to the output folder, right click on the file, select Open with LibreOffice Calc, click OK (it should automatically of select Tab as the delimiter).

We can simply run SPEAR on all the ```.fa``` files in the current directory ```.``` by using this command:

```
spear consensus --pangolin=none . spear_output
```

Although we can examine the same files as above (qc, annotation summary etc), we can also open up the HTML report that SPEAR creates to analyse the samples collectively. This is useful for comparing all the sequences on a single run, or if you have collated a set of related sequences together:


```
firefox spear_output/report/report.html 
```

**Question:** what percentage of the samples contain the N501Y mutation in Spike (surface glycoprotein)? 

### 2.3: Extra work

If you have reached this point and have time to spare you could:

1. Try running pangolin and/or SPEAR on the consensus sequences generated in the MinION and Illumina practicals
2. Try using the [Pangolin web tool](https://pangolin.cog-uk.io) rather than the command line

## 3: SARS-CoV-2 Phylogenetics

In the previous sessions of this course you have learnt how to create sequence alignments and construct phylogenetic trees, using tools such as mafft and iqtree, these tools can (and are) readily applied to SARS-CoV-2 datasets. However, this session introduces a couple of additional phylogenetic tools that are available for SARS-CoV-2 (although though in general they are can be applied to other viruses), one focussed on placing your query genome sequences in the context of the millions already available and another focussed on epidemiological investigations. 

### 3.1: UShER

UShER stands for **U**ltrafast **S**ample placement on **E**xisting t**R**ees. UShER is a program for rapid, accurate placement of samples to existing phylogenies. It essentially enables real-time phylogenetic analyses of the SARS-CoV-2 pandemic which currently (August 2022) has over 11 million genome sequences publicly available.

* [UShER paper](https://www.nature.com/articles/s41588-021-00862-7)
* [UShER manual](https://usher-wiki.readthedocs.io/en/latest/)
* [UShER GitHub](https://github.com/yatisht/usher)
* [UShER Web Tool](https://genome.ucsc.edu/cgi-bin/hgPhyloPlace)

UShER works by adding your sequence onto an exsiting phylogenetic tree. It is rapid because it does not need to re-build the phylogenetic tree from scratch - it is essentially just adding your sequence(s) to certain clusters within the existing tree (the clusters that contain the closest sequences in terms of mutations to the query sequence). However, creating an existing phylogenetic tree of millions of SARS-CoV-2 sequences to add to is by no means trivial - but UShER provide a daily update of all fully public SARS-CoV-2 sequences available for download [here](http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/). 

**NB:** An important point here is that GISAID sequences are NOT fully public due to the GISAID license agreements to protect the rights of sequence submitters. Therefore, this 'fully public' UShER tree only contains sequences from GenBank, COG-UK, and the China National Center for Bioinformation - although this does represent around two thirds of all sequences available. 

In this session, we will be using the [UShER Web Tool](https://genome.ucsc.edu/cgi-bin/hgPhyloPlace) (part of the [UCSC Genome Browser](https://genome.ucsc.edu/util.html)) to phylogenetically place  our sequences into the global tree of the SARS-CoV-2 pandemic. One advantage of this web tool is that it can utilise a phylogenetic tree which includes all available SARS-CoV-2 genomes (including those found only on GISAID).

* On the VM use the Firefox browser to vist the [UShER Web Tool](https://genome.ucsc.edu/cgi-bin/hgPhyloPlace).
* Click the **Browse** button at the top (depending on the browser, this is sometimes called **Choose File**)
* Use the pop-up window to navigate to the ~/SARS-CoV-2/Variants folder and select the **Scotland\_CAMC-14DE972\_2021.fa** file
* Under the 'Phylogenetic tree version' select the GISAID tree
* Click the **Upload** button
* **Wait for UShER to run** (it may take a while!)

***

![](https://github.com/WCSCourses/ViralBioinfAsia2022/blob/main/course_data/SARS-CoV-2_workflows/usher.png)
**Figure 3.1.1:** UShER Web Tool Screenshot: [https://genome.ucsc.edu/cgi-bin/hgPhyloPlace](https://genome.ucsc.edu/cgi-bin/hgPhyloPlace).

***

When UShER has finished it will automatically redirect to a results page. A screenshot of the results you should obtain for the sequence Scotland\_CAMC-14DE972\_2021.fa is below.


***

![](https://github.com/WCSCourses/ViralBioinfAsia2022/blob/main/course_data/SARS-CoV-2_workflows/usher2.png)

**Figure 3.1.2:** UShER results page for genome sequence Scotland\_CAMC-14DE972\_2021.fa 
***

The UShER output contains alot of useful information such as lineage assignments and mutations for each sequence analysed. To actually examine the cluster our sequence has been placed in within the glonal tree we can follow the automatic links to [NextStrain](https://nextstrain.org). You should see a tree similar to this (I changed the "Color By" option to "Pango lineage assigned by UShER").

**NB:** As the global tree is updated daily, there is a chance the tree you get may be slightly different depending to the one belon when it is run.

***

![](https://github.com/WCSCourses/ViralBioinfAsia2022/blob/main/course_data/SARS-CoV-2_workflows/nextstrain.png)
**Figure 3.1.3:** NextStrain visualisation of subtree for sequence Scotland\_CAMC-14DE972\_2021.fa. The colour scheme was changed using the "Color By" option to "Pango lineage assigned by UShER". Our uploaded seqeunce is highlighted in yellow. Note that our sequence was actually already in the gloabl tree - you can see the original Scotland/CAMC-14DE972/2021 sequence about 10 nodes above our uploaded on. 

***


As UShER may take a while to run (because this is a free and public resource), I have prerun each sequence and the NextStrain view of the subtree's can be found here:

* [Scotland/LSPA-3E62463/2022](https://nextstrain.org/fetch/genome.ucsc.edu/trash/ct/subtreeAuspice1_genome_12d2a_5fd460.json)
* [Scotland/LSPA-3E3B700/2022](https://nextstrain.org/fetch/genome.ucsc.edu/trash/ct/subtreeAuspice2_genome_12d2a_5fd460.json)
* [Scotland/QEUH-3DC59CA/2022](https://nextstrain.org/fetch/genome.ucsc.edu/trash/ct/subtreeAuspice3_genome_12d2a_5fd460.json)
* [Scotland/QEUH-3D90306/2022](https://nextstrain.org/fetch/genome.ucsc.edu/trash/ct/subtreeAuspice4_genome_12d2a_5fd460.json)
* [Scotland/QEUH-3D43C43/2022](https://nextstrain.org/fetch/genome.ucsc.edu/trash/ct/subtreeAuspice5_genome_12d2a_5fd460.json)
* [Scotland/CVR14531/2022](https://nextstrain.org/fetch/genome.ucsc.edu/trash/ct/subtreeAuspice6_genome_12d2a_5fd460.json)
* [Scotland/QEUH-36897FC/2022](https://nextstrain.org/fetch/genome.ucsc.edu/trash/ct/subtreeAuspice7_genome_12d2a_5fd460.json)
* [Scotland/QEUH-36491AF/2022](https://nextstrain.org/fetch/genome.ucsc.edu/trash/ct/subtreeAuspice8_genome_12d2a_5fd460.json)
* [Scotland/QEUH-2D86BD1/2021](https://nextstrain.org/fetch/genome.ucsc.edu/trash/ct/subtreeAuspice9_genome_12d2a_5fd460.json)
* [Scotland/QEUH-2D7F704/2021](https://nextstrain.org/fetch/genome.ucsc.edu/trash/ct/subtreeAuspice10_genome_12d2a_5fd460.json)
* [Scotland/QEUH-1BA3933/2021](https://nextstrain.org/fetch/genome.ucsc.edu/trash/ct/subtreeAuspice11_genome_12d2a_5fd460.json)
* [Scotland/QEUH-1B0246E/2021](https://nextstrain.org/fetch/genome.ucsc.edu/trash/ct/subtreeAuspice12_genome_12d2a_5fd460.json)
* [Scotland/QEUH-1725CCB/2021](https://nextstrain.org/fetch/genome.ucsc.edu/trash/ct/subtreeAuspice13_genome_12d2a_5fd460.json)
* [Scotland/QEUH-1585B0A/2021](https://nextstrain.org/fetch/genome.ucsc.edu/trash/ct/subtreeAuspice14_genome_12d2a_5fd460.json)
* [Scotland/QEUH-158D786/2021](https://nextstrain.org/fetch/genome.ucsc.edu/trash/ct/subtreeAuspice15_genome_12d2a_5fd460.json)
* [Scotland/CAMC-14DE972/2021](https://nextstrain.org/fetch/genome.ucsc.edu/trash/ct/subtreeAuspice16_genome_12d2a_5fd460.json)
* [Scotland/QEUH-147E6F5/2021](https://nextstrain.org/fetch/genome.ucsc.edu/trash/ct/subtreeAuspice17_genome_12d2a_5fd460.json)
* [Scotland/QEUH-13ADEF6/2021](https://nextstrain.org/fetch/genome.ucsc.edu/trash/ct/subtreeAuspice18_genome_12d2a_5fd460.json)
* [Scotland/QEUH-138F944/2021](https://nextstrain.org/fetch/genome.ucsc.edu/trash/ct/subtreeAuspice19_genome_12d2a_5fd460.json)

The overall UShER report on all sequences is available [here]()

If you have a set of potentially related sequences you would upload them together into an UShER run and then examine if they are being assigned to the same subtree cluster by UShER - this is similar to what we will be doing next with CIVET.

[Old Link](https://nextstrain.org/fetch/genome.ucsc.edu/trash/ct/subtreeAuspice1_genome_1dd2a_56d490.json)

### 3.2: CIVET

[CIVET](https://github.com/artic-network/civet) stands for **C**luster **I**nvestigation & **V**irus **E**pidemiology **T**ool, its full documentation is available [here](https://cov-lineages.org/resources/civet.html). CIVET was developed to aid SARS-CoV-2 outbreak investigations, it puts new sequences into the context of known background diversity, and can summarise the background diversity based on the users input.

We will using CIVET to investigate a small cluster of 7 hypothetic SARS-CoV-2 sequences. Lets imagine these have all been obtained from a single epidemiological incident such as an outbreak at a school, and we are wanting to examine how related the sequences are to one another and also place them all into the context of the background sequences **available** at the time. All 7 hypothetical sequences are from Scotland and Pangolin assigned them to linage BA.5.3:

* This lineage was simply chosen as an example of a current Omicron lineage with a small number of sequences (89 in total) to ease the burden of running things on the VM
* May and June 2022 as those were the most recent months (downloaded 7th July 2022)
* Scotland because the [MRC-University of Glasgow Centre for Virus Research](https://www.gla.ac.uk/research/az/cvr/) is in Scotland
* Importantly, both the sequence and metadata is publicly available to download from the [COG-UK website](https://www.cogconsortium.uk/priority-areas/data-linkage-analysis/public-data-analysis/).
* [cov-lineages BA.5.3 link](https://cov-lineages.org/lineage.html?lineage=BA.5.3)
* [outbreak.info BA.5.3 link](https://outbreak.info/situation-reports?pango=BA.5.3)

Here is a simple phylogenetic tree of all 89 BA.5.3 samples (created using iqtree2), and using the aligned sequences downloaded from the [COG-UK website](https://www.cogconsortium.uk/priority-areas/data-linkage-analysis/public-data-analysis/):

![](https://github.com/WCSCourses/ViralBioinfAsia2022/blob/main/course_data/SARS-CoV-2_workflows/ba53_tree.png)

**Figure 3.2.1:** Phylogenetic tree of Scottish BA.5.3 sequences from May and June 2022. A single BA.3 sequence (Scotland/QEUH-3E1D6AD/2022) is used to the root the tree. 

First lets create a directory to work in

```
mkdir ~/SARS-CoV-2/Phylo/civet
```

and then move into in

```
cd ~/SARS-CoV-2/Phylo/civet
```

and then download the data using the ```wget``` command:

```
wget https://tinyurl.com/vgba2022/civet_data.tar
```

Now untar (unpack) the data

```
tar -xvf civet_data.tar
```

If if we list the contents of the directory we should see a **samples.fasta** and **samples_metadata.csv** file, along with a folder called **BA53**.

```
ls
```

We have 7 SARS-CoV-2 genome sequences to investigates, these are simply called Sample1-7:

```
grep ">" samples.fasta
```

The **samples_metadata.csv** file is a simple comma separate file (csv) with two columns:

1. **name:** the name of the sample - must match the names used in the FASTA file
2. **sample_date:** the date the sample was taken in YYYY-MM-DD format

This is the minimal metadata needed for CIVET to run. If you were using your own data and don't know the sample date you could simply estimate it. Many additional columns can be added to the metadata file and used for tip labels/colours in the tree: see the [CIVET documentation](https://cov-lineages.org/resources/civet.html) for full details.

The CIVET Conda environment has already been installed on the VM, but as before we first need to activate it when we want to use it:

```
conda activate civet
```

To run a CIVET analysis we need:

* background data alignment file: **BA53/ba.5.3_aligned.fasta** 
* background data metadata: **BA53/ba.5.3_metadata.csv**
* input sample sequences: **samples.fasta**
* input sample metadata: **samples_metadata.csv**

To create the CIVET report we need to type this command

```
civet -i samples_metadata.csv -f samples.fasta -d BA53 -o civet_output
```

Breaking this command down:

* **civet:** the name of the program
* **-i samples_metadata.csv**: the input metadata file in csv format
* **-f samples.fasta**: the input sequence file in FASTA format
* **-d BA53**: the background data directory - with the background sequence alignment and metadata
* **-o civet_output**: the name of the folder to create and output the results into

After finishing CIVET should create the output folder and results:

```
ls civet_output
```

Let's open up the main CIVET report file:

```
firefox civet_output/civet.html
```
**NB:** You could navigate to the output via the Files browser and double click to open

**CIVET Report Table 1**

This is a summary of the input sequences and their input data along with what **catchment** tree the sequence is contained in - all our samples should be within a single tree called **catchment_1**.

Essentially, CIVET searches through the background sequences and finds all sequences that fall within a customisable SNP distance of each input sequence (I believe the default is 2 down and 2 up i.e. sequences that have all mutations and up to two more, and those that have  up to 2 of the mutations missing). All equally distant targets are included in the catchment. For a given input sequence, if no background sequences are found within the SNP distance cut off, the CIVET algorithm increases the SNP distance in all directions and attempts to get at least one sequence per category (up, down or side). This results in a set of background sequences for each input sequence, and any input sequences with overlapping targets have their catchments merged together.

**CIVET Report Table 2**

This is a table of tall the sequences that passed the CIVET quality control (QC) checks such as minimum length, and maximum N content.

**CIVET Catchments**

The report then contains details of each catchment find - our report only contains a single catchment. First there is a table summary of the catchment including information on the number of input sequences that are assigned to it and the earliest and latest dates of the background sequences in the catchment (from the background metadata file).

One of the main results of the CIVET report is the catchment tree. We can expand the tree to ease visualisation by:
* Clicking the + symbol on the Tree Option line
* Sliding the Expansion slider along to expand the gaps between sequences in the tree

You should see a tree like this:

![](https://github.com/WCSCourses/ViralBioinfAsia2022/blob/main/course_data/SARS-CoV-2_workflows/civet_tree.png)

**Figure 3.2.2:** CIVET catchment tree containing all 7 of all hypothetical SARS-CoV-2 input sequences. By default, input sequences are coloured with turquoise circles and sequences from the background data in purple.

It is important to emphasise that this type of analysis is not able to infer direct transmission between two samples; even identical sequences may be unrelated as
SARS-CoV-2 is relatively slow evolving for an RNA virus. Furthermore, the completeness of the background will obviously be an important factor. All our samples form part of the same catchment and are relatively close to one another. Samples 6 and 7 are identical and are a 1 SNP distance from many other samples (e.g.Scotland/LSPA-3E7E789/2022). Sample 5 is identical to the background sequence Scotland/LSPA-3E8210C, whereas sample 1-4 branch of from that, with samples 1 and 2 being identical, and sample 4 having a number of additional mutations.

The nucleotide mutations that each sample contains can be examined and compared using the [snipit](https://github.com/aineniamh/snipit) plot in the CIVET report. As can be seen, Sample 1-4 all contain a T to G mutation at position 3,154 that the other 3 samples do not contain, whilst samples 6 and 7 contain A to G at position 24,595 not seen in the other 5 samples, and sample 4 contains a number of unique mutations:

![](https://github.com/WCSCourses/ViralBioinfAsia2022/blob/main/course_data/SARS-CoV-2_workflows/civet_snipit.png)

**Figure 3.2.3:** snipit plot of the nucleotide mutations contained within each input sequence in the catchment; coloured by base: A (dark blue), C (red), G (light blue), T (green). The x-axis is the genome position the mutation is located, the y-axis are the samples. The grey bar at the bottom represented the SARS-CoV-2 genome and highlights where on the genome each mutation is located.

### 3.3: SARS-CoV-2 Alignments

Just some quick notes on alignment (not part of this session). Aligning millions of SARS-CoV-2 sequences together can present problems computationally (not least the amount of time takes). Another issue is the large N tracts in sequences which have failed ARTIC amplicons, these can cause numerous issues during alignments resulting in spurious indels and an overall poor alignment. 

Different groups have taken different approaches to solve this. COG-UK uses the [grapevine](https://github.com/COG-UK/grapevine) pipeline which utilises minimap2 to align each read individually against the SARS-CoV-2 reference sequence (similar to aligning each read in a FASTQ read to the reference), insertions are then trimmed and the data is outputted from the BAM file to create a FASTA alignment whose length is the size of original reference sequence; an obvious downside here is that you loose all the insertions from the alignment. This is often combined with aggressive filtering of sequences that are too short and have too many N bases.

GISAID

## 4: SARS-CoV-2 Group Practical

In this session, we will be working on some more Illumina paired end read data. The FASTQ data was downloaded from the [European Nucleotide Archive](https://www.ebi.ac.uk/ena/browser) (ENA), and there are 4 samples in total (the samples are not related to one another), with R1 and R2 FASTQ files for each:

* ERR9105817 - ARTIC primer version 3
* ERR9731990 - ARTIC primer version 3
* ERR9761275 - ARTIC primer version 3
* ERR9788433 - ARTIC primer version 3

Your task is to work as a group in the breakout rooms to analyse these samples. Initial read QC (with trim_galore) is not required (but you could add it if you wanted).  You should:

* align the reads to the Wuhan-Hu-1 reference sequence
* Report the number of mapped reads
* Trim the ARTIC primers
* Call a consensus sequence
* Use Pangolin to assign a lineage
* Use SPEAR to call the mutations

This is a flexible session, and a chance to collate all the steps that you have learnt onto a single sample(s).

As a group you could:

* Analyse a sample each and collate the results. As there are only 4 samples (and groups will likely be larger than 4) - multiple people could analyse a single sample and check you get the same results
* Write a bash script to process the sample automatically. Remember all the steps to analyse a sample are the same, it is just the input/output names that are changing. Completed example bash scripts will be uploaded here after the session.


## 5: Warnings

I would consider this VM a good place to learn BUT not necessarily a good place to conduct 'real' analyses. The reason being is that many of the SARS-CoV-2 tools and datasets are updated very frequently which means many will be out of date on the VM already (many of the tools were installed a few months ago). Tools such as Pangolin and SPEAR do however have good update functions.
