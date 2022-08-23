# Coverage plots and statistics

Here we are going to use the software packages [SAMtools](http://samtools.sourceforge.net/) and [Qualimap](http://qualimap.conesalab.org/) to explore how the genome sequence reads are distributed across the viral genome. In other words, we are going to generate a coverage plot. In addition to the coverage plot, the Qualimap software also provides several other reports and metrics that can be used for quality control (QC). You can read about Qualimap in this paper:

 - Okonechnikov K, Conesa A, García-Alcalde F. (2016) Qualimap 2: advanced
   multi-sample quality control for high-throughput sequencing data.
   *Bioinformatics*. **32**:292-4. doi:[10.1093/bioinformatics/btv566](https://doi.org/10.1093/bioinformatics/btv566).

 The SAMtools package contains many really useful tools for analysing and manipulating BAM files. You can read more about SamTools here: 
 
 - Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N.,    Marth, G., Abecasis, G., Durbin, R., & 1000 Genome Project
   Data          Processing Subgroup (2009). The Sequence Alignment/Map
   format and SAMtools. Bioinformatics, **25**, 2078–2079.       
   https://doi.org/10.1093/bioinformatics/btp352.      

Qualimap and SAMtools are already  installed on the virtual machine. However, it you need to install the Qualimap software, you can obtain it from: http://qualimap.conesalab.org/ and SAMtools from  http://samtools.sourceforge.net/.

We can run the Qualimap software either:

 - in a graphical mode or
 - on the command line.

During this workshop activity, we will try both modes of running the software.

We will use SAMtools only on the command line; there is not graphical interface for this software package.

## Input data
Qualimap examines coverage of the reference genome sequence by aligned sequence reads. Therefore, it requires that the reads have already been aligned against the reference genome sequence and that the alignment is presented in BAM format.

We already discussed BAM format in the session "NGS file formats and data and QC".  You can find a formal description of the BAM file format here: https://samtools.github.io/hts-specs/SAMv1.pdf.

The specific viral sequence datasets that we will examine are:
 - [ERR8261957](https://www.ebi.ac.uk/ena/browser/view/ERR8261957) (Oxford Nanopore)
 - [ERR8261968](https://www.ebi.ac.uk/ena/browser/view/ERR8261968) (llumina MiSeq)

These datasets are already downloaded and provided to you in the virtual machine. Now, navigate to the appropriate directory with this command:

    cd /home/manager/ViralBioinfAsia2022/course_data/Coverage_Plots_Stats

Now, when you list the contents of the current directory with `ls -lh`, you should see something like this:

    total 16M
    -rw-rw-r-- 1 manager manager 2.9K Jul 25 16:20 ARTIC_amplicons.bed
    -rw-rw-r-- 1 manager manager 8.0M Jul 25 16:20 ERR8261957.bam
    -rw-rw-r-- 1 manager manager  152 Jul 25 16:20 ERR8261957.bam.bai
    -rw-rw-r-- 1 manager manager 7.6M Jul 25 16:20 ERR8261968.bam
    -rw-rw-r-- 1 manager manager  152 Jul 25 16:20 ERR8261968.bam.bai
    drwxrwxr-x 2 manager manager 4.0K Jul 25 16:20 images
    -rw-rw-r-- 1 manager manager 1.2K Jul 25 16:20 readME.md

So, you can see we have the two BAM files, along with their accompanying index files. These BAM files consist of genomic sequence reads aligned against Alignments of sequence reads against the [Wuhan-Hu-1 reference genome]([https://github.com/WCSCourses/ViralBioinfAsia2022/raw/main/course_data/Coverage_Plots_Stats/images/Screenshot%202022-07-04%20at%2016.25.00.png](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2)) sequence.

## SAMtools flagstat

Let's gather some information about these alignments. First, let's use the rather basic *flagstat* tool from the SamTools package:

    samtools flagstat ERR8261957.bam
    samtools flagstat ERR8261968.bam

Examine the output from these two commands. Can you figure out:

 - How many sequence reads are in ERR8261957?
 - How many sequence reads are in ERR8261968?
 - Which of the two datasets contains paired reads?

## Qualimap (on the command line)

Now, let's use Qualimap to gather some information about these same two alignments. Execute the following command lines:

    qualimap bamqc -bam ERR8261957.bam
    qualimap bamqc -bam ERR8261968.bam

Now, if you list the directory contents with `ls -lh`, you will notice new directories called `ERR8261957_stats` and `ERR8261968_stats`. In each of these directories, you will find an HTML file that contains a detailed report on the BAM alignment. 
Execute the command `ls -lh *_stats` to see the contents of these directories. Now, let's open the HTML files in the Firefox web browser:

    firefox *_stats/*.html

You should now see something like this:

![Qualimap reports opened in Firefox web browser](https://github.com/WCSCourses/ViralBioinfAsia2022/raw/main/course_data/Coverage_Plots_Stats/images/Screenshot%202022-07-29%20at%2021.17.12.png)

Note that the web browser contains two tabs; one contains the report for ERR8261957.bam and the other contains report for ERR8261968.bam.

Take a look at the reports to see the range of information provided. Then close the Firefox web browser.
 

## Qualimap (graphical interface)

As an alternative to the command-line, we can launch a graphical interface for Qualimap. 
   
Start the Qualimap graphical interface by executing the command `qualimap &` in the Terminal. Don't worry about the warning message concerning missing R packages; just click the *Ok* button.  You should now see something like this:

![enter image description here](https://github.com/WCSCourses/ViralBioinfAsia2022/raw/main/course_data/Coverage_Plots_Stats/images/Screenshot%202022-07-30%20at%2011.07.42.png)

Notice the File menu item at the top left corner of the Qualimap window. Select *File* -> *New analysis* -> *BAM QC* and then select one of the BAM files. When the analysis is complete you will see the results in the Qualimap window and you can click on the various sections of the results report in the menu on the left:

![enter image description here](https://github.com/WCSCourses/ViralBioinfAsia2022/raw/main/course_data/Coverage_Plots_Stats/images/Screenshot%202022-07-30%20at%2011.12.25.png)

Note that under the File menu, we can export the results report as an HTML or PDF file. 

You can find some good information about interpreting the Qualimap results here: http://qualimap.conesalab.org/doc_html/analysis.html#output.

For our datasets, try to answer these questions:

 - What is the depth of coverage?
 - Is sequencing coverage uniform across the whole genome?
 - How long are the sequence reads?
 - For paired-read data, what is the insert size (i.e. the length of the sequenced fragments)?

### Visualising the data using IGV
If you have any spare time then optionally, let's visually inspect the alignment data using the Integrative Genome Viewer (IGV). This will provide us with an overview of the data. (Alternatively, you can use Tablet instead of IGV, if you pr).

To do this, we need the BAM files, obviously; we also need files containing the reference genome sequence. Lets download those now. First, make sure that we are in the correct directory:

    cd /home/manager/ViralBioinfAsia2022/course_data/Coverage_Plots_Stats

Then use *wget* to download the reference-genome files:

    wget https://github.com/WCSCourses/ViralBioinfAsia2022/raw/main/reference_genome.fasta https://github.com/WCSCourses/ViralBioinfAsia2022/raw/main/reference_genome.fasta.fai https://github.com/WCSCourses/ViralBioinfAsia2022/raw/main/reference_genome.gff3
    --2022-07-30 12:07:19--  https://github.com/WCSCourses/ViralBioinfAsia2022/raw/main/reference_genome.fasta

Now, when you execute `ls -lh` , you will notice we now have the following files: `reference_genome.fasta`,  `reference_genome.fasta.fai` and `reference_genome.gff3`. We are ready to run the IGV genome browser software by executing the command `igv`.

Click on *Genomes* -> *Load genome from file* and select `reference_genome.fasta`.
Click on *File* -> *Load from file* and select `reference_genome.gff3`.
Click on  *File* -> *Load from file* and select  `ERR8261957.bam` and `ERR8261968.bam`. 

You can now explore the two alignments in the IGV browser, which should look something like this:

![Alignments of sequence reads against the Wuhan-Hu-1 reference genome, visualised using IGV.](https://github.com/WCSCourses/ViralBioinfAsia2022/raw/main/course_data/Coverage_Plots_Stats/images/Screenshot%202022-07-04%20at%2016.25.00.png)
Alignments of sequence reads against the Wuhan-Hu-1 reference genome, visualised using IGV.

What do you notice about the patterns of coverage of the reference genome? Are there any differences between the Illumina sequencing data versus the Oxford Nanopore sequencing data?

