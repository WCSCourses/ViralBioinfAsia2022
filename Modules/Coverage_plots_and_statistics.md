# Coverage plots and statistics

Here we are going to use the [Qualimap software package](http://qualimap.conesalab.org/) to explore how the genome sequence reads are distributed across the viral genome. In other words, we are going to generate a coverage plot. In addition to the coverage plot, the Qualimap software also provides several other reports and metrics that can be used for quality control (QC). You can read about Qualimap in this paper:

 - Okonechnikov K, Conesa A, García-Alcalde F. (2016) Qualimap 2: advanced
   multi-sample quality control for high-throughput sequencing data.
   *Bioinformatics*. **32**:292-4. doi:[10.1093/bioinformatics/btv566](https://doi.org/10.1093/bioinformatics/btv566).

Qualimap should already be installed on the virtual machine. However, it you need to install the software, you obtain it from: http://qualimap.conesalab.org/.

We can run the Qualimap software either:

 - in a graphical mode or
 - on the command line.

During this workshop activity, we will try both modes of running the software.

## Input data
Qualimap examines coverage of the reference genome sequence by aligned sequence reads. Therefore, it requires that the reads have already been aligned against the reference genome sequence and that the alignment is presented in BAM format.

We already discussed BAM format in the session "NGS file formats and data and QC".  You can find a formal description of the BAM file format here: https://samtools.github.io/hts-specs/SAMv1.pdf.

The specific viral sequence datasets that we will examine are:
 - [ERR8261957](https://www.ebi.ac.uk/ena/browser/view/ERR8261957) (Oxford Nanopore)
 - [ERR8261968](https://www.ebi.ac.uk/ena/browser/view/ERR8261968) (llumina MiSeq)
 
### Visualising the data using IGV
Before we calculate the QC metrics on out data, first, let's visually inspect it using the Integrative Genome Viewer (IGV). This is will provide us with an overview of the data. Alternatively, we can use Tablet instead of IGV, if you prefer.

![Alignments of sequence reads against the Wuhan-Hu-1 reference genome, visualised using IGV.](https://github.com/WCSCourses/ViralBioinfAsia2022/raw/main/course_data/Coverage_Plots_Stats/images/Screenshot%202022-07-04%20at%2016.25.00.png)
Alignments of sequence reads against the Wuhan-Hu-1 reference genome, visualised using IGV.

What do you notice about the patterns of coverage of the reference genome? Are there any differences between the Illumina sequencing data versus the Oxford Nanopore sequencing data?
