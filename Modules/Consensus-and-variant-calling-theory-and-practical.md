**9 Consensus and variant calling theory and practical**

*9.1 Overview*

![Consensus-image1.png](https://github.com/WCSCourses/ViralBioinfAsia2022/blob/main/Modules/images/Consensus-image1.png)

>*9.2 Creating a consensus with BCFtools*

>Let's start by doing some quick prep/housekeeping:

-----------------------------------------------------------------------
unzip /home/manager/course_data/Consensus_and_variant_calling/Consensus_and_variant_calling.zip
-----------------------------------------------------------------------
mv /home/manager/course_data/Consensus_and_variant_calling/Consensus_and_variant_calling/09-dengue-lofreq/ /home/manager/course_data/Consensus_and_variant_calling/
-----------------------------------------------------------------------
mv /home/manager/course_data/Consensus_and_variant_calling/Consensus_and_variant_calling/09-chikv-lofreq/ /home/manager/course_data/Consensus_and_variant_calling/
-----------------------------------------------------------------------
rm -r /home/manager/course_data/Consensus_and_variant_calling/Consensus_and_variant_calling/
-----------------------------------------------------------------------

>We need the genome files from the Reference_alignment section so let's copy those over:

-----------------------------------------------------------------------
cp /home/manager/course_data/Reference_alignment/07-dengue_align/dengue-genome.fa /home/manager/course_data/Consensus_and_variant_calling/09-dengue-lofreq/
-----------------------------------------------------------------------
cp /home/manager/course_data/Reference_alignment/07-chikv-align/chikv-genome.fasta /home/manager/course_data/Consensus_and_variant_calling/09-chikv-lofreq/
-----------------------------------------------------------------------
>We will also need the annotation files so let's bring those over as well:
-----------------------------------------------------------------------
cp /home/manager/course_data/Reference_alignment/07-dengue_align/annotation.txt /home/manager/course_data/Consensus_and_variant_calling/09-dengue-lofreq/
-----------------------------------------------------------------------
cp /home/manager/course_data/Reference_alignment/07-chikv-align/annotation.txt /home/manager/course_data/Consensus_and_variant_calling/09-chikv-lofreq/
-----------------------------------------------------------------------

>I've already run this analysis for you so you can check your results.  In order to not overwrite these files, let's rename them:

-----------------------------------------------------------------------
mv /home/manager/course_data/Consensus_and_variant_calling/09-dengue-lofreq/dengue-consensus-bcftools.fa /home/manager/course_data/Consensus_and_variant_calling/09-dengue-lofreq/OMS-dengue-consensus-bcftools.fa
-----------------------------------------------------------------------
mv /home/manager/course_data/Consensus_and_variant_calling/09-chikv-lofreq/chikv-consensus-bcftools.fa /home/manager/course_data/Consensus_and_variant_calling/09-chikv-lofreq/OMS-chikv-consensus-bcftools.fa 
-----------------------------------------------------------------------
mv /home/manager/course_data/Consensus_and_variant_calling/09-dengue-lofreq/dengue-direct-low-requenncy-variants.vcf /home/manager/course_data/Consensus_and_variant_calling/09-dengue-lofreq/OMS-dengue-direct-low-requenncy-variants.vcf
-----------------------------------------------------------------------
mv /home/manager/course_data/Consensus_and_variant_calling/09-chikv-lofreq/chikv-direct-low-requenncy-variants.vcf /home/manager/course_data/Consensus_and_variant_calling/09-chikv-lofreq/OMS-chikv-direct-low-requenncy-variants.vcf
-----------------------------------------------------------------------

>Finally, one of the programs we are going to use isn't installed so let's do that now:

-----------------------------------------------------------------------
sudo apt install bcftools
-----------------------------------------------------------------------

>You'll be using the 'sudo' (super user) command to install this so your system will ask you for a password.  Type in your login (which should be 'manager')
>
>Ok, all the prep work is done!  Now let's navigate to the proper folder so we can get started:

-----------------------------------------------------------------------
cd /home/manager/course_data/Consensus_and_variant_calling/09-dengue-lofreq/
-----------------------------------------------------------------------

>Now that we have our alignment file, we can now generate a consensus
>sequence. We will explore two separate ways to generate this consensus.
>In the first example we will use **bcftools mpileup** and **bcftools
>call** to generate a "majority rules" consensus and the output format
>will be in **.fasta** file format. For many downstream applications this
>will be the output you will likely want. In the second example, we will
>also be using use **bcftools mpileup** and **bcftools call** but will
>also add in the **vcfutils.pl vcf2fq** script to generate a **.fastq**
>file that allows for ambiguities at each position if there is a mixture
>of variants present in the sample.

>Example 1 - Using **bcftools** to generate a **.fasta**
>\"majority rules\" consensus file

>We will be using the following options for bcftools mpileup and bcftools
>call respectively:

>**bcfools mpileup**

>-f, \--fasta-ref FILE faidx indexed reference sequence file

>**bcftools call**

>-m, \--multiallelic-caller alternative model for multiallelic and
>rare-variant calling (conflicts with -c)

>-v, \--variants-only output variant sites only

>Here is how we will call this command:

  -----------------------------------------------------------------------
  bcftools mpileup -f dengue-genome.fa dengue-direct_10p.bam | bcftools call -mv -Oz --ploidy-file annotation.txt -o calls.vcf.gz
  -----------------------------------------------------------------------

  -----------------------------------------------------------------------
  tabix calls.vcf.gz
  -----------------------------------------------------------------------

  -----------------------------------------------------------------------
  cat dengue-genome.fa \| bcftools consensus calls.vcf.gz \> dengue-cns.fa
  -----------------------------------------------------------------------

>Finally, we will use a quick hack to rename the header of the .fasta
>file we just generated:

  -----------------------------------------------------------------------
  sed \'s/\>MN566112.1 Dengue virus 2 isolate New Caledonia-2018-AVS127,
  complete genome/\>dengue-consensus-bcftools/\' dengue-cns.fa \>
  dengue-consensus-bcftools.fa && rm dengue-cns.fa
  -----------------------------------------------------------------------

  -----------------------------------------------------------------------

>We created (and removed) several files here. The one that contains the
>consensus genome you are looking for is **dengue-consensus-bcftools.fa**

>Example 2 - Using **samtools/bcftools** along with vcfutils.pl vcf2fq to
>generate a **.fastq** consensus file that allows for ambiquities

>Here we will be using the following options for samtools mpileup and
>bcftools call respectively:

>**bcftools mpileup**

>-f, \--fasta-ref FILE faidx indexed reference sequence file

>**bcftools call**

>-c, \--consensus-caller the original calling method (conflicts with -m)

>Here\'s how we will call this command:

  -----------------------------------------------------------------------
  bcftools mpileup -f dengue-genome.fa dengue-direct_10p.bam \| bcftools call -c
  \--ploidy-file annotation.txt \| vcfutils.pl vcf2fq \>
  dengue-consensus-full.fq
  -----------------------------------------------------------------------

>Finally, we will again use a quick hack to rename the header of your
>output file:

  -----------------------------------------------------------------------
  sed \'s/@MN566112.1/\>dengue-consensus-vcf2fq/\'
  dengue-consensus-full.fq \> dengue-consensus-vcf2fq.fq && rm
  dengue-consensus-full.fq
  -----------------------------------------------------------------------

>We created (and removed) several files here. The one that contains the
>consensus genome you are looking for is **dengue-consensus-vcf2fq.fq**

>*9.3 Visualizing the difference between these files*

>So what effectively is the difference between the
>**dengue-consensus-bcftools.fa** and **dengue-consensus-vcf2fq.fq**
>files we generated? If we were to align these back to the reference
>genome dengue-genome.fa we can see the ambiquity present in the
>**dengue-consensus-vcf2fq.fq** whereas the
>**dengue-consensus-bcftools.fa** merely shows a complete consensus
>change at this position.

![Consensus-image2.png](https://github.com/WCSCourses/ViralBioinfAsia2022/blob/main/Modules/images/Consensus-image2.png)

>*9.4 Variant calling*

>There are a number of different variant callers available including
>FreeBayes, LoFreq, VarDict, and VarScan2. In this tutorial, we will be
>using LoFreq\* (i.e. LoFreq version 2) - a fast and sensitive
>variant-caller for inferring SNVs and indels from next-generation
>sequencing data. LoFreq\* makes full use of base-call qualities and
>other sources of errors inherent in sequencing (e.g. mapping or
>base/indel alignment uncertainty), which are usually ignored by other
>methods or only used for filtering.

>LoFreq\* can run on almost any type of aligned sequencing data (e.g.
>Illumina, IonTorrent or Pacbio) since no machine- or
>sequencing-technology dependent thresholds are used. It automatically
>adapts to changes in coverage and sequencing quality and can therefore
>be applied to a variety of data-sets e.g. viral/quasispecies, bacterial,
>metagenomics or somatic data.

>LoFreq\* is very sensitive; most notably, it is able to predict variants
>below the average base-call quality (i.e. sequencing error rate). Each
>variant call is assigned a p-value which allows for rigorous false
>positive control. Even though it uses no approximations or heuristics,
>it is very efficient due to several runtime optimizations and also
>provides a (pseudo-)parallel implementation. LoFreq\* is generic and
>fast enough to be applied to high-coverage data and large genomes. On a
>single processor it takes a minute to analyze Dengue genome sequencing
>data with nearly 4000X coverage, roughly one hour to call SNVs on a 600X
>coverage E.coli genome and also roughly an hour to run on a 100X
>coverage human exome dataset.

>For more details on the original version of LoFreq see Wilm et al.
(2012).

>LoFreq comes with a variety of subcommands. By just typing lofreq you
>will get a list of available commands. The command we will be using
>today is lofreq call which, will call variants in our dengue alignment
>from the previous section. The output for our variant calls will be in
>variant call format or, vcf.

>*9.5 Running LoFreq\**

>Let\'s start by checking to make sure we are working with the right
>files. Usually this isn\'t necessary if you\'ve just made these files
>but it\'s a good sanity check if you have thousands of samples you are
>working with at any given time:

  -----------------------------------------------------------------------
  lofreq checkref dengue-consensus-bcftools.fa dengue-direct_10p.bam
  -----------------------------------------------------------------------

>When you run this, you should see an "OK".  Great, so now that we got that out of the way, let\'s build an index for
>our reference genome:

  -----------------------------------------------------------------------
  lofreq faidx dengue-consensus-bcftools.fa
  -----------------------------------------------------------------------

>Now let\'s see what all **lofreq call** can do. Enter the open command
>into the prompt and see what options you have:

  -----------------------------------------------------------------------
  lofreq call
  -----------------------------------------------------------------------

>There\'s a lot here and depending on your particular project or sample
>type, you may want to play with some of these variables. For now, let\'s
>just get on to the fun part and use the standard presets:

  -----------------------------------------------------------------------
  lofreq call -f dengue-consensus-bcftools.fa -o
  dengue-direct-low-frequency-variants.vcf dengue-direct_10p.bam
  -----------------------------------------------------------------------

>This will take a few minutes to run depending on your system, but you
>should soon see a new file in your directory called:
>**dengue-low-frequency-variants.vcf**

Done!

>*9.5 Group practical*

>Now let's take what you've learned and try it out on another dataset! In
>this example, we will be looking for variants from a Chikungunya virus
>alignment.

>To start, navigate to:

-----------------------------------------------------------------------
cd /home/manager/course_data/Consensus_and_variant_calling/09-chikv-lofreq/
-----------------------------------------------------------------------

>You should see four files:

>OMS-chikv-consensus-bcftools.fa
>
>chikv-direct_10p.bam
>
>annotation.txt
>
>chikv-genome.fasta

**Question 1: What mutation is present at position 269?**

**Question 2: What is the allele frequency at position 2318?**
