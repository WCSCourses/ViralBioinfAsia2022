# Course software list:
## Software is due 30 June 2022
Please list the software you require on the course virtual machine. 

| Software  Name | Version | Link | Session |
|----------------|---------|------|----------|
| wget | 1.2.1 | https://www.gnu.org/software/wget/manual/wget.html | intro to unix | 
| artic-ncov19 (conda env) | latest | https://github.com/artic-network/artic-ncov2019.git | sars-cov-2 |
| pangolin (conda env) | latest | https://github.com/cov-lineages/pangolin | sars-cov-2 |
| civet (conda env) | latest | https://cov-lineages.org/resources/civet/updating.html | sars-cov-2 |
| spear (conda env) | latest | https://github.com/m-crown/SPEAR | sars-cov-2 |
| usher (conda-env) | latest | https://usher-wiki.readthedocs.io/en/latest/QuickStart.html#quick-install | sars-cov-2 |
| ncov2019-artic-nf (nextflow) | latest | https://github.com/connor-lab/ncov2019-artic-nf | sars-cov-2 |
| trim_galore | latest (0.6.7) | https://github.com/FelixKrueger/TrimGalore | sars-cov-2 |
| minimap2 | latest (2.2.4) | https://github.com/lh3/minimap2 | sars-cov-2, BASH_scripting |
| samtools | latest | http://www.htslib.org/download/ | sars-cov-2, BASH_scripting Coverage_plots NGS_file_formats|
| ivar | latest (1.3.1) | https://github.com/andersen-lab/ivar | sars-cov-2 |
| lofreq | lofreq3 | https://github.com/andreas-wilm/lofreq3 | sars-cov-2 |
| snpeff | latest | http://pcingola.github.io/SnpEff/ | sars-cov-2 |
| mafft | latest | https://mafft.cbrc.jp/alignment/software/ | sars-cov-2 |
| iqtree2 | latest (2.2.0) | http://www.iqtree.org/#download | sars-cov-2 |
| R | latest | https://cran.r-project.org | sars-cov-2 |
| python3 | latest | https://www.python.org/downloads/ | sars-cov-2 |
| FigTree | latest | https://github.com/rambaut/figtree/releases | sars-cov-2 |
| seqtk | latest | https://github.com/lh3/seqtk | sars-cov-2 |
| bbmap | latest | https://sourceforge.net/projects/bbmap/ | sars-cov-2 NGS_file_formats |
| weeSAM | latest | https://github.com/centre-for-virus-research/weeSAM | sars-cov-2 |
| snipit | latest | https://github.com/aineniamh/snipit | sars-cov-2 |
| Trimmomatic-0.39 | 0.39 | https://github.com/usadellab/Trimmomatic/files/5854859/Trimmomatic-0.39.zip | BASH_scripting |
| Qualimap | latest | http://qualimap.conesalab.org/ | Coverage_plots | 
| fastqc| latest | https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ |NGS_file_formats|
|tablet| latest | https://ics.hutton.ac.uk/tablet/download-tablet/| NGS_file_formats Coverage_plots |
| IGV | latest | https://ics.hutton.ac.uk/tablet/download-tablet/ | NGS_file_formats Coverage_plots | 
| SRA toolkit | latest | https://github.com/ncbi/sra-tools/wiki/ and https://github.com/ncbi/sra-tools/wiki/03.-Quick-Toolkit-Configuration | NGS_file_formats | 
| Quast | latest | http://quast.sourceforge.net/quast.html | Coverage_plots | 
| example | example | example | example | 
| **PREVIOUS** | **VM** | **TOOLS** | **BELOW** | 
| Aliview |||| 
| Tempest||||
| Figtree |||| 
| Beast | v1.10.4 |||
| fastqc| latest |||
| spades||||
| idba_ud||||
| abyss||||
| soap-denovo||||
| garm||||
| blastn||||
| quast|||
| scaffold-builder||||
| mummer||||
| diamond||||
| Krona||||
| centrifuge||https://ccb.jhu.edu/software/centrifuge/||
|sra-toolkit ||||
| trim_galore||||
|prinseq-lite||||
|bwa||||
|bowtie2||||
|minimap2||||
|samtools||||
|bcftools||||
|htslib||||
|WeeSAM||||
|tablet|| https://ics.hutton.ac.uk/tablet/download-tablet/||
|ivar||||
|lofreq||||
|SnpEff||||
|VPhaser2||||
|VarScan2||||
|Java JDK||||
|R||||
| example | example | example | example | 
| example | example | example | example | 




 NB: SARS-CoV-2 session - spear will ask if you want to install a pangolin env - options would be to selecte Yes (and not install pangolin separately) - or install pangolin first and then select No when installing spear. Also for ncov2019-artic-nf nextflow - I think the easiest will be to install via the conda route "-profile conda"
 
Diamond and Centrifuge - diamond db of viral proteins and the centrifuge db of bacteria, virus, archea and human.  

krona and sra-toolkit - taxonomy and accession scripts for krona and sra-toolkit.

