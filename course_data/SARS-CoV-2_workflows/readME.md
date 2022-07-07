## ~/SARS-CoV-2/

* MinION: folder with an example MinION run data - fast5, fastq, seq_sum etc - 3 samples to keep size down ~5GB
* Illumina: folder with a small example Illumina run data - fastq - 3 samples ~0.5GB
* Variants: folder with various sample genome consensus seqs for varaint calling - size is ~10MB
* Phylo: folder with sample genome sequences along with pre-built alignments, trees and metadata size is ~100MB
* Group: folder with various illumina sample data for full analysis ~0.5GB
* Extra: extra data, simulated, size is ~100MB

## Modular directory for data

Upload corresponding data in the course materials directory- course_data under your specific module 

This will be the course_data folder present in the home directory of the VM
### Upload limits
Files under 25MB can be uploaded via web server, files up to 100MB can be uploaded via git command line, and if files are bigger than this, you will need to use the Git large file storage system https://git-lfs.github.com (we have an account for up to 50GB for all files on our repo)
### how to use git lfs
* Install git lfs - instructions at https://git-lfs.github.com
* Clone this repository
* Copy your large files to the respective course_data directory
* Track the suffix of large files with the following command :
```
git lfs track "*.sam"
```
Then add the .gitattributes, commit and push 
```
git add .gitattributes
git commit -m "now tracking sam files"
git push
```
Now add your files then commit and push
```
git add course_data/dir/file.sam
git commit -m "Adding big sam file"
git push
```
NB - git LFS does NOT track changes in binary files - use different names for modded files of this type

