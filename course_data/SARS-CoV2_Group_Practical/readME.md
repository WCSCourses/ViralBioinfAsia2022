## Delete

The data for this is being included in SARS-CoV2 workflows - in a group subfolder there

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

