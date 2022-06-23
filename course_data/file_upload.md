Modular directory for data

Check the corresponding data in the course materials directory- course_data
This will be the course_data folder present in the home directory of the VM
Please try to remove files that are not relevant to your module from the respective course data module folder (some may remain from previous iterations of this course) 
Files under 25MB can be uploaded via web server, files up to 100MB can be uploaded via git command line, and if files are bigger than this, you will need to use the Git large file storage system https://git-lfs.github.com (we have an account for up to 50GB for all files on our repo)
Install git lfs 
Clone the repository
Copy your large files to the respective course_data directory
Track the suffix of large files with the following command :
git lfs track "*.sam"
Then add the .gitattributes, commit and push 
git add .gitattributes
git commit -m "now tracking sam files"
git push
Now add your files then commit and push
git add course_data/dir/file.sam
git commit -m "Adding big sam file"
git push
NB - git LFS does NOT track changes in binary files - use different names for modded files

