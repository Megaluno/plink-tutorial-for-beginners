#!/usr/bin/env bash
# Commands from the chapter: Installing PLINK and running your first command
# Run from a folder that holds the input files named in the commands.

cd C:\genomics
plink --version
cd ~/genomics
chmod +x plink
./plink --version
mkdir -p ~/genomics && cd ~/genomics
unzip plink_linux_x86_64_*.zip
chmod +x plink
./plink --version
# optional: make it available everywhere
sudo cp plink /usr/local/bin/
# 1. once only: tell conda where to look (Bioconda needs conda-forge)
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# 2. a separate environment for this course, with both versions
conda create -n plink_course plink plink2
conda activate plink_course

# 3. check
plink --version
plink2 --version
which plink          # e.g. /home/you/miniforge3/envs/plink_course/bin/plink
plink --bfile cattle --cow --geno 0.05 --make-bed --out cattle_clean
#      \___ read ___/ \_______ do ________/ \___ write ___/
