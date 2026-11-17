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
plink --bfile cattle --cow --geno 0.05 --make-bed --out cattle_clean
#      \___ read ___/ \_______ do ________/ \___ write ___/
