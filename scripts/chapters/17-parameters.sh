#!/usr/bin/env bash
# Commands from the chapter: How each setting changes your results
# Run from a folder that holds the input files named in the commands.

plink --bfile cattle --cow --maf 0.01 --hwe 1e-6 --geno 0.05 --mind 0.05 --make-bed --out order_demo
