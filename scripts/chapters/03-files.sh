#!/usr/bin/env bash
# Commands from the chapter: Genotype files: what PLINK reads and writes
# Run from a folder that holds the input files named in the commands.

plink --file cattle --cow --make-bed --out cattle
