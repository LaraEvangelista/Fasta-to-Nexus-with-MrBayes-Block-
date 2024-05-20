## Fasta to Nexus with MrBayes Block -
This repository contains a Python script that converts a FASTA file to a NEXUS file format with a MrBayes block for phylogenetic analysis.
Project developed in 'Análise de Sequências Biológicas' (Analysis of Biological Sequences), in context of pursuing a Bioinformatics Degree. I was motived to create this simple script after being assigned a homework in said class to write it, and also due to my neccessity of a reliable way to convert fasta files to nexus without using a lot of outside software or bloating my code.

## Files
- `Fasta_to_Nexus.py`: Main script to convert FASTA to NEXUS, appending a MrBayes block and printing it to the console. 
- `Unit_Testing_Fasta_to_Nexus.py`: Unit tests for the functions in `Fasta_to_Nexus.py`.

## Usage
The Fasta_to_Nexus script takes three arguments, the name of the fasta file to convert, the number of generations for the MrBayes block, and the outgroup also for the MrBayes block. Only the fasta file is requred, the ngen has the default value of 5000, and the outgroup defaults to "Placeholder Outgroup". For Linux, an example of how to run the script is:

python3 Fasta_to_Nexus.py file.fasta 4000 Outgroup

The Unit_Testing_Fasta_to_Nexus script takes no arguments. It's just for the purpose of testing the code. For Linux, an example of how to run the script is:

python3 Unit_Testing_Fasta_to_Nexus.py
