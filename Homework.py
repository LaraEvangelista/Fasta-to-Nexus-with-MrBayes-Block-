import sys
import os
from Bio import SeqIO

def parse_fasta(fasta_file):

    """
    Parses a FASTA file and returns a dictionary with sequence names as the keys and sequences as the values.
    
    Args:
        file (str): Path to the FASTA file.
    
    Returns:
        dict: A dictionary where keys are sequence names and values are sequences.
    """

    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        name = record.id[:99]  
        sequences[name] = str(record.seq)
    return sequences

def nexus_data_header(sequences):
    
    """
    Generates the NEXUS DATA header from the given sequences.
    
    Args:
        sequences (dict): A dictionary where the keys are the sequence names and the values are the sequences.
    
    Returns:
        str: The NEXUS DATA header as a string.
    """

    taxa_count = len(sequences)
    seq_length = len(next(iter(sequences.values())))
    header = (
        "#NEXUS\n\n"
        "begin data;\n"
        f"dimensions ntax={taxa_count} nchar={seq_length};\n"
        "format datatype=dna missing=N gap=-;\n"
        "matrix\n"
    )
    return header

def nexus_matrix(sequences):

    """
    Generates the NEXUS MATRIX block from the given sequences.
    
    Args:
        sequences (dict): A dictionary where the keys are the sequence names and the values are the sequences.
    
    Returns:
        str: The NEXUS MATRIX block as a string.
    """

    matrix = ""
    for name, seq in sequences.items():
        matrix += f"{name} {seq}\n"
    matrix += ";\nend;\n"
    return matrix

def generate_mrbayes_block(fasta_file, ngen, outgroup):
    
    """
    Generates a MrBayes block for the NEXUS file.
    
    Args:
        fasta_file (str): name of the fasta file, 
        ngen (int): Number of generations for MrBayes.
        outgroup (str): Name of the outgroup.
    
    Returns:
        str: The MrBayes block as a string.
    """

    base_filename = os.path.splitext(fasta_file)[0]

    block = (
        "begin mrbayes;\n"
        "  set autoclose=yes;\n"
        f"  outgroup {outgroup};\n"
        f"  mcmcp ngen={ngen} printfreq=1000 samplefreq=100 diagnfreq=1000 nchains=4 "
        f"savebrlens=yes filename={base_filename};\n"
        "  mcmc;\n"
        f"  sumt filename={base_filename};\n"
        "end;\n"
    )
    return block

if __name__ == "__main__":

    print(len(sys.argv))

    if len(sys.argv) < 2 or len(sys.argv) > 3:
        raise Exception("Error, please write the arguments as follows: fasta_filename(required) ngen(optional) outgroup(optional)")

    fasta_file = sys.argv[1]
    ngen = sys.argv[2] if len(sys.argv) < 3 else 5000
    outgroup = sys.argv[3] if len(sys.argv) < 4 else "Placeholder Outgroup"

    sequences = parse_fasta(fasta_file)
    nexus_data_header = nexus_data_header(sequences)
    nexus_matrix = nexus_matrix(sequences)
    mrbayes_block = generate_mrbayes_block(fasta_file, ngen, outgroup)

    nexus_output = nexus_data_header + nexus_matrix + mrbayes_block
    print(nexus_output)
