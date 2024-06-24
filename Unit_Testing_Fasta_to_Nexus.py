import os
import unittest
import Fasta_to_Nexus

class Test_Fasta_To_Nexus(unittest.TestCase):
    def test_parse_fasta(self):
        """
        Test the parse_fasta function to ensure it correctly reads a FASTA file
        and returns a dictionary with sequence names as keys and sequences as values.

        Steps:
        1. Create a string with two sequences.
        2. Write the test string to a fasta file named 'test.fasta'.
        3. Define the expected output as a dictionary with sequence names as keys and sequences as values.
        4. Assert that the output of parse_fasta matches the expected dictionary.
        5. Remove the test fasta file from the system to clean up.
        """
        test_fasta = ">seq1\nATGC\n>seq2\nA--C\n"
        with open("test.fasta", "w") as f:
            f.write(test_fasta)
        expected_output = {"seq1": "ATGC", "seq2": "A--C"}
        self.assertEqual(Fasta_to_Nexus.parse_fasta("test.fasta"), expected_output)
        os.remove("test.fasta")

    def test_nexus_data_header(self):
        """
        Test the nexus_data_header function to ensure it correctly generates
        the NEXUS data header for a given set of sequences.

        Steps:
        1. Define a dictionary of sequences.
        2. Define the expected output as the NEXUS data header string.
        3. Assert that the output of nexus_data_header matches the expected string.
        """
        sequences = {"seq1": "ATGC", "seq2": "A--C"}
        expected_output = (
            "#NEXUS\n\n"
            "BEGIN DATA;\n"
            "DIMENSIONS NTAX=2 NCHAR=4;\n"
            "FORMAT DATATYPE=DNA MISSING=N GAP=-;\n"
            "MATRIX\n"
        )
        self.assertEqual(Fasta_to_Nexus.nexus_data_header(sequences), expected_output)

    def test_nexus_matrix(self):
        """
        Test the nexus_matrix function to ensure it correctly generates
        the NEXUS matrix block for a given set of sequences.

        Steps:
        1. Define a dictionary of sequences.
        2. Define the expected output as the NEXUS matrix block string.
        3. Assert that the output of nexus_matrix matches the expected string.
        """
        sequences = {"seq1": "ATGC", "seq2": "A--C"}
        expected_output = (
            "seq1\tATGC\n"
            "seq2\tA--C\n"
            ";\nEND;\n"
        )
        self.assertEqual(Fasta_to_Nexus.nexus_matrix(sequences), expected_output)

    def test_mrbayes_block(self):
        """
        Test the mrbayes_block function to ensure it correctly generates
        the MrBayes block for given parameters.

        Steps:
        1. Define the number of generations and the outgroup.
        2. Define the name of the test file.
        3. Define the expected output as the MrBayes block string.
        4. Assert that the output of mrbayes_block matches the expected string.
        """
        ngen = 5000
        outgroup = "Placeholder Outgroup"
        test_file = "Test.fasta"

        expected_output = (
            "\nbegin mrbayes;\n"
            "  set autoclose=yes;\n"
            f"  outgroup {outgroup};\n"
            f"  mcmcp ngen={ngen} printfreq=1000 samplefreq=100 diagnfreq=1000 nchains=4 "
            f"savebrlens=yes filename=Test;\n"
            "  mcmc;\n"
            f"  sumt filename=Test;\n"
            "end;\n"
        )
        self.assertEqual(Fasta_to_Nexus.mrbayes_block(test_file, ngen, outgroup), expected_output)

if __name__ == "__main__":
    unittest.main()
