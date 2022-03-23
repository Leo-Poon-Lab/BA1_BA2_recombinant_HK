import sys
import gzip
import re
from functools import partial
from Bio import SeqIO

# modified from https://biopython.org/wiki/Sequence_Cleaner

def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'

def sequence_cleaner(fasta_file, por_n, mut_site_filepath):
    # Create our hash table to add the sequences
    _open = partial(gzip.open, mode='rt') if is_gz_file(fasta_file) else open
    # Using the Biopython fasta parse we can read our fasta input

    with open(mut_site_filepath, "r") as mut_site_file:
        mut_sites = mut_site_file.read().split("\n")
        mut_sites = list(filter(None, mut_sites))
        mut_sites = list(map(int, mut_sites))

    with _open(fasta_file) as handle:
        for seq_record in SeqIO.parse(handle, "fasta"):
            # print(seq_record.id)
            sequence = str(seq_record.seq)
            
            # Check if the current sequence is according to the user parameters
            if (
                # float(sequence.count("N"))==0 # only keep sequences with no "N"s
                float(sequence.count("N")) / float(len(sequence)) * 100 <= por_n
            ):
                pos_n = [m.end(0) for m in re.finditer("N", sequence)]
                check = [pos_n_i in mut_sites for pos_n_i in pos_n]
                if not any(check): # no N on mutsites of the interest
                    sys.stdout.write(">" + seq_record.id + "\n" + sequence + "\n")


userParameters = sys.argv[1:]

try:
    sequence_cleaner(userParameters[0], float(userParameters[1]), userParameters[2])
except:
    print("There is a problem while running!")