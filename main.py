import os
import itertools
from Bio import SeqIO
import argparse

# TODO:
# make cli interface
# write orfs to a fasta file -- this is kinda done, it needs more work

'''
    THIS IS DEPRECATED
    Input: None
    Output: returns list of files containing '*.fasta'
'''
def get_list_of_files():
    files = []
    for filename in os.listdir():
        if filename.endswith('.fasta'):
            files.append(filename)
    return files

'''
    Input: file name
    Output: returns list of record objects
'''
def read_fasta_file(file):
    records = []
    for record in SeqIO.parse(file, 'fasta'):
        records.append(record)
    return records

'''
    Input: a single record
    Output: returns list of nucleotides, strands, frames of each record
'''
def find_orfs(record):
    min_nucleotide_len = 100
    all_orfs = []
    for strand, nucleotide in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
        for frame in range(3):
            length = 3 * ((len(record)-frame) // 3 )     # multiples of 3 
            for nucleotide in nucleotide[frame:frame+length].split('*'):
                local_orfs = []
                if len(nucleotide) >= min_nucleotide_len:
                    local_orfs.append(record.description)
                    local_orfs.append(nucleotide)
                    local_orfs.append(strand)
                    local_orfs.append(frame)
                    all_orfs.append(local_orfs)
    return all_orfs

'''
    Input: record_desc, nucleotide, strand, frame
    Output: None
'''
def write_to_fasta_file(record_desc, nucleotide, strand, frame):
    with open(output_filename, 'a') as f:
        f.write(f'{record_desc} - ORF on strand {strand}, frame {frame} - length {len(nucleotide)} \n')
        f.write(f"{nucleotide}\n")

parser = argparse.ArgumentParser(description="Finds all possible ORFs in a DNA sequence from a fasta file.")
parser.add_argument("filename", help="Path to fasta file", metavar="<Input FASTA file>", type=str)
parser.add_argument("output_filename", nargs='?', help="Path to fasta file output", 
                    metavar="<Output FASTA file>", type=str, default=None)

args = parser.parse_args()
file = args.filename
global output_filename
if args.output_filename is None:
    output_filename = f'orfs_{file}'
    print(f'No output file provided, defaulting to {output_filename}')
else:
    output_filename = args.output_filename

if __name__ == '__main__':
    records = read_fasta_file(file)
    all_orfs = []
    for record in records:
        record_orfs = find_orfs(record)
        all_orfs.append(record_orfs)

    for orf in all_orfs:
        for record_desc, nucleotide, strand, frame in orf:
            #print(f'{record_desc} - {nucleotide} - {strand} - {frame} \n')
            write_to_fasta_file(record_desc, nucleotide, strand, frame)