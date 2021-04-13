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
    nucleotides, strands, frames = [], [], []
    for strand, nucleotide in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
        for frame in range(3):
            length = 3 * ((len(record)-frame) // 3 )     # multiples of 3 
            for nucleotide in nucleotide[frame:frame+length].split('*'):
                if len(nucleotide) >= min_nucleotide_len:
                    nucleotides.append(nucleotide)
                    strands.append(strand)
                    frames.append(frame)
                    #print(f'{pro[:30]} {pro[-3:]} - length {len(pro)}, strand {strand}, frame {frame}')
    return nucleotides, strands, frames

'''
    Input: record_desc, nucleotide, strand, frame
    Output: None
'''
def write_to_fasta_file(record_desc, nucleotide, strand, frame):
    with open('orfs_'+file, 'a') as f:
        f.write(f'{record_desc} - ORF on strand {strand}, frame {frame} - length {len(nucleotide)} \n')
        f.write(f"{nucleotide}\n")

parser = argparse.ArgumentParser(description="Finds all possible ORFs in a DNA sequence from a fasta file.")
parser.add_argument("filename", help="Path to fasta file", metavar="<FASTA file>", type=str)

args = parser.parse_args()
file = args.filename

if __name__ == '__main__':
    records = read_fasta_file(file)
    for record in records:
        record_description = record.description
        nucleotides, strands, frames = find_orfs(record)
        for nucleotide, strand, frame in zip(nucleotides, strands, frames):
            write_to_fasta_file(record_description, nucleotide, strand, frame)
            #print(f'{pro[:30]} {pro[-3:]} - length {len(pro)}, strand {strand}, frame {frame}')
