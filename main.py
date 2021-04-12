import os
import itertools
from Bio import SeqIO

# TODO:
# make cli interface
# write orfs to a fasta file

'''
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
    Output: returns list of pros, strands, frames of each record
'''
def find_orfs(record):
    table = 11
    min_pro_len = 100
    pros, strands, frames = [], [], []
    for strand, nucleotide in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
        for frame in range(3):
            length = 3 * ((len(record)-frame) // 3 )     # multiples of 3 
            for pro in nucleotide[frame:frame+length].translate(table).split('*'):
                if len(pro) >= min_pro_len:
                    pros.append(pro)
                    strands.append(strand)
                    frames.append(frame)
                    #print(f'{pro[:30]} {pro[-3:]} - length {len(pro)}, strand {strand}, frame {frame}')
    return pros, strands, frames

if __name__ == '__main__':
    fasta_files = get_list_of_files()
    
    for file in fasta_files:
        records = read_fasta_file(file)
        for record in records:
            pros, strands, frames = find_orfs(record)
            for pro, strand, frame in zip(pros, strands, frames):
                print(f'{pro[:30]} {pro[-3:]} - length {len(pro)}, strand {strand}, frame {frame}')
