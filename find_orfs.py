import os
from Bio import SeqIO
import argparse

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
    for strand, nucleotide in [(+1, record.seq), (-1, record.seq.reverse_complement())]:    # reads in each direction of sequence
        for frame in range(3):  # breaks into 3 frames
            length = 3 * ((len(record)-frame) // 3 )
            for nucleotide_seq in nucleotide[frame:frame+length].split('*'): # reads nucleotide sequece in slices according to orf
                local_orfs = []
                if len(nucleotide_seq) >= min_nucleotide_len:
                    local_orfs.append(record.description) # record description
                    local_orfs.append(nucleotide_seq)   # appends orf nucleotide sequence
                    local_orfs.append(strand)   # appends the strand direction (1 or -1)
                    local_orfs.append(frame)    # appends the frame (0, 1, 2)
                    all_orfs.append(local_orfs) # appends local_orfs of a single record to all_orfs of all records
    return all_orfs

'''
    Input: record_desc, nucleotide, strand, frame
    Output: None
'''
def write_to_fasta_file(record_desc, nucleotide, strand, frame):
    with open(output_filename, 'a') as f:
        f.write(f'{record_desc} - ORF on strand {strand}, frame {frame} - length {len(nucleotide)} \n')
        f.write(f"{nucleotide}\n")

'''
    Validates Input and Output file formats
    Input: None
    Output: input_file_name, output_file_name
'''
def check_args():
    input_file = ""
    output_filename = ""

    # validates input file name
    if args.filename.endswith(".fasta"):
        input_file = args.filename
    else:
        print("Input file is not a FASTA file")
        quit()

    # validates output name
    if args.output_filename is None:
        output_filename = f'orfs_{input_file}'
        print(f'No output file provided, defaulting to {output_filename}')
    elif args.output_filename.endswith(".fasta"):
        output_filename = args.output_filename
    else:
        print("Output file is not a FASTA file")
        quit()
    return input_file, output_filename

parser = argparse.ArgumentParser(description="Finds all possible ORFs in a DNA sequence from a fasta file.")
parser.add_argument("filename", help="Path to fasta file", metavar="<Input FASTA file>", type=str)
parser.add_argument("output_filename", nargs='?', help="Path to fasta file output", 
                    metavar="<Output FASTA file>", type=str, default=None)
args = parser.parse_args()

input_file, output_filename = check_args()

if __name__ == '__main__':
    records = read_fasta_file(input_file)
    all_orfs = []
    for record in records:
        record_orfs = find_orfs(record)
        all_orfs.append(record_orfs)

    for orf in all_orfs:
        for record_desc, nucleotide, strand, frame in orf:
            #print(f'{record_desc} - {nucleotide} - {strand} - {frame} \n')
            write_to_fasta_file(record_desc, nucleotide, strand, frame)