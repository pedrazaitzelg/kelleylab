from Bio import SeqIO
import sys
import os
import shutil
import csv

### STDIN: list of genome names
### run in a directory containing 

for f in os.listdir():
    genome = f[:-4]      ## remove .faa from the filename
    print(genome)
    old_faa = genome + ".faa"
    new_faa = genome + ".prepped.faa"
    record_list = list(SeqIO.parse(old_faa, 'fasta'))

    with open(new_faa, 'w') as corrected:
        for i in range(len(record_list)):
            record = record_list[i]
            record.id = record.id + '_' + genome ## rename records to have genome name in them
            record.description = ''
            if '*' in record.seq:
                if record.seq[-1] == '*': ## remove stop codon from end of sequences
                    record.seq = record.seq[:-1]
                    SeqIO.write(record, corrected, 'fasta')
            else:
                SeqIO.write(record, corrected, 'fasta')
