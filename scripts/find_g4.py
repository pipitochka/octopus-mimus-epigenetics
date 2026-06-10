import re
import sys

def find_g4(fasta_file, bed_file):
    pattern = re.compile(r'([gG]{3,5}[ATGCatgc]{1,7}){3,}[gG]{3,5}')
    
    with open(fasta_file, 'r') as f_in, open(bed_file, 'w') as f_out:
        chrom = ""
        seq = []
        
        for line in f_in:
            line = line.strip()
            if line.startswith(">"):
                # Если уже накопили последовательность - ищем паттерн
                if chrom:
                    full_seq = "".join(seq)
                    for match in pattern.finditer(full_seq):
                        f_out.write(f"{chrom}\t{match.start()}\t{match.end()}\tG4\n")
                
                chrom = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
                
        # Не забываем последний скаффолд
        if chrom:
            full_seq = "".join(seq)
            for match in pattern.finditer(full_seq):
                f_out.write(f"{chrom}\t{match.start()}\t{match.end()}\tG4\n")

if __name__ == "__main__":
    find_g4(sys.argv[1], sys.argv[2])
