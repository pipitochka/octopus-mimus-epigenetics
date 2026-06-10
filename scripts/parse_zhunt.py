import sys
import glob
import os

def parse_zhunt_folder(in_dir, out_bed):
    # Ищем все файлы Zhunt
    files = glob.glob(f"{in_dir}/*.Z-SCORE")
    
    with open(out_bed, 'w') as f_out:
        for filepath in files:
            # Вырезаем имя хромосомы (убираем .fna.Z-SCORE)
            chrom = os.path.basename(filepath).replace(".fna.Z-SCORE", "")
            
            with open(filepath, 'r') as f_in:
                lines = f_in.readlines()
                if not lines: continue
                
                # Читаем данные
                for i in range(1, len(lines)):
                    parts = lines[i].strip().split()
                    if len(parts) >= 3:
                        try:
                            # 3й столбец - это Z-score
                            z_score = float(parts[2])
                            
                            if z_score > 400:
                                # Поскольку Zhunt идет окном по порядку, start = i
                                start = i
                                end = i + len(parts[3])
                                f_out.write(f"{chrom}\t{start}\t{end}\tZhunt\n")
                        except:
                            continue

if __name__ == "__main__":
    parse_zhunt_folder("mimus_split", "mimus_zhunt.bed")
    parse_zhunt_folder("rubescens_split", "rubescens_zhunt.bed")
