#!/usr/bin/env python3
"""Generate README.md for a species. Usage: <species> <latin> <ru> <accession> <pubmed>"""
import sys, os, subprocess
from collections import defaultdict

species, latin, ru, acc, pubmed = sys.argv[1:6]

HOME = os.path.expanduser("~")
REPO = f"{HOME}/genome_project/team_repo"
DATA = f"{REPO}/{species}/data"
DA = f"{HOME}/genome_project/downstream_analysis"

# Genome stats
sizes = []
with open(f"{DATA}/{species}.chrom.sizes") as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            sizes.append(int(parts[1]))
total_len = sum(sizes)
n_scaf = len(sizes)
cum, n50 = 0, 0
for s in sorted(sizes, reverse=True):
    cum += s
    if cum >= total_len/2:
        n50 = s
        break

def gc_from_fasta(path):
    gc, tot = 0, 0
    with open(path) as f:
        for line in f:
            if line.startswith('>'): continue
            s = line.strip().upper()
            tot += len(s) - s.count("N")
            gc += s.count('G') + s.count('C')
    return f"{100.0*gc/tot:.2f}" if tot else "N/A"

fna = f"{DA}/{species}_clean.fna"
if not os.path.exists(fna):
    fna = f"{DA}/{species}_filtered.fna"
gc = gc_from_fasta(fna) if os.path.exists(fna) else "N/A"

def linecount(p):
    with open(p) as f: return sum(1 for _ in f)

n_genes = linecount(f"{DATA}/{species}_genes.bed")

tmp_g4 = f"/tmp/{species}_g4_merged.bed"
subprocess.run(f"sort -k1,1 -k2,2n {DATA}/{species}_g4.bed | bedtools merge -i - > {tmp_g4}",
               shell=True, check=True)

g4_tot = linecount(tmp_g4)
zh_tot = linecount(f"{DATA}/{species}_zhunt_merged.bed")
zd_tot = linecount(f"{DATA}/{species}_zdnabert_merged.bed")

def ix(a, b):
    r = subprocess.run(f"bedtools intersect -u -a {a} -b {b} | wc -l",
                       shell=True, capture_output=True, text=True)
    return int(r.stdout.strip())

def pct(n, t):
    return f"{100.0*n/t:.2f}" if t else "0.00"

regions = [("exons","Exons"),("introns","Introns"),
           ("promoters","Promoters (1000 bp up)"),
           ("downstream","Downstream (200 bp)"),
           ("intergenic","Intergenic")]

S = {"g4": tmp_g4,
     "zh": f"{DATA}/{species}_zhunt_merged.bed",
     "zd": f"{DATA}/{species}_zdnabert_merged.bed"}

region_len = {}
for r, _ in regions:
    region_len[r] = 0
    with open(f"{DATA}/{species}_{r}.bed") as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) >= 3:
                try:
                    region_len[r] += int(cols[2]) - int(cols[1])
                except: pass

C, RT, RC = {}, {}, {}
for r, _ in regions:
    rbed = f"{DATA}/{species}_{r}.bed"
    RT[r] = linecount(rbed)
    for k, p in S.items():
        C[(k, r)] = ix(p, rbed)
        RC[(k, r)] = ix(rbed, p)

target_pfams = []
prof = f"{REPO}/profiles/target_profiles.txt"
if os.path.exists(prof):
    with open(prof) as f:
        target_pfams = [ln.strip() for ln in f if ln.strip()]

family_hits = defaultdict(list)
family_name = {}
total_hits = 0
with open(f"{DATA}/{species}_epigenetics.txt") as f:
    for line in f:
        if line.startswith('#'): continue
        cols = line.split()
        if len(cols) < 6: continue
        try:
            evalue = float(cols[4])
        except: continue
        family_hits[cols[3]].append((cols[0], evalue))
        family_name[cols[3]] = cols[2]
        total_hits += 1

family_summary = {}
for pfam, hits in family_hits.items():
    best = {}
    for gid, ev in hits:
        if gid not in best or ev < best[gid]:
            best[gid] = ev
    family_summary[pfam] = sorted(best.items(), key=lambda x: x[1])

m = []
m.append(f"# {latin} ({ru})\n")
m.append(f"NCBI accession: [{acc}](https://www.ncbi.nlm.nih.gov/datasets/genome/{acc}/)\n")
m.append("## Описание организма\nTODO: среда обитания, температура, основные признаки\n")
m.append("## Параметры сборки\n")
m.append("| Параметр | Значение |\n|---|---|")
m.append(f"| Длина генома (scaffolds >1 Mb) | {total_len:,} bp |")
m.append(f"| Число скаффолдов | {n_scaf} |")
m.append(f"| N50 (scaffold) | {n50:,} bp |")
m.append(f"| GC% | {gc} |")
m.append(f"| Число генов | {n_genes:,} |")
m.append(f"| Статей в PubMed | {pubmed} |\n")

m.append("## Таблица 1. Распределение вторичных структур по регионам\n")
m.append("| Участок | G4 (n) | G4 (%) | Zhunt (n) | Zhunt (%) | ZDNABERT (n) | ZDNABERT (%) |")
m.append("|---|---|---|---|---|---|---|")
for r, lb in regions:
    m.append(f"| {lb} | {C[('g4',r)]} | {pct(C[('g4',r)],g4_tot)} | {C[('zh',r)]} | {pct(C[('zh',r)],zh_tot)} | {C[('zd',r)]} | {pct(C[('zd',r)],zd_tot)} |")
m.append(f"| **Всего** | **{g4_tot}** | | **{zh_tot}** | | **{zd_tot}** | |\n")

m.append("## Таблица 2. Доля регионов с хотя бы одной структурой\n")
m.append("| Участок | Всего регионов | % с G4 | % с Zhunt | % с ZDNABERT |")
m.append("|---|---|---|---|---|")
for r, lb in regions:
    m.append(f"| {lb} | {RT[r]} | {pct(RC[('g4',r)],RT[r])} | {pct(RC[('zh',r)],RT[r])} | {pct(RC[('zd',r)],RT[r])} |")
m.append("")

m.append("## Сравнение с фоном\n")
m.append("Фон = доля длины региона в геноме. Обогащение = (% структур в регионе) / (% длины региона в геноме).\n")
m.append("| Участок | Длина (фон, %) | G4 (%) | Zhunt (%) | ZDNABERT (%) | G4 enrich | Zhunt enrich | ZDNABERT enrich |")
m.append("|---|---|---|---|---|---|---|---|")
for r, lb in regions:
    bg = 100.0*region_len[r]/total_len if total_len else 0
    g4p = 100.0*C[('g4',r)]/g4_tot if g4_tot else 0
    zhp = 100.0*C[('zh',r)]/zh_tot if zh_tot else 0
    zdp = 100.0*C[('zd',r)]/zd_tot if zd_tot else 0
    e1 = g4p/bg if bg else 0
    e2 = zhp/bg if bg else 0
    e3 = zdp/bg if bg else 0
    m.append(f"| {lb} | {bg:.2f} | {g4p:.2f} | {zhp:.2f} | {zdp:.2f} | {e1:.2f}× | {e2:.2f}× | {e3:.2f}× |")
m.append("\nОбогащение >1 = регион обогащён структурой относительно фона; <1 = обеднён.\n")

m.append("## Эпигенетические гены (HMMER + Pfam)\n")
m.append(f"Проверено {len(target_pfams)} семейств. Всего хитов: {total_hits}.\n")
m.append("| Семейство | Pfam ID | Уник. генов | Лучший ген (E-value) |")
m.append("|---|---|---|---|")
for pfam in target_pfams:
    if pfam in family_summary:
        gns = family_summary[pfam]
        fn = family_name.get(pfam, "?")
        bg_, be_ = gns[0]
        m.append(f"| {fn} | {pfam} | {len(gns)} | {bg_} ({be_:.1e}) |")
    else:
        m.append(f"| — | {pfam} | 0 | не найдено |")
m.append(f"\nПолные результаты — `data/{species}_epigenetics.txt`.\n")

m.append("## Файлы")
m.append(f"- `data/{species}_g4.bed` — G-квадруплексы")
m.append(f"- `data/{species}_zhunt.bed.gz` — Zhunt z-score >400 (сжато)")
m.append(f"- `data/{species}_zhunt_merged.bed` — Zhunt после merge")
m.append(f"- `data/{species}_zdnabert.bed` — ZDNABERT (threshold 0.5)")
m.append(f"- `data/{species}_zdnabert_merged.bed` — ZDNABERT после merge")
m.append(f"- `data/{species}_epigenetics.txt` — hmmsearch")
m.append(f"- `data/{species}_proteins.fasta.gz` — протеом")
m.append(f"- `data/{species}_<region>.bed` — gene/exon/intron/promoter/downstream/intergenic\n")
m.append("## Использованные команды\nСм. `/scripts/` в корне репозитория.")

with open(f"{REPO}/{species}/README.md", "w") as f:
    f.write("\n".join(m))
os.remove(tmp_g4)
print(f"Создан {REPO}/{species}/README.md")
