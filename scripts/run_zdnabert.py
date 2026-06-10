import torch
from transformers import BertTokenizer, BertForTokenClassification
import numpy as np
from Bio import SeqIO
import scipy.ndimage
from tqdm import tqdm
import sys

# Параметры из оригинала
model_confidence_threshold = 0.5
minimum_sequence_length = 10

def seq2kmer(seq, k):
    return [seq[x:x+k] for x in range(len(seq)+1-k)]

def split_seq(seq, length=512, pad=16):
    res = []
    for st in range(0, len(seq), length - pad):
        end = min(st+512, len(seq))
        res.append(seq[st:end])
    return res

def stitch_np_seq(np_seqs, pad=16):
    res = np.array([])
    for seq in np_seqs:
        res = res[:-pad] if len(res) > 0 else res
        res = np.concatenate([res, seq])
    return res

def predict_zdna(fasta_in, bed_out):
    print(f"Загрузка модели для {fasta_in}...")
    tokenizer = BertTokenizer.from_pretrained('ZDNABERT_model/')
    model = BertForTokenClassification.from_pretrained('ZDNABERT_model/')
    
    # Используем CPU - у нас 32 крутых ядра!
    device = torch.device('cpu') 
    model.to(device)
    model.eval()

    print(f"Начинаем анализ: {fasta_in}")
    
    with open(bed_out, 'w') as fh_out:
        for seq_record in SeqIO.parse(fasta_in, 'fasta'):
            chrom = seq_record.name
            print(f"Обработка хромосомы: {chrom}")
            
            seq_str = str(seq_record.seq).upper()
            
            # ЧТОБЫ ИЗБЕЖАТЬ OOM (ошибки памяти), берем геном кусками
            chunk_size = 1000000 
            
            for chunk_start in range(0, len(seq_str), chunk_size):
                chunk_seq = seq_str[chunk_start : chunk_start + chunk_size + 100]
                
                kmer_seq = seq2kmer(chunk_seq, 6)
                if not kmer_seq: continue
                
                seq_pieces = split_seq(kmer_seq)
                
                with torch.no_grad():
                    preds = []
                    for seq_piece in seq_pieces:
                        input_ids = torch.LongTensor(tokenizer.encode(' '.join(seq_piece), add_special_tokens=False)).to(device)
                        outputs = torch.softmax(model(input_ids.unsqueeze(0))[0], axis=-1)[0, :, 1]
                        preds.append(outputs.numpy())
                        
                stitched_preds = stitch_np_seq(preds)
                labeled, max_label = scipy.ndimage.label(stitched_preds > model_confidence_threshold)
                
                for label in range(1, max_label+1):
                    candidate = np.where(labeled == label)[0]
                    if candidate.shape[0] > minimum_sequence_length:
                        abs_start = chunk_start + candidate[0]
                        abs_end = chunk_start + candidate[-1]
                        # Пишем в BED формат!
                        fh_out.write(f"{chrom}\t{abs_start}\t{abs_end}\tZDNABERT\n")

if __name__ == "__main__":
    predict_zdna(sys.argv[1], sys.argv[2])
