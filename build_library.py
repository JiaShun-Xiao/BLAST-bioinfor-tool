from numba import jit
import numpy as np
import re
from multiprocessing import Pool
from math import ceil

@jit
def BaseToNum(chr_seq):
    chr_seq = re.sub(r'A', '1', chr_seq)
    chr_seq = re.sub(r'C', '2', chr_seq)
    chr_seq = re.sub(r'G', '3', chr_seq)
    chr_seq = re.sub(r'T', '4', chr_seq)
    return chr_seq

@jit
def BaseToIndex(word,word_len):
    tmp = 0
    for i,v in enumerate(word):
        tmp += (int(v)-1)*4**(word_len-i)
    return tmp    

@jit
def GenSeek(library,word_len):
    seeks = np.zeros((4**word_len,2),dtype=int)
    tmp = 0
    for i,l in enumerate(library):
        seeks[i,0] = tmp
        seeks[i,1] = len(l)
        tmp += len(l)
    return seeks

def BuildLibrary(chr_name):
    word_len = 11
    chr_seq = chrom_dict[chr_name]
    chr_seq = BaseToNum(chr_seq)
    chr_len = len(chr_seq)
    library = np.zeros(4**word_len,dtype=str).tolist()
    ii = 0
    while ii<chr_len-word_len:
        w = chr_seq[ii:ii+word_len]
        ii += 1
        if 'N' in w:
            continue
        try:
            library[BaseToIndex(w,word_len-1)] += str(ii)+","
        except:
            pass
    
    seeks = GenSeek(library,word_len)
    lib_seq = ''.join(library)
    with open('/home/jxiaoae/class/blast/chromosome_{}_library.txt'.format(chr_name), 'w') as f:
        f.write(lib_seq)
        f.close()
    np.save('/home/jxiaoae/class/blast/chromosome_{}_library_seeks.npy'.format(chr_name),seeks)
    

if __name__ == '__main__':
    hg19 = open("/home/share/GRCh37/human_g1k_v37.fasta")
    head = True
    chrom_dict = {}
    head_line = []
    chr_names = []
    for line in hg19:
        if re.match(r">[1-9X-Y]|[12][0-9]",line):
            head_line.append(line)
            if head:
                head = False
            else:
                chr_seq = re.sub(r'\n', '', chr_seq)
                chr_seq = chr_seq.upper()
                chrom_dict[chr_name] = chr_seq
            chr_name = line.split()[0][1:]
            chr_names.append(chr_name)
            chr_seq = ''
            print(chr_name,end=",")
        else:
            chr_seq += line
    chr_seq = re.sub(r'\n', '', chr_seq)
    chr_seq = chr_seq.upper()
    chrom_dict[chr_name] = chr_seq
    chrom_seek_index = np.array([[int(line.split(":")[-2]),len(line)] for line in head_line])
    for i in range(1,24):
        chrom_seek_index[i,1]=chrom_seek_index[i,1]+chrom_seek_index[i-1,1]+chrom_seek_index[i-1,0]+ceil(chrom_seek_index[i-1,0]/60)
    np.save('/home/jxiaoae/class/blast/GRCh37_chrom_seek_index.npy',chrom_seek_index)
    np.save('/home/jxiaoae/class/blast/GRCh37_chr_names.npy',np.array(chr_names))
    print(chr_names)
    # reset multiprocessing num according to your server
    with Pool(10) as p:
        p.map(BuildLibrary, chr_names)
