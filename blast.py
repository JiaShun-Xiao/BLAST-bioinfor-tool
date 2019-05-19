__author__ = 'Jiashun'
import re
import numpy as np
from collections import Counter
from math import ceil
from math import floor

# compare single base
def SingleBaseCompare(seq1,seq2,i,j):
    if seq1[i] == seq2[j]:
        return 2
    else:
        return -1
    
# Smith–Waterman Alignment 
def SMalignment(seq1, seq2):
    m = len(seq1)
    n = len(seq2)
    g = -3
    matrix = []
    for i in range(0, m):
        tmp = []
        for j in range(0, n):
            tmp.append(0)
        matrix.append(tmp)
    for sii in range(0, m):
        matrix[sii][0] = sii*g
    for sjj in range(0, n):
        matrix[0][sjj] = sjj*g
    for siii in range(1, m):
        for sjjj in range(1, n):
            matrix[siii][sjjj] = max(matrix[siii-1][sjjj] + g, matrix[siii - 1][sjjj - 1] + SingleBaseCompare(seq1,seq2,siii, sjjj), matrix[siii][sjjj-1] + g)
    sequ1 = [seq1[m-1]]
    sequ2 = [seq2[n-1]]
    while m > 1 and n > 1:
        if max(matrix[m-1][n-2], matrix[m-2][n-2], matrix[m-2][n-1]) == matrix[m-2][n-2]:
            m -= 1
            n -= 1
            sequ1.append(seq1[m-1])
            sequ2.append(seq2[n-1])
        elif max(matrix[m-1][n-2], matrix[m-2][n-2], matrix[m-2][n-1]) == matrix[m-1][n-2]:
            n -= 1
            sequ1.append('-')
            sequ2.append(seq2[n-1])
        else:
            m -= 1
            sequ1.append(seq1[m-1])
            sequ2.append('-')
    sequ1.reverse()
    sequ2.reverse()
    align_seq1 = ''.join(sequ1)
    align_seq2 = ''.join(sequ2)
    align_score = 0.
    for k in range(0, len(align_seq1)):
        if align_seq1[k] == align_seq2[k]:
            align_score += 1
    align_score = float(align_score)/len(align_seq1)
    return align_seq1, align_seq2, align_score

# Display BlAST result
def Display(seque1, seque2):
    le = 60
    while len(seque1)-le >= 0:
        print('sequence1: ',end='')
        for a in list(seque1)[le-40:le]:
            print(a,end='')
        print("\n")
        print('           ',end='')
        for k in range(le-40, le):
            if seque1[k] == seque2[k]:
                print('|',end='')
            else:
                print(' ',end='')
        print("\n")
        print('sequence2: ',end='')
        for b in list(seque2)[le-40:le]:
            print(b,end='')
        print("\n")
        le += 40
    if len(seque1) > le-40:
        print('sequence1: ',end='')
        for a in list(seque1)[le-40:len(seque1)]:
            print(a,end='')
        print("\n")
        print('           ',end='')
        for k in range(le-40, len(seque1)):
            if seque1[k] == seque2[k]:
                print('|',end='')
            else:
                print(' ',end='')
        print("\n")
        print('sequence2: ',end='')
        for b in list(seque2)[le-40:len(seque2)]:
            print(b,end='')
        print("\n")

# transform base to numeric value
def WordToNum(word):
    tmp = []
    trans = {'A':1,'C':2,'G':3,'T':4}
    for w in word:
        tmp.append(trans[w])
    return tmp

# transform word with 11 bases to its index
def WordToIndex(word,word_len):
    tmp = 0
    word_num = WordToNum(word)
    for i,v in enumerate(word_num):
        tmp += (v-1)*4**(word_len-i)
    return tmp   

# Get word's postion in genome from library
def GetWordPos(word):
    assert len(word)== 11
    seek_index = WordToIndex(word,11-1)
    positions = []
    for chr_name in chr_names:
        chr_seq = open('/home/jxiaoae/class/blast/chromosome_{}_library.txt'.format(chr_name),'r')
        seeks = np.load("chromosome_{}_library_seeks.npy".format(chr_name))
        chr_seq.seek(seeks[seek_index,0])
        position = chr_seq.read(seeks[seek_index,1])
        try:
            positions.append(list(map(int, position[:-1].split(","))))
        except:
            positions.append([])
    return positions

# Extract subsequence from GRCh37 file
def ExtractSeq(chr_index,pos,length):
    pos = pos+floor(pos/60)
    hg19.seek(chrom_seek_index[chr_index,1]+pos-1)
    return re.sub(r'\n', '', hg19.read(length))

# main blast function
def Blast(query_seq):
    i = 0
    query_words = []
    query_seq_length = len(query_seq)
    words_length = query_seq_length-11+1
    while i < words_length:
        query_words.append(query_seq[i:i+11])
        i += 1
    words_positions = []
    for word in query_words:
        words_positions.append(GetWordPos(word))
    for chr_index in range(24):
        for word_index in range(words_length):
            for pos in range(len(words_positions[word_index][chr_index])):
                words_positions[word_index][chr_index][pos] += words_length - word_index - 1
        
        words_positions_corrects = []
        for word_index in range(words_length):
            words_positions_corrects += words_positions[word_index][chr_index]
        
        words_positions_corrects_count = Counter(words_positions_corrects)
        finded_postions = []
        for count_ in words_positions_corrects_count:
            # we can select the bigger threshold of words_positions_corrects_count[count_] just 
            # like we select the highly similar sequence in NCBI BLAST
            if words_positions_corrects_count[count_] > 5:
                finded_postions.append(count_)
        if finded_postions:
            chr_seq = open('/home/jxiaoae/class/blast/chromosome_{}_library.txt'.format(chr_names[chr_index]),'r')
            chr_seq = chr_seq.read().strip()
            for finded_postion in finded_postions:
                candidate_seq_pos = finded_postion - query_seq_length + 11 - 5
                candidate_seq_length = query_seq_length + 11
                candidate_sequence = ExtractSeq(chr_index,candidate_seq_pos,candidate_seq_length)
                i_start_indexs = []
                for i_start in range(15):
                    _,_,score = SMalignment(candidate_sequence[i_start:],query_seq)
                    i_start_indexs.append(score)
                i_start = np.array(i_start_indexs).argmax()
                i_end_indexs = []
                for i_end in range(1,16):
                    _,_,score = SMalignment(candidate_sequence[:-i_end],query_seq)
                    i_end_indexs.append(score)
                i_end = np.array(i_end_indexs).argmax()+1
                candidate_sequence = candidate_sequence[i_start:-i_end]
                align_seq1,align_seq2,align_score = SMalignment(candidate_sequence,query_seq)
                if align_score>0.7:
                    print("find in chromosome "+chr_names[chr_index]+": "+str(candidate_seq_pos+i_start)+' ---> '+str(candidate_seq_pos+i_start+len(candidate_sequence)-1)+", align score: "+str(align_score))
                    Display(align_seq1, align_seq2)
    return None

if __name__ == "__main__":
    chr_names = np.load('/home/jxiaoae/class/blast/GRCh37_chr_names.npy')
    chrom_seek_index = np.load('/home/jxiaoae/class/blast/GRCh37_chrom_seek_index.npy')
    hg19 = open("/home/share/GRCh37/human_g1k_v37.fasta")
    query_sequence = 'GTATCGGAACTTCCAACTTGTAGGCAAAATAGATATGCTTCATATTCTTAAAAACCACAAGAAA'
    Blast(query_sequence)