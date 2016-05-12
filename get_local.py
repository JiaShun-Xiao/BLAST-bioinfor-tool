__author__ = 'Jiashun'
import re
import time
import os
import string
from collections import Counter

starts = time.clock()


def s(si_, sj_):
    if qu[si_] == seq2[sj_]:
        return 2
    else:
        return -1


def alignment(seq1, seq2):
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
            matrix[siii][sjjj] = max(matrix[siii-1][sjjj] + g, matrix[siii - 1][sjjj - 1] + s(siii, sjjj), matrix[siii][sjjj-1] + g)
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
    seque1 = string.join(sequ1, '')
    seque2 = string.join(sequ2, '')
    global score
    score = 0
    for k in range(0, len(seque1)):
        if seque1[k] == seque2[k]:
            score += 1
    score = float(score)/len(seque2)
    return seque1, seque2


def display(seque1, seque2):
    le = 40
    while len(seque1)-le >= 0:
        print 'sequence1: ',
        for a in list(seque1)[le-40:le]:
            print a,
        print "\n"
        print '           ',
        for k in range(le-40, le):
            if seque1[k] == seque2[k]:
                print '|',
            else:
                print ' ',
        print "\n"
        print 'sequence2: ',
        for b in list(seque2)[le-40:le]:
            print b,
        print "\n"
        le += 40
    if len(seque1) > le-40:
        print 'sequence1: ',
        for a in list(seque1)[le-40:len(seque1)]:
            print a,
        print "\n"
        print '           ',
        for k in range(le-40, len(seque1)):
            if seque1[k] == seque2[k]:
                print '|',
            else:
                print ' ',
        print "\n"
        print 'sequence2: ',
        for b in list(seque2)[le-40:len(seque2)]:
            print b,
        print "\n"


def split_w(w):
    loc = dict()
    loc[w] = []
    if os.path.exists('/2_disk/xiaojs/python/blast/library/'+w):
        fin = open('/2_disk/xiaojs/python/blast/library/'+w, 'r')
        for ii in range(1, 26):
            loc[w].append([])
        for line in fin:
            if re.match(r"[\d]:(.*)", line):
                ch = line[0:1]
            elif re.match(r"[\d]{2}:(.*)", line):
                ch = line[0:2]
            else:
                line = map(int, (line.split()))
                loc[w][int(ch)-1] = line
        fin.close()
        return loc[w]
    else:
        for i in range(0, 25):
            loc[w].append([])
        return loc[w]


def chromosome(words):
    loc = dict()
    for w in words:
        loc[w] = split_w(w)
    lenw = len(loc)
    global qu
    lenq = len(qu)
    qu_ = qu[0:6]
    qu__ = qu[-6:len(qu)]
    for j in range(0, 25):
        for ww in range(0, lenw-1):
            for jj in range(0, len(loc[words[ww]][j])):
                loc[words[ww]][j][jj] += lenw - ww - 1
        chj = []
        for www in words:
            chj += loc[www][j]
        chjj = Counter(chj)
        local = []
        for chjj_ in chjj:
            # we can select the bigger threshold of chjj[chjj_] just like we select the highly similar sequence in NCBI BLAST
            if chjj[chjj_] > 5:
                local.append(chjj_)
        if local:
            chromo = open('/2_disk/xiaojs/python/blast/chrom/chromosome'+str(change[j+1])+".txt", 'r')
            chromoso = chromo.read().strip()
            for local_ in local:
                temp = chromoso[local_-lenq:local_+15]
                for iii in range(0, lenq):
                    if temp[iii:iii+5].upper() in qu_:
                        break
                for iiii in range(0, 15):
                    if temp[-5-iiii:len(temp)-iiii].upper() in qu__:
                        break
                global seq2
                seq2 = temp[iii:len(temp)-iiii].upper()
                (seque1, seque2) = alignment(qu, seq2)
                if score > 0.5:
                    print "find in chromosome"+str(j+1)
                    print str(local_-lenq+iii)+' ---> '+str(local_+lenq-iiii)
                    print temp[iii:len(temp)-iiii]+"\n"
                    display(seque1, seque2)
                    print "the score is "+str(score)+"\n"
            chromo.close()
        if not local:
            continue

#  qu is the query sequence
qu = 'CTAAAACCAAGAGCGGGAGGGGACGGGGCTGCCGCAGCCCTCCCAGA'
#qu = 'ACTGCCAGACACATTCATGACCTAATCCCTACATTAGCATAAGAGGGCACACTCTCCTCCTATGGGGGAAACTGAGGTACGAAGAget_local.pyACTAAAGTGACTTTCCCACAGCTGGTGGGAGGCAGACGGGAAATTCACA'

i = 0
word = []
while i <= len(qu)-11:
    word.append(qu[i:i+11])
    i += 1
change = {1:1, 10:2, 11:3, 12:4, 13:5, 14:6, 15:7, 16:8, 17:9, 18:10, 19:11, 2:12, 20:13, 21:14, 22:15, 3:16,
          4:17, 5:18, 6:19, 7:20, 8:21, 9:22, 23:23, 24:24, 25:25 }

seq2 = ''
score = 0
chromosome(word)

end = time.clock()
print "using time: %fs" % (end - starts)
print '----------------------------------------------------------------------------------------------'


