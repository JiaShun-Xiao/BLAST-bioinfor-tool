__author__ = 'Jiashun'
import re
import time

starts = time.clock()

for ch in range(1, 25):
    file1 = open("/2_disk/xiaojs/python/blast/chrom/chromosome"+str(ch)+".txt", 'r')
    #file1 = open("EGFR.fasta", 'r')
    #file2 = open("EGFR.fasta", 'r')
    lib = dict()
    line = file1.read().strip()
    ii = 0
    le = len(line)
    while ii < le-11:
        w = line[ii:ii+11]
        w = w.upper()
        if lib.has_key(w):
            lib[w] += (" "+str(ii))
        else:
            lib[w] = str(ii)
        ii += 1
    for l in lib:
        fin = open('/2_disk/xiaojs/python/blast/library/'+l, 'a')
        #fin = open('D:/python/PythonProject/lib/'+l+'.txt', 'a')
        fin.write(str(ch)+": \n")
        fin.write(lib[l]+"\n")
        fin.close()
    file1.close()

end = time.clock()
print "using time: %fs" % (end - starts)

