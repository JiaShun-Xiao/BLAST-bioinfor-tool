__author__ = 'Jiashun'
import re
import time

starts = time.clock()
hg19 = open("/2_disk/xiaojs/python/blast/hg19.fa")

i = 0
for line in hg19:
    if re.match(r">(.*)", line):
        if i > 0:
            fin = open('/2_disk/xiaojs/python/blast/chrom/chromosome'+str(i)+'.txt', 'w')
            chrom = re.sub(r'\n', '', chrom)
            fin.write(chrom)
            fin.close()
        i += 1
        chrom = ''
    else:
        chrom += line

fin = open('/2_disk/xiaojs/python/blast/chrom/chromosome'+str(i)+'.txt', 'w')
chrom = re.sub(r'\n', '', chrom)
fin.write(chrom)
fin.close()

hg19.close()

end = time.clock()
print "using time: %fs" % (end - starts)