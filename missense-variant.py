#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import re
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser
from collections import defaultdict
from matplotlib import cm

script_path = os.path.split(os.path.abspath(sys.argv[0]))[0]
usage = script_path+'''/%prog [options]

------------------------------------------------------------
Get the missense variant type(HOM or HET) and draw a picture
------------------------------------------------------------
'''
parser = OptionParser(usage=usage)
parser.add_option("-i", "--input", dest="input_dir", help="the input directory, which only \
contains the information files of the samples")
(options, args) = parser.parse_args()

if not options.input_dir:
    parser.print_help()
    sys.exit("\nError: the input directory is required...")
try:
    os.mkdir("missense_variant")
except OSError:
    pass
inputdir = os.path.abspath(os.path.expanduser(options.input_dir))
samples = os.listdir(inputdir)
os.chdir("missense_variant")
samplelist = []
for name in samples:
    sample = re.findall(r"(^.*?)\.", name)[0]
    samplelist.append(sample)
samplelist = list(set(samplelist))
samplelist.sort()

with open("samples.list", "w") as file1:
    for i in samplelist:
        print >> file1, i,
        for j in samples:
            if re.search(r"^{0}\..*\.vcf".format(i), j):
                vcf = "{0}/{1}".format(inputdir, j)
            if re.search(r"^{0}\..*\.txt".format(i), j):
                txt = "{0}/{1}".format(inputdir, j)
        print >> file1, vcf, txt

info = defaultdict(lambda: defaultdict(lambda: "WT"))
genelist = []
with open("samples.list") as file2:
    for eachline in file2:
        eachline = eachline.split()
        with open(eachline[1]) as file3:
            for eachline1 in file3:
                if re.match(r"^#", eachline1):
                    pass
                else:
                    if re.search(r"\|missense_variant\|", eachline1):
                        gene = re.findall(r"\|missense_variant\|(?:MODERATE|HIGH|LOW)\|(.*?)\|", eachline1)[0]
                        genelist.append(gene)
                        mutationtype = re.findall(r";(HET|HOM);", eachline1)
                        if "HOM" in mutationtype:
                            mutationtype = "HOM"
                        else:
                            mutationtype = "HET"
                        if info[eachline[0]][gene] == "HOM":
                            pass
                        else:
                            info[eachline[0]][gene] = mutationtype
        with open(eachline[2]) as file4:
            file4.readline()
            for eachline2 in file4:
                eachline2 = eachline2.split()
                genelist.append(eachline2[0])
                if float(eachline2[5]) == 0:
                    info[eachline[0]][eachline2[0]] = "NSC"
                else:
                    info[eachline[0]][eachline2[0]]

genelist = list(set(genelist))
genelist.sort()
with open("mutation.info.xls", "w") as file5:
    print >> file5, "#NSC--No Sequence Coverage"
    print >> file5, "#WT--Wild Type"
    print >> file5, "#HET--Heterozygous"
    print >> file5, "#HOM--Homozygous"
    print >> file5, "gene"+"\t"+"\t".join(samplelist)
    for gene in genelist:
        file5.write(gene+"\t")
        for sample in samplelist:
            file5.write(info[sample][gene]+"\t")
        file5.write("\n")

array = []
with open("mutation.info.xls") as file6:
    for i in range(5):
        file6.readline()
    for eachline3 in file6:
        eachline3 = eachline3.split()[1:]
        for i in eachline3:
            if i == "NSC":
                array.append(0)
            elif i == "WT":
                array.append(1)
            elif i == "HET":
                array.append(2)
            elif i == "HOM":
                array.append(3)

ncol = len(samplelist)
nrow = len(genelist)
data = np.array(array).reshape(nrow, ncol)
plt.figure()
ax = plt.subplot(111)
ax.imshow(data, interpolation="nearest", aspect="auto", cmap=cm.Purples)
ax.yaxis.set_ticks_position('right')
ax.yaxis.set_ticks_position('none')
ax.xaxis.set_ticks_position('none')
plt.xticks(range(ncol), samplelist)
plt.yticks(range(nrow), genelist)
plt.axes([0.065, 0.2, 0.05, 0.6])
plt.yticks([0, 1, 2, 3], ["$NSC$", "$WT$", "$HET$", "$HOM$"])
plt.xticks([0], '')
plt.imshow(np.arange(4).reshape(4, 1), interpolation="nearest", aspect="auto", cmap=cm.Purples)
plt.savefig("mutation.pdf")
