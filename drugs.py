#!/usr/bin/env python
# coding: utf-8

"""the parameter one is the information of deleterious snps, 
the two is the path of impact data, always '/mnt/ilustre/app/medical/.data/impact'
the third is the prefix of the output files
example: drugs.py snp.af.deleterious.txt /mnt/ilustre/app/medical/.data/impact snp"""

import re
import sys
import os
from collections import defaultdict

print sys.version
os.system("which python")

def linktag(dgs, links):
    newdgs = []
    for dg in dgs:
        if dg in links.keys():
            newdgs.append('<a href="{0}" target="_blank">{1}</a>'.format(links[dg], dg))
        else:
            newdgs.append(dg)
    return ",".join(newdgs)


def tablerow(*args):
    if len(args) == 4:
        return """
        <tr>
            <td>{0}</td>
            <td class="col-break">{1}</td>
            <td class="col-break">{2}</td>
            <td>{3}</td>
        </tr>""".format(args[0], args[1], args[2], args[3])
    elif len(args) == 5:
        return """
        <tr>
            <td>{0}</td>
            <td>{1}</td>
            <td>{2}</td>
            <td>{3}</td>
            <td>{4}</td>
        </tr>""".format(args[0], args[1], args[2], args[3], args[4])


def avoidshort(le):
    if len(le.split("\t")) >= 16:
        return le.split("\t")[14]
    else:
        return le.split("\t")[9]


data = os.path.abspath(os.path.expanduser(sys.argv[2]))

mdainfo = dict()
kieoinfo = dict()
urls = dict()
snps = defaultdict(lambda: [])

with open(sys.argv[1]) as infile:
    for line in infile:
        snps[line.split("\t")[6]].append(line.strip("\n"))

with open("{0}/MDA.txt".format(data)) as mda:
    for line in mda:
        line = line.strip("\n").split(",")
        mdainfo.setdefault(line[0], line[1:])

with open("{0}/kieo.txt".format(data)) as kieo:
    for line in kieo:
        line = line.strip("\n").split(",")
        kieoinfo.setdefault(line[0], line[1:])

with open("{0}/dsig_hyper.txt".format(data)) as link:
    for line in link:
        line = line.strip("\n").split("\t")
        urls.setdefault(line[0], line[1])

outfile = open("{0}.drugs-prediction.txt".format(sys.argv[3]), "w")
html = open("{0}.drugs-prediction.html".format(sys.argv[3]), "w")

# level 1
print "Level 1 drugs prediction..."
outfile.write("\n\n{0}\n\nLEVEL 1: Actionable Therapeutics\n\n{0}\n\n\n".format("-" * 80))
outfile.write("NCI Match Clinical Trials\n")
outfile.write("Gene\tVariant-Information\tActionable-Therapeutic(s)\tAllele-Frequency\n")
outfile.write("{0}\t{1}\t{2}\t{3}\n".format("-" * 4, "-" * 19, "-" * 25, "-" * 16))

html.write("""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Drugs Prediction</title>
    <style type="text/css">
        h3 {
            font-weight: normal;    
        }
        table {
            width: 100%;
        }
        thead {
            text-align: center;
        }
        .col-break {
            word-wrap: break-word;
            word-break: break-all;
        }
    </style>
</head>
<body>
    <h1 align="center">IMPACT Drug Prediction Analysis</h1>
    <br><hr><br>
    <h2>LEVEL 1: Actionable Therapeutics</h2>
    <h3>NCI Match Clinical Trials</h3>
    <table border="1">
        <thead>
            <tr>
                <th>GENE</th>
                <th>VARIANT</th>
                <th>ACTIONABLE THERAPEUTIC(S)</th>
                <th>ALLELE FREQUENCY</th>
            </tr>
        </thead>""")

if "BRAF" in snps.keys():
    match = 0
    unmatch = []
    for line in snps["BRAF"]:
        if re.search(r"V600E|V600K", line):
            match = 1
            outfile.write("BRAF\t{0}\tDabrafenib,Trametinib\t{1}\n".format(avoidshort(line), line.split()[-1].strip()))
            html.write(tablerow("BRAF", avoidshort(line), linktag(["Dabrafenib", "Trametinib"], urls),
                                line.split()[-1].strip()))
        else:
            unmatch.append(line)
    if match == 0:
        for line in unmatch:
            outfile.write("BRAF\t{0}\tTrametinib\t{1}\n".format(avoidshort(line), line.split()[-1].strip()))
            html.write(tablerow("BRAF", "Deleterious", (linktag(["Trametinib"], urls)), line.split()[-1].strip()))

if "EGFR" in snps.keys():
    for line in snps["EGFR"]:
        if re.search(r"T790M", line):
            outfile.write("EGFR\t{0}\tAZD9291\t{1}\n".format(avoidshort(line), line.split()[-1].strip()))
            html.write(tablerow("EGFR", avoidshort(line), linktag(["AZD9291"], urls), line.split()[-1].strip()))
        else:
            outfile.write("EGFR\t{0}-Deleterious\tAfatinib\t{1}\n".format(avoidshort(line), line.split()[-1].strip()))
            html.write(tablerow("EGFR", "{0}-Deleterious".format(avoidshort(line)), linktag(["Afatinib"], urls),
                                line.split()[-1].strip()))

if "HER2" in snps.keys():
    for line in snps["HER2"]:
        if re.search(r"SNV", line):
            outfile.write("HER2\t{0}-Deleterious\tAfatinib\t{1}\n".format(avoidshort(line), line.split()[-1].strip()))
            html.write(tablerow("HER2", "{0}-Deleterious".format(avoidshort(line)), linktag(["Afatinib"], urls),
                                line.split()[-1].strip()))
        elif re.search(r"AMP", line):
            outfile.write("HER2\t{0}-AMPLIFICATION\tAdo-trastuzumab_emtansine\t{}\n".format(avoidshort(line),
                                                                                           line.split()[-1].strip()))
            html.write(tablerow("HER2", "{0}-AMPLIFICATION".format(
                avoidshort(line)), linktag(["Ado-trastuzumab_emtansine"], urls), line.split()[-1].strip()))

if "KIT" in snps.keys():
    for line in snps["KIT"]:
        outfile.write("KIT\t{0}-Deleterious\tSunitinib\t{1}\n".format(avoidshort(line), line.split()[-1].strip()))
        html.write(tablerow("KIT", "{0}-Deleterious".format(avoidshort(line)), linktag(["Sunitinib"], urls),
                            line.split()[-1].strip()))

if "NF2" in snps.keys():
    for line in snps["NF2"]:
        if re.search(r"DEL", line):
            outfile.write("NF2\t{0}-DELETION\tVS6063\t{1}\n".format(avoidshort(line), line.split()[-1].strip()))
            html.write(tablerow("NF2", "{0}-DELETION".format(avoidshort(line)), linktag(["VS6063"], urls),
                                line.split()[-1].strip()))

for info in [mdainfo, kieoinfo]:
    if info == mdainfo:
        outfile.write("\nMD Anderson Personalized Cancer Therapy\n")
        outfile.write("Gene\tVariant-Information\tActionable-Therapeutic(s)\tAllele-Frequency\n")
        outfile.write("{0}\t{1}\t{2}\t{3}\n".format("-" * 4, "-" * 19, "-" * 25, "-" * 16))
        html.write("""
    </table>
    <h3>MD Anderson Personalized Cancer Therapy</h3>
    <table border="1">
        <thead>
            <tr>
                <th>GENE</th>
                <th>VARIANT</th>
                <th>ACTIONABLE THERAPEUTIC(S)</th>
                <th>ALLELE FREQUENCY</th>
            </tr>
        </thead>""")
    else:
        outfile.write("\nDsigDB FDA Approved Kinase Inhibitors\n")
        outfile.write("Gene\tVariant-Information\tActionable-Therapeutic(s)\tAllele-Frequency\n")
        outfile.write("{0}\t{1}\t{2}\t{3}\n".format("-" * 4, "-" * 19, "-" * 25, "-" * 16))
        html.write("""
    </table>
    <h3>DSigDB FDA Approved Kinase Inhibitors</h3>
    <table border="1">
        <thead>
            <tr>
                <th>GENE</th>
                <th>VARIANT</th>
                <th>ACTIONABLE THERAPEUTIC(S)</th>
                <th>ALLELE FREQUENCY</th>
            </tr>
        </thead>""")
    for gene in info.iterkeys():
        if gene in snps.keys():
            for line in snps[gene]:
                if re.search(r"p\.\w+\d+\w+\s", line):
                    outfile.write("{0}\t{1}\t{2}\t{3}\n".format(gene, avoidshort(line), ",".join(info[gene]),
                                                            line.split()[-1].strip()))
                    html.write(tablerow(gene, avoidshort(line), linktag(info[gene], urls), line.split()[-1].strip()))

# level 2
print "Level 2 drugs prediction..."
outfile.write("\n\n{0}\n\nLEVEL 2: Actionable Therapeutics from DsigDB Database\n\n{0}\n\n".format("-" * 80))
outfile.write("Potential Gene Targets({0})\n{1}\n".format(len(snps.keys()), " ".join(sorted(snps.keys()))))
html.write("""
    </table>
    <h2>LEVEL 2: Actionable Therapeutics from DSigDB</h2>
    <p>Potential Gene Targets({0}): {1}</p>
    <table border="1">
        <thead>
            <tr>
                <th>DRUG</th>
                <th>TARGET HIT</th>
                <th>POTENTIAL TARGETS</th>
                <th>P-value(Hypergeometric Test)</th>
                <th>P-value(Permutation Test)</th>
            </tr>
        </thead>""".format(len(snps.keys()), " ".join(sorted(snps.keys()))))
genes = [x for x in snps.keys() if x in [y.strip("\n") for y in open("{0}/D1_geneList.txt".format(data)).readlines()]]
drugs = defaultdict(lambda: {"genes": None, "info": None})
unique = []
with open("{0}/DSigDB_D1_data_set.txt".format(data)) as dsigdb:
    for line in dsigdb:
        line = line.strip("\n").split("\t")
        for gene in line[2:]:
            if gene in genes:
                drugs[line[0]]["genes"] = line[2:]
                unique.extend(line[2:])

unique = list(set(unique))

outfile.write("Drug\tTargets-Hit\tPotential-Targets\tP-value(hypergeometric test)\tP-value(Permutation test)\n")

hypertest = open("{0}.hyperTest.txt".format(sys.argv[3]), "w")

druglist = drugs.keys()
druglist = sorted(druglist, key=lambda j: j.lower())
for drug in druglist:
    match = 0
    for gene in drugs[drug]["genes"]:
        if gene in genes:
            match += 1
    hypertest.write("{0}\t{1}\t{2}\n".format(match, len(genes), len(drugs[drug]["genes"])))
    drugs[drug]["info"] = (match, len(drugs[drug]["genes"]))
hypertest.close()

print "Calculate p values..."
os.system(
    "sh /mnt/ilustre/app/medical/tools/script/PermRHypertest.sh {0}.hyperTest.txt > {0}.p_value.txt".format(sys.argv[3])
)

permp = []
hyper = []
with open("{0}.p_value.txt".format(sys.argv[3])) as pvalue:
    for line in pvalue:
        hyper.append(line.split("\t")[-2])
        permp.append(line.split("\t")[-1].strip("\n"))

druglist2 = [k.strip("\n") for k in open("{0}/druglist2.txt".format(data)).readlines()]

count = 1
for drug in druglist:
    if drug in druglist2:
        outfile.write("{0}\t{1}\t{2}\t{3:.3f}\t{4}\n".format(
            drug, drugs[drug]["info"][0], drugs[drug]["info"][1], float(hyper[count]), permp[count]
        ))
        html.write(tablerow(linktag([drug], urls), drugs[drug]["info"][0], drugs[drug]["info"][1],
                            "{0:.3f}".format(float(hyper[count])), permp[count]))
    count += 1

html.write("""
    </table>
</body>
</html>""")
outfile.close()
html.close()
