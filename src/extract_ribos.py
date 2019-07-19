import os
from tqdm import tqdm
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
import pandas
from Bio import Entrez
import operator
from ete3 import NCBITaxa

ncbi = NCBITaxa()

with open("/home/moritz/dbs/UniRef90.tax") as handle:
    uni_tax = {l.split()[0] : l.split()[1] for l in tqdm(handle)  if len(l.split())>1}

def get_taxa(i):
    try:
        lineage = ncbi.get_lineage(i)
    except:
        return None
    names = ncbi.get_taxid_translator(lineage)
    names = {k : v for k , v in names.items() if v not in ['root', 'cellular organisms']}

    ranks = ncbi.get_rank(lineage)
#    ranks = {k : v for k , v in ranks.items() if v  != "no rank"}

    return [names[l] for l in lineage if names.get(l)]

def sum_tax(vect, level, clean = True):
    taxa = [v[level] if len(v) > level else None for v in vect ]
    if clean:
        taxa = [t for t in taxa if t not in ['unclassified sequences', 'metagenomes', 'ecological metagenomes']]
    counts = {t : taxa.count(t) for t in set(taxa)}
    top = max(counts, key=lambda l : counts[l])
    return (top, counts[top]/len(taxa))


def process_fasta(fasta, specific = 0, threads=20):
    if specific == 0:
        switch = 'normal'
    elif specific == 1:
        switch = 'more-sensitive'
    else :
        switch = 'very-sensitive'
    file = fasta + "." + switch + ".uniref.diamond"
    if not os.path.exists(file):
        print("running diamond")
        os.system("diamond blastp --query {fasta} --db $DIAMOND_UNIREF90 --threads 20  --out {file} --outfmt 6 {switch} --threads {threads} > /dev/null".format(threads = threads, file = file, fasta = fasta , switch = ("--" + switch) if switch  != 'normal' else ""))
    diamond = pandas.read_csv(file, sep="\t", header=None)
    if len(diamond) == 0:
        [('hit_rate', 0)]
    diamond = diamond.loc[diamond[10] < 10**-10]
    best_hits = {f[0] : list(f[1].loc[f[1][2] == max(f[1][2])][1])[0] for f in diamond.groupby(0)}
    to_tax = set(best_hits.values())
    uni_tax_map = {k : get_taxa(uni_tax.get(k)) for k in tqdm(to_tax) if uni_tax.get(k)}
    uni_tax_map = {k : v for k, v in uni_tax_map.items() if v}
    tax_facts = [uni_tax_map[v] for k,v in best_hits.items() if uni_tax_map.get(v) ]
    bestclasses = [sum_tax(tax_facts, i) for i in range(14) ]
    hit_rate = len(tax_facts)/sum([1 for s in SeqIO.parse(fasta,"fasta")])
    print(bestclasses + [('hit_rate', hit_rate)])
    return bestclasses + [('hit_rate', hit_rate)]
