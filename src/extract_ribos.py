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
        if len(taxa) == 0:
            return None
    counts = {t : taxa.count(t) for t in set(taxa)}
    top = max(counts, key=lambda l : counts[l])
    return (top, counts[top]/len(taxa))

def run_diamond(fasta, specific = 0, threads = 20, sub = "ibosomal"):
    if specific == 0:
        switch = 'normal'
    elif specific == 1:
        switch = 'more-sensitive'
    else :
        switch = 'very-sensitive'

    if sub != "":
        run_fasta = "temp.fasta"
    else :
        sub = "full"
        run_fasta = fasta

    file = fasta + "." + sub + "." + switch + ".uniref.diamond"

    if not os.path.exists(file + ".done"):
        print("running diamond")
        if sub != "":
            SeqIO.write([s for s in SeqIO.parse(fasta, "fasta") if sub in s.description], "temp.fasta", "fasta")
            if os.stat(run_fasta).st_size == 0:
                return [('hit_rate', -1)]
        os.system("diamond blastp --query {fasta} --db $DIAMOND_UNIREF90 --threads 20  --out {file} --outfmt 6 {switch} --threads {threads} 2> /dev/null > /dev/null".format(threads = threads, file = file, fasta = run_fasta , switch = ("--" + switch) if switch  != 'normal' else ""))
        os.system("touch {file}.done".format(file = file))
    if os.stat(file).st_size != 0:
        return pandas.read_csv(file, sep="\t", header=None)
    else :
        return [('hit_rate', 0)]

def process_single_fasta(fasta, specific = 0, threads=20, sub = "ibosomal"):
    diamond = run_diamond(fasta, specific, threads, sub )
    if len(diamond == 1):
        return diamond

    diamond = diamond.loc[diamond[10] < 10**-10]
    best_hits = {f[0] : list(f[1].loc[f[1][2] == max(f[1][2])][1])[0] for f in diamond.groupby(0)}
    to_tax = set(best_hits.values())
    uni_tax_map = {k : get_taxa(uni_tax.get(k)) for k in tqdm(to_tax) if uni_tax.get(k)}
    uni_tax_map = {k : v for k, v in uni_tax_map.items() if v}
    tax_facts = [uni_tax_map[v] for k,v in best_hits.items() if uni_tax_map.get(v) ]
    if len(tax_facts) > 0:
        bestclasses = [sum_tax(tax_facts, i) for i in range(14) ]
        hit_rate = len(tax_facts)/sum([1 for s in SeqIO.parse(run_fasta,"fasta")])
    #print(bestclasses + [('hit_rate', hit_rate)])
        return bestclasses + [('hit_rate', hit_rate)]
    return [('hit_rate', 0)]

def process_composite_fasta(fasta, specific = 0, threads=20, sub = "ibosomal"):
    diamond = run_diamond(fasta, specific, threads, sub)
    if len(diamond == 1):
        return diamond

    diamond = diamond.loc[diamond[10] < 10**-10]
    best_hits = {f[0] : list(f[1].loc[f[1][2] == max(f[1][2])][1])[0] for f in diamond.groupby(0)}
    to_tax = set(best_hits.values())
    uni_tax_map = {k : get_taxa(uni_tax.get(k)) for k in tqdm(to_tax) if uni_tax.get(k)}
    uni_tax_map = {k : v for k, v in uni_tax_map.items() if v}
    mags = { "_".join(t.split("_")[:-1]) for t in best_hits}
    mags2class = {}
    for mag in tqdm(mags):
        hits = {k : v for k,v in best_hits.items() if "_".join(k.split("_")[:-1]) == mag}
        tax_facts = [uni_tax_map[v] for k,v in hits.items() if uni_tax_map.get(v) ]
        if len(tax_facts) > 0:
            bestclasses = [sum_tax(tax_facts, i, clean = True) for i in range(14) ]
            bestclasses = [t for t in bestclasses if t]
            hit_rate = len(tax_facts)/len(hits)
            mags2class[mag] = bestclasses + [('hits', len(hits), len(tax_facts), hit_rate)]
        else :
            mags2class[mag] = [('hits', len(hits), len(tax_facts), 0)]

    return  mags2class
