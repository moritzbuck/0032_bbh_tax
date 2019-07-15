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
ncbi.update_taxonomy_database()

mags = [ f for f in os.listdir() if f.endswith(".faa")]

ibosomals = {m[:-4] : [(" ".join(s.description.split()[1:]), s.seq) for s in SeqIO.parse(m,"fasta") if "S ribosomal protein " in s.description]  for m in tqdm(mags)}

ibosomals = {k : { vv : [q[1] for q in v if q[0] == vv] for vv in {p[0] for p in v} } for k, v in ibosomals.items()}

raw_fasta = { k : [SeqRecord(seq=s,id = k + "_" + kk.replace(" ","_") + "_" + str(i) ) for kk, vv in v.items() for i,s in enumerate(vv) ]  for k,v in ibosomals.items()}
for k, v in raw_fasta.items():
    SeqIO.write(v, k + "_ibosomal.fasta", "fasta")

# annot_set = set(sum([ list(v.keys()) for k,v in ibosomals.items()],[]))
# annot_set = {k for k in annot_set if "lase" not in k and "type" not in k and "Altern" not in k and "transf" not in k }
#
#
# ordered_ibosomals = {g : { k : ibosomals[k][g] for k in ibosomals if ibosomals[k].get(g)} for g in annot_set}
# ordered_ibosomals = {k : v for k,v in ordered_ibosomals.items() if len(v) > 3 }
# ordered_ibosomals = {k.replace(" ","_") : v for k, v in ordered_ibosomals.items()}
# ordered_ibosomals = {k.replace("/","_") : v for k, v in ordered_ibosomals.items()}
#
# SSC_ribos = { k : len(v) for k, v in ordered_ibosomals.items() if max([len(vv) for vv in v.values()]) == 1}
#
# sc_ordered_ibosomals = {k : ordered_ibosomals[k] for k,v in SSC_ribos.items() if v > 6 }
#
# aligned_scs = {}
# for k, v in sc_ordered_ibosomals.items():
#     if not os.path.exists("ibosomals/" + k):
#         os.makedirs("ibosomals/" + k )
#     recs = [SeqRecord(id = kk + "_{numb}_{count}".format(numb=i+1, count=len(vv)), seq = vvv, description = "") for kk, vv in v.items() for i,vvv in enumerate(vv)]
#     fasta = "ibosomals/" + k + "/seqs.faa"
#     SeqIO.write(recs, fasta, "fasta")
#     os.system("muscle -in {inp} -out {out}".format(inp = fasta, out = fasta.replace(".faa", ".ali.faa")))
#     os.system("fasttree < {inp} > {out}".format(inp = fasta.replace(".faa", ".ali.faa"), out = fasta.replace(".faa", ".tree")))
#     aligned_scs[k] = { "_".join(s.id.split("_")[:-2]) : s for s in SeqIO.parse(fasta.replace(".faa", ".ali.faa"), "fasta")}
#
# sc_mags = list(set(sum([list(v.keys()) for v in sc_ordered_ibosomals.values()],[])))
# sc_genes = list(sc_ordered_ibosomals.keys())

# concat_ali = {s : Seq("") for s in sc_mags}
# for g in sc_genes:
#     glen = len(list(aligned_scs[g].values())[0])
#     for m in sc_mags:
#         gen = aligned_scs[g].get(m)
#         if gen:
#             gen = gen.seq
#         else :
#             gen = Seq("-"*glen)
#         concat_ali[m] += gen
#
# ali_len = len(list(concat_ali.values())[0])
#mask = [ [s[i] for s in concat_ali.values()].count("-")/len(concat_ali) < 0.5 for i in range(ali_len) ]
#masked_ali = { k : Seq("".join([ss for ss,m in zip(str(s), mask) if m ])) for k,s in concat_ali.items()}
#masked_ali_len = len(list(masked_ali.values())[0])
# SeqIO.write([SeqRecord(id = k, seq=v, description ="") for k, v in concat_ali.items()],"concat.ali.fasta", "fasta")
# os.system("fasttree < {inp} > {out}".format(inp = "concat.ali.fasta", out = "concat.tree"))


SeqIO.write(sum(raw_fasta.values(),[]), "ribosomal_proteins.faa", "fasta")
os.system("diamond blastp --query all_ibosomal.faa --db $DIAMOND_UNIREF90 --threads 20 --more-sensitive --out all_ibosomal.uniref.more-sensitive.diamond --outfmt 6")

#files = [d + "/" + d.replace("run_", "full_table_") + ".tsv" for d in os.listdir() if d.startswith("run_")]

#diamond = pandas.read_csv("all_ibosomal.env-nr.more-sensitive.diamond", sep="\t", header=None)
#diamond2 = pandas.read_csv("all_ibosomal.more-sensitive.diamond", sep="\t", header=None)
diamond = pandas.read_csv("all_ibosomal.uniref.more-sensitive.diamond", sep="\t", header=None)

#diamond = pandas.concat([diamond, diamond2, diamond3])
diamond = diamond.loc[diamond[10] < 10**-10]
diamond['genome'] = [t.split("0S_r")[0][:-2] for t in diamond[0]]
best_hits = {f[0] : list(f[1].loc[f[1][2] == max(f[1][2])][1])[0] for f in diamond.groupby(0)}
to_tax = set(best_hits.values())
genomes = set([t.split("0S_r")[0][:-2]  for t in best_hits.keys()])

unis = {u for u in best_hits.values() if u.startswith('UniRef')}
to_tax = to_tax.difference(unis)

with open("/home/moritz/dbs/UniRef90.tax") as handle:
    uni_tax = {l.split()[0] : l.split()[1] for l in tqdm(handle) if l.split()[0] in unis}

def get_taxa(i):
    lineage = ncbi.get_lineage(i)
    names = ncbi.get_taxid_translator(lineage)
    names = {k : v for k , v in names.items() if v not in ['root', 'cellular organisms']}

    ranks = ncbi.get_rank(lineage)
#    ranks = {k : v for k , v in ranks.items() if v  != "no rank"}

    return [names[l] for l in lineage if names.get(l)]

uni_tax_map = {k : get_taxa(v) for k, v in uni_tax.items()}


Entrez.email = "murumbii@gmail.com"
entry =[]
handle = Entrez.efetch(db="protein", id=list(to_tax), rettype = 'gb', retmode='text')
entry += [s for s in tqdm(SeqIO.parse(handle,"genbank"))]
tax_map = {s.id : s.annotations['taxonomy'] for s in entry }
tax_map.update(uni_tax_map)

tax_facts = {g : [tax_map[v] for k,v in best_hits.items() if k.startswith(g + "_")] for g in genomes}

def sum_tax(vect, level, clean = True):
    taxa = [v[level] if len(v) > level else None for v in vect ]
    if clean:
        taxa = [t for t in taxa if t not in ['unclassified sequences', 'metagenomes', 'ecological metagenomes']]
    counts = {t : taxa.count(t) for t in set(taxa)}
    top = max(counts, key=lambda l : counts[l])
    return (top, counts[top]/len(taxa))


bestclasses = { g :  [sum_tax(tax_facts[g], i) for i in range(8) ] for g in genomes}
tax_out = {k : ";".join([ "{name}({nub:.2f})".format(name = vv[0], nub=vv[1]) for vv in v if vv[0] and vv[1] >0.4]) for k,v in bestclasses.items()}
