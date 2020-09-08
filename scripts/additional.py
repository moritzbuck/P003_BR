from numpy import log10

with open('amino-acid-sequences.fa') as handle:
  aa2ctg = { l.split('|')[0][1:] : l.split('|')[1].split(':')[1] for l in handle if l.startswith('>')}
ctg2aa = {ctg : [] for ctg in set(aa2ctg.values())}
for aa, ctg in aa2ctg.items():
  ctg2aa[ctg] += [aa]
with open('sourmash.tax') as handle:
  taxo = handle.readlines()
taxo_out = ['gene_callers_id\tt_domain\tt_phylum\tt_class\tt_order\tt_family\tt_genus\tt_species\n']
tax = {t.split(",")[0] :  t.split(",")[2:-1] for t in taxo[1:] if "nomatch" not in t }
taxo_out += [ "\t".join([kk] + t) + "\n"  for k,t in tax.items() for kk in ctg2aa[k]
with open('sourmash4anvio.tax', 'w') as handle:
  handle.writelines(taxo_out)
with open('abricate/anvi_genes.abrivate.all.tsv') as handle:
  abr_lines = [l[:-1].split('\t') for l in  handle.readlines()]

abricate = [ [l[1], 'abricate_simple:' + l[11], l[12], 'abricate_simple ::' + l[5], '0' if float(l[9])== 100 else str(log10((100-float(l[9]))/100))] for l in abr_lines if not l[0].startswith('#')]
abricate += [ [l[1], 'abricate_long:' + l[11], l[12], 'abricate_long ::' + l[13], '0' if float(l[9])== 100 else  str(log10((100-float(l[9]))/100))] for l in abr_lines if not l[0].startswith('#')]
abricate += [ [l[1], 'abricate_res:' + l[11], l[12], 'abricate_res ::' + l[14], '0' if float(l[9])== 100 else  str(log10((100-float(l[9]))/100))] for l in abr_lines if not l[0].startswith('#') if l[14] != '']
abricate = [[ 'gene_callers_id', 'source' , 'accession', 'function', 'e_value']]+ abricate
with open('abricate4anvio.tsv', 'w') as handle:
  handle.writelines(['\t'.join(l) + '\n' for l in abricate])

import os, pandas

tt2 = []
for f in [l for l in os.listdir("abricate/") if 'summary' in l]:
    tt = pandas.read_csv("abricate/" + f, sep ="\t", index_col = 0)
    del tt['NUM_FOUND']
    db = f.split(".")[2]
    tt.columns = [db + ":" + c for c in tt.columns]
    tt.index = [i.split(".")[0] for i in tt.index]
    tt2 += [tt.transpose()]

 full_amr_table = pandas.concat(tt2)

per_gene_table = pandas.read_csv("abricate/all_bins.abrivate.all.tsv", sep ="\t")
id2res = { i[1]['DATABASE'] + ":" + i[1]['GENE'] : i[1]['RESISTANCE'] for i in per_gene_table.iterrows()}
id2prod = { i[1]['DATABASE'] + ":" + i[1]['GENE'] : i[1]['PRODUCT'] for i in per_gene_table.iterrows()}
full_amr_table['resistance'] = [ "" if id2res.get(i,nan) is nan else id2res[i] for i in full_amr_table.index]
full_amr_table['product'] = [ "" if id2prod.get(i,nan) is nan else id2prod[i] for i in full_amr_table.index]
full_amr_table.to_csv("full_amr_presence_table.csv")
reses = set([(i[1]['#FILE'].split("/")[-1][:-4], chem.lower()) for i in per_gene_table.iterrows() if not i[1]['RESISTANCE'] is nan for chem in i[1]['RESISTANCE'].split(";")])

with open("resistances.csv", "w") as handle :
    handle.writelines([f[0] + "," +f[1] + "\n" for f in reses])
