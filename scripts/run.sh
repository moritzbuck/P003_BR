cd /home/moritz/projects/P003_BR
db_end=/home/moritz/dbs/gtdb.realease89/gtdbtk_sourmash.lca.json

conda activate fastn
fastp -i rawdata/all_reads.fastq -o intermediate/all_reads.fastp.fastq --dont_overwrite
mv fastp.* intermediate/
conda deactivate

conda activate filtlong
filtlong --min_length 1000  intermediate/all_reads.fastp.fastq > intermediate/all_reads.filtlong.fastq
conda deactivate

conda activate fyre
flye -t 20 --nano-raw intermediate/all_reads.filtlong.fastq --meta --out-dir intermediate/flye_coassembly_all -g 100M --plasmids --trestle
conda deactivate

conda activate canu
canu -p canu_assembly_all_reads -d intermediate/canu_assembly_all -nanopore-raw  intermediate/all_reads.filtlong.fastq corOutCoverage=10000 corMhapSensitivity=high corMinCoverage=0 redMemory=32 oeaMemory=32 batMemory=200 correctedErrorRate=0.16
conda deactivate


conda activate medaka
medaka_consensus -i rawdata/all_reads.fastq -d intermediate/flye_coassembly/assembly.racon3.fasta -o intermediate/medaka/ -t 20 -m r941_min_fast_grle
conda deactivate

conda activate minimap2
#minimap2 intermediate/medaka/consensus.fasta intermediate/all_reads.filtlong.fastq -t 20 > intermediate/medaka/selfmap.paf
minimap2 intermediate/medaka/consensus.fasta intermediate/all_reads.fastp.fastq -a   -t 20 > intermediate/medaka/selfmap.sam
for f in `ls rawdata`;
do
  minimap2 intermediate/medaka/consensus.fasta rawdata/$f -a   -t 20 > intermediate/medaka/${f%%.fastq}.sam

  samtools view -u intermediate/medaka/${f%%.fastq}.sam | samtools sort -@ 20 > intermediate/medaka/${f%%.fastq}.bam
  samtools index intermediate/medaka/${f%%.fastq}.bam
  samtools flagstat intermediate/medaka/${f%%.fastq}.bam > intermediate/medaka/${f%%.fastq}.stats
done
conda deactivate

conda activate abricate
  for db in `abricate --list | cut -f1  | tail -n +2`
  do
    for bin in `ls intermediate/medaka/bins`
    do
    abricate --db $db  intermediate/medaka/bins/$bin  > intermediate/medaka/abricate/${bin%%.fna}.abrivate.${db}.tsv
    done
  done

  for db in `abricate --list | cut -f1  | tail -n +2`
  do
    abricate --summary *.abrivate.${db}.tsv | grep -v anvi_genes > summary.abricate.${db}.tsv
  done

deactivate




conda activate anvio-master
source ~/virtual-envs/anvio-master/bin/activate
anvi-script-reformat-fasta consensus.fasta -o contigs4anvio.fa -l 0 --simplify-names
anvi-gen-contigs-database -f contigs4anvio.fa -o anvio/contigs.db -n 'Bioreactor pilot metagenomes'
anvi-run-hmms -c anvio/contigs.db  --also-scan-trnas -T 20
#anvi-run-kegg-kofams -c anvio/contigs.db  -T 20
anvi-get-sequences-for-gene-calls --report-extended-deflines -c anvio/contigs.db --get-aa-sequences -o amino-acid-sequences.fa
emapper.py  -o amino-acid-sequences.fa.emapper.annotations --cpu 20  -m diamond -i amino-acid-sequences.fa
anvi-script-run-eggnog-mapper -c anvio/contigs.db --annotation amino-acid-sequences.fa.emapper.annotations --use-version 2.0.1
sed 's/^\([^#]\)/g\1/' amino-acid-sequences.fa.emapper.annotations > amino-acid-sequences.fa.emapper.annotations.fixed

conda activate sourmash
sourmash compute -k 31 --singleton --scaled 10000 -p 20 contigs4anvio.fa
sourmash lca classify --db $db_end --query  contigs4anvio.fa.sig

for f in `ls maps | grep ".bam$"`
do
  anvi-init-bam maps/$f -o maps/${f%%.bam}.anvio.bam -T20
done

for f in `ls maps | grep ".anvio.bam$"`
do
  anvi-profile -i maps/$f -c anvio/contigs.db --skip-INDEL-profiling --skip-SNV-profiling --sample-name ${f%%.anvio.bam} -o maps/${f%%.anvio.bam} -T 20
done

anvi-merge maps/*/PROFILE.db -o anvio/merged_profiles -c anvio/contigs.db

#make taxo python
"with open('amino-acid-sequences.fa') as handle:
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
from numpy import log10
abricate = [ [l[1], 'abricate_simple:' + l[11], l[12], 'abricate_simple ::' + l[5], '0' if float(l[9])== 100 else str(log10((100-float(l[9]))/100))] for l in abr_lines if not l[0].startswith('#')]
abricate += [ [l[1], 'abricate_long:' + l[11], l[12], 'abricate_long ::' + l[13], '0' if float(l[9])== 100 else  str(log10((100-float(l[9]))/100))] for l in abr_lines if not l[0].startswith('#')]
abricate += [ [l[1], 'abricate_res:' + l[11], l[12], 'abricate_res ::' + l[14], '0' if float(l[9])== 100 else  str(log10((100-float(l[9]))/100))] for l in abr_lines if not l[0].startswith('#') if l[14] != '']
abricate = [[ 'gene_callers_id', 'source' , 'accession', 'function', 'e_value']]+ abricate
with open('abricate4anvio.tsv', 'w') as handle:
  handle.writelines(['\t'.join(l) + '\n' for l in abricate])
"

anvi-import-taxonomy-for-genes -c anvio/contigs.db -i sourmash4anvio.tax -p default_matrix

anvi-run-scg-taxonomy -T 20 -c anvio/contigs.db
anvi-run-kegg-kofams -c anvio/contigs.db  -T 20
anvi-export-collection -p anvio/merged_profiles/PROFILE.db -C bins
for bin in `cut -f2 collection-bins.txt  | sort | uniq`
do
  grep $bin collection-bins.txt | cut -f1 | cut -f1-3 -d_ | sort | uniq > tt
  anvi-export-contigs -c anvio/contigs.db --contigs-of-interest tt -o bins/${bin}.fna
  rm tt
done


anvi-interactive -p anvio/merged_profiles/PROFILE.db  -c anvio/contigs.db
