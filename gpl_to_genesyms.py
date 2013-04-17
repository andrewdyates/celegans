#!/usr/bin/python
"""Convert GPL200 to probe, gene symbol list.
Account for _s_ probe multiples.

ID_REF	GB_ACC	GI	ORF	SPOT_ID	Species Scientific Name	Annotation Date	Sequence Type	Sequence Source	Representative Public ID	Target Description	Gene Title	Gene Symbol	ENTREZ_GENE_ID	RefSeq Transcript ID	WormBase	Gene Ontology Biological Process	Gene Ontology Cellular Component	Gene Ontology Molecular Function
0:ID_REF: probe name
12:gene symbol; format "lin-26 /// lir-1 /// WBGene00003012"

USE:
python gpl_to_genesyms.py ../GSE2180_GPL200.probes.tab > ../GSE2180_GPL200.probe2genesym.tab
"""
import sys

def main(fname):
  fp = open(fname)
  gene_to_probes = {}
  fp.next() # skip header
  for line in fp:
    d = line.split('\t')
    try:
      probe, genes = d[0], d[12].split(' /// ')
    except:
      print >> sys.stderr, line
      sys.exit(1)
    if "_x_" in probe:
      continue
    gene = ";".join([g for g in genes if not g.startswith("WBGene")])
    if len(gene) == 0:
      continue
    gene_to_probes.setdefault(gene, set()).add(probe)

  for g in sorted(gene_to_probes):
    for p in gene_to_probes[g]:
      print "%s\t%s" % (p, g)

if __name__ == "__main__":
  main(sys.argv[1])
