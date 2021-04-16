""" +165 lines from collecting basic_data_set
diff = 90 lines """

from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib
import matplotlib.pyplot as plot

basic_set = "basicset.fasta"
related_set = "OPN1LW_refseq_transcript.fasta"
basicset_alignment_output = "basicaligned.clustal"
relatedset_alignment_output = "relatedaligned.clustal"

def read_fastafile(fastafile):
    return list(SeqIO.parse(fastafile, "fasta"))

def clustalW_alignment(fastafile, alignment_output):
    basic_data_set = read_fastafile(fastafile)
    related_data_set = read_fastafile("OPN1LW_refseq_transcript.fasta")

    cline = ClustalwCommandline(r"C:\Program Files (x86)\ClustalW2\clustalw2",
                                infile=fastafile, outfile=alignment_output)
    cline()
    print(cline)

def convertClustal_toFASTA(clustalalignment, fastafilename):
    records = SeqIO.parse(clustalalignment, "clustal")
    count = SeqIO.write(records, fastafilename, "fasta")

def makePhyloTree(alignedfile, outfile):
    with open(alignedfile,'r') as a:
        align = AlignIO.read(a,'clustal')
        print(type(align))
        calculator = DistanceCalculator('identity')
        dist_matrix = calculator.get_distance(align)
        c = DistanceTreeConstructor(calculator)
        phylotree = c.build_tree(align)
        phylotree.rooted = True
        print(phylotree)
        Phylo.write(phylotree, outfile, "phyloxml")
        return phylotree

def displayPhyloTree(phylotree):
    f = Phylo.draw(phylotree)

if __name__ == "__main__":
    basicfilename = "convertedbasicalignment.fasta"
    relatedfilename = "convertedrelatedalignment.fasta"
    #clustalW_alignment(basic_set, basicset_alignment_output)
    #clustalW_alignment(related_set, relatedset_alignment_output)
    with open(basicset_alignment_output, "r") as align:
        alignment = AlignIO.read(align, "clustal")
    i = 0
    for rec in alignment:
        print(rec, i )
        i = i + 1

    convertClustal_toFASTA(basicset_alignment_output, basicfilename)
    convertClustal_toFASTA(relatedset_alignment_output, relatedfilename)
    basic_phylotree = makePhyloTree(basicset_alignment_output, "basictree.xml")
    displayPhyloTree(basic_phylotree)

    related_phylotree = makePhyloTree(relatedset_alignment_output, "relatedtree.xml")
    displayPhyloTree(related_phylotree)