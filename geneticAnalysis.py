""" +165 lines from collecting basic_data_set
diff = 90 lines """

from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib
import matplotlib.pyplot as plot
from readBasicFasta import hitinfo


def read_fastafile(fastafile):
    return list(SeqIO.parse(fastafile, "fasta"))

def clustalW_alignment(fastafile, alignment_output):
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
    fig = plot.figure(figsize=(13, 5), dpi=100)  # create figure & set the size
    matplotlib.rc('font', size=12)  # fontsize of the leaf and node labels
    matplotlib.rc('xtick', labelsize=10)  # fontsize of the tick labels
    matplotlib.rc('ytick', labelsize=10)  # fontsize of the tick labels
    # turtle_tree.ladderize()		   # optional way to reformat your tree
    axes = fig.add_subplot(1, 1, 1)
    f = Phylo.draw(phylotree, axes = axes)

def assignname_toalignment(alignments):
    i = 0
    for rec in alignments:
        sciname = hitinfo(rec.id)
        rec.name = rec.id
        rec.id = sciname
        print(rec)
        i = i + 1

"""Swaps species name with Accession number for better visualization"""
def write_speciestoalign(alignments, phylotree_filename):
    with open(alignments, "r") as align:
        alignment = AlignIO.read(align, "clustal")
        assignname_toalignment(alignment)
        count = SeqIO.write(alignment, phylotree_filename, "clustal")

def choose25aligned(rawalignment):
    seq = read_fastafile(rawalignment)
    i = 0
    print(seq[:10])
    for hit in seq:
        #print(hit, i)
        i = i+1


if __name__ == "__main__":
    basicfilename = "convertedbasicalignment.fasta"
    relatedfilename = "convertedrelatedalignment.fasta"
    basic_species_alignedfile = "basicspeciesaligned.aln"
    related_species_alignedfile = "relatedspecies.aln"
    basicset_alignment_output = "basicaligned.clustal"
    relatedset_alignment_output = "relatedaligned.clustal"
    basic_set = "basicset.fasta"
    related_set = "OPN1LW_refseq_transcript.fasta"
    basicphylout = "basictree.xml"
    relatedphylout = "relatedtree.xml"
    #clustalW_alignment(basic_set, basicset_alignment_output)
    #clustalW_alignment(related_set, relatedset_alignment_output)
    choose25aligned(related_set)
    #write_speciestoalign(basicset_alignment_output, basic_species_alignedfile)
    convertClustal_toFASTA(basicset_alignment_output, basicfilename)
    #convertClustal_toFASTA(relatedset_alignment_output, relatedfilename)
    basic_phylotree = makePhyloTree(basic_species_alignedfile, basicphylout)
    #displayPhyloTree(basic_phylotree)


    #related_phylotree = makePhyloTree(related_species_alignedfile, relatedphylout)
    #displayPhyloTree(related_phylotree)
    read_fastafile(relatedfilename)