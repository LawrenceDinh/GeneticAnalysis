""" +165 lines from collecting basic_data_set
270 lines without comments
"""


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
        #print(phylotree)
        Phylo.write(phylotree, outfile, "phyloxml")
        a.close()
        return phylotree


def displayPhyloTree(phylotree):
    fig = plot.figure(figsize=(13, 7), dpi=100)  # create figure & set the size
    matplotlib.rc('font', size=12)  # fontsize of the leaf and node labels
    axes = fig.add_subplot(1, 1, 1)
    Phylo.draw(phylotree, axes = axes)

def assignname_toalignment(alignments):
    i = 0
    for rec in alignments:
        sciname = hitinfo(rec.id)
        rec.name = rec.id
        rec.id = sciname
        #print(rec)
        i = i + 1

"""Swaps species name with Accession number for better visualization"""
def write_speciestoalign(alignments, phylotree_filename):
    with open(alignments, "r") as align:
        alignment = AlignIO.read(align, "clustal")
        assignname_toalignment(alignment)
        count = SeqIO.write(alignment, phylotree_filename, "clustal")
    align.close()

"""Retrieve last 25 sequences from Ortholog file"""
def choose25aligned(rawalignment, newalignment):
    seq = read_fastafile(rawalignment)
    newseq = seq[-25:]
    newseq.append(seq[0])

    count = SeqIO.write(newseq, newalignment, "fasta")

def readphylofromfile(phylofile):
    tree = Phylo.read(phylofile, 'phyloxml')
    tree.ladderize() # Deeper clades are displayed at top
    return tree

if __name__ == "__main__":
    basicfilename = "convertedbasicalignment.fasta"
    relatedfilename = "convertedrelatedalignment.fasta"
    basic_species_alignedfile = "basicspeciesaligned.aln"
    related_species_alignedfile = "relatedspecies.aln"
    basicset_alignment_output = "basicaligned.clustal"
    relatedset_alignment_output = "relatedaligned.clustal"
    basic_set = "basicset.fasta"
    related_set = "OPN1LW_refseq_transcript.fasta"
    related_set_25 = "relatedset25.fasta"
    basicphylout = "basictree.xml"
    relatedphylout = "relatedtree.xml"
    """Used to choose 25 sequences from 44 orthologs"""
    #choose25aligned(related_set, related_set_25)

    """Perform ClustalW2"""
    #clustalW_alignment(basic_set, basicset_alignment_output)
    #clustalW_alignment(related_set_25, relatedset_alignment_output)

    """Swap accession number and species name for readability"""
    #write_speciestoalign(basicset_alignment_output, basic_species_alignedfile)
    #write_speciestoalign(relatedset_alignment_output, related_species_alignedfile)

    """Convert clustal to fasta format """
    #convertClustal_toFASTA(basicset_alignment_output, basicfilename)
    #convertClustal_toFASTA(relatedset_alignment_output, relatedfilename)

    """Create phylogenetic tree from .aln file"""
    #basic_phylotree = makePhyloTree(basic_species_alignedfile, basicphylout)
    basic_phylotree = readphylofromfile(basicphylout)
    Phylo.draw_ascii(basic_phylotree)
    displayPhyloTree(basic_phylotree)

    #related_phylotree = makePhyloTree(related_species_alignedfile, relatedphylout)
    related_phylotree = readphylofromfile(relatedphylout)
    displayPhyloTree(related_phylotree)
    Phylo.draw_ascii(related_phylotree)