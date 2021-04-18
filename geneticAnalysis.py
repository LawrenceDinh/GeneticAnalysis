""" +165 lines from collecting basic_data_set
300+ lines """
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib
import matplotlib.pyplot as plot
from readBasicFasta import hitinfo
from Bio.Seq import Seq

def read_fastafile(fastafile):
    return list(SeqIO.parse(fastafile, "fasta"))

def clustalW_alignment(fastafile, alignment_output):
    cline = ClustalwCommandline(r"C:\Program Files (x86)\ClustalW2\clustalw2",
                                infile=fastafile, outfile=alignment_output)
    cline()

def convertClustal_toFASTA(clustalalignment, fastafilename):
    records = SeqIO.parse(clustalalignment, "clustal")
    count = SeqIO.write(records, fastafilename, "fasta")

def makePhyloTree(alignedfile, outfile, matrixfile):
    with open(alignedfile,'r') as a:
        align = AlignIO.read(a,'clustal')
        print(type(align))
        calculator = DistanceCalculator('identity')
        dist_matrix = calculator.get_distance(align)
        c = DistanceTreeConstructor(calculator)
        print(dist_matrix)
        f = open(matrixfile, "w")
        f.write(str(dist_matrix))
        f.close()
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

    count = SeqIO.write(newseq, newalignment, "fasta")

"""Retrieve last 25 sequences from Ortholog file"""
def choose25aligned_withroot(rawalignment, newalignment):
    seq = read_fastafile(rawalignment)
    newseq = seq[-25:]
    newseq.append(seq[0])
    count = SeqIO.write(newseq, newalignment, "fasta")

def readphylofromfile(phylofile):
    tree = Phylo.read(phylofile, 'phyloxml')
    tree.ladderize() # Deeper clades are displayed at top
    return tree

"""Selects 80% of sequences for training set and 20% for query set"""
def choose_hmm_set(filename, training, query):
    seq = read_fastafile(filename)
    length = len(seq)
    max_seq = int(length * 0.8)
    min_seq = int(length - max_seq)
    training_seq = seq[:max_seq]
    query_seq = seq[-min_seq:]

    print(len(training_seq), len(query_seq))
    print(filename, max_seq, min_seq)
    print(type(training_seq))


    count = SeqIO.write(training_seq, training, "fasta")
    count = SeqIO.write(query_seq, query, "fasta")

"""Crops interesting section of DNA from Clustal MSA file"""
def choose_interesting_dna(training, query, start, end, trainingout, queryout):
    seq = read_fastafile(training)
    que = read_fastafile(query)
    seq = cropsequence(seq, start, end)
    que = cropsequence(que, start, end)
    count = SeqIO.write(seq, trainingout, "fasta")
    count = SeqIO.write(que, queryout, "fasta")

def cropsequence(seq, start, end):
    for s in seq:
        crop = str(s.seq)
        crop = crop[start:end]
        interesting = Seq(crop)
        #print(crop)
        s.seq = interesting
        print(type(interesting))
    return seq


if __name__ == "__main__":
    hmm_basic_training_set = "hmm_basic_training_set.fasta"
    hmm_basic_query_set = "hmm_basic_query_set.fasta"
    hmm_related_training_set = "hmm_related_training_set.fasta"
    hmm_related_query_set = "hmm_related_query_set.fasta"
    croppedhmm_basic_training_set = "cropped_hmm_basic_training_set2.fasta"
    croppedhmm_basic_query_set = "cropped_hmm_basic_query_set2.fasta"
    croppedhmm_related_training_set = "cropped_hmm_related_training_set.fasta"
    croppedhmm_related_query_set = "cropped_hmm_related_query_set.fasta"
    basicfilename = "convertedbasicalignment.fasta"
    relatedfilename = "convertedrelatedalignment.fasta"
    rootrelatedfilename = "root_convertedrelatedalignment.fasta"
    basic_species_alignedfile = "basicspeciesaligned.aln"
    related_species_alignedfile = "relatedspecies.aln"
    related_species_alignedfile_root = "rooted_relatedspecies.aln"
    basicset_alignment_output = "basicaligned.clustal"
    relatedset_alignment_output = "relatedaligned.clustal"
    relatedset_alignment_output_root = "rooted_relatedaligned.clustal"
    basic_set = "basicset.fasta"
    related_set = "OPN1LW_refseq_transcript.fasta"
    related_set_25 = "relatedset25.fasta"
    related_set_25_root = "rooted_relatedset25.fasta"
    basicphylout = "basictree.xml"
    relatedphylout = "relatedtree.xml"
    rootedrelatedphylout = "rooted_relatedtree.xml"
    basicmatrix = "basicmatrix.txt"
    relatedmatrix = "relatedmatrix.txt"
    related25matrix = "related25matrix.txt"
    """Used to choose 25 sequences from 44 orthologs"""
    choose25aligned_withroot(related_set, related_set_25_root)
    choose25aligned(related_set, related_set_25)
    """Perform ClustalW2"""
    #clustalW_alignment(basic_set, basicset_alignment_output)
    #clustalW_alignment(related_set_25, relatedset_alignment_output)
    #clustalW_alignment(related_set_25_root, relatedset_alignment_output_root)

    """Swap accession number and species name for readability"""
    print('W')
    #write_speciestoalign(basicset_alignment_output, basic_species_alignedfile)
    #write_speciestoalign(relatedset_alignment_output, related_species_alignedfile)
    #write_speciestoalign(relatedset_alignment_output_root, related_species_alignedfile_root)

    """Convert clustal to fasta format """
    #convertClustal_toFASTA(basicset_alignment_output, basicfilename)
    #convertClustal_toFASTA(relatedset_alignment_output, relatedfilename)
    #convertClustal_toFASTA(relatedset_alignment_output_root, rootrelatedfilename)

    """Create phylogenetic tree from .aln file"""
    #basic_phylotree = makePhyloTree(basic_species_alignedfile, basicphylout, basicmatrix)
    basic_phylotree = readphylofromfile(basicphylout)
    Phylo.draw_ascii(basic_phylotree)
    #displayPhyloTree(basic_phylotree)

    #related_phylotree = makePhyloTree(related_species_alignedfile, relatedphylout, relatedmatrix)
    related_phylotree = readphylofromfile(relatedphylout)
    #displayPhyloTree(related_phylotree)
    Phylo.draw_ascii(related_phylotree)

    root_related_phylotree = makePhyloTree(related_species_alignedfile_root, rootedrelatedphylout, related25matrix)
    root_related_phylotree = readphylofromfile(rootedrelatedphylout)
    # displayPhyloTree(related_phylotree)
    Phylo.draw_ascii(root_related_phylotree)
    basic = read_fastafile("relatedset25.fasta")
    i = 1
    for set in basic:
        print(set.id)
    """Prepare dataset for HMM"""
    #choose_hmm_set(basicfilename, hmm_basic_training_set, hmm_basic_query_set)
    #choose_hmm_set(relatedfilename, hmm_related_training_set, hmm_related_query_set)

    """Choose chunk of interesting sequence to be processed by HMM"""
    choose_interesting_dna(hmm_basic_training_set, hmm_basic_query_set, 920, 1420, croppedhmm_basic_training_set, croppedhmm_basic_query_set)
    #choose_interesting_dna(hmm_related_training_set, hmm_related_query_set, 400, 900, croppedhmm_related_training_set, croppedhmm_related_query_set)