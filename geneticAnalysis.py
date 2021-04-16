""" +165 lines from collecting basic_data_set
diff = 90 lines """

from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline

basic_set = "basicset.fasta"
related_set = "OPN1LW_refseq_transcript.fasta"
basicset_alignment_output = "basicaligned.clustal"
relatedset_alignment_output = "relatedaligned.aln"

def read_fastafile(fastafile):
    return list(SeqIO.parse(fastafile, "fasta"))

def clustalW_alignment(fastafile, alignment_output):
    basic_data_set = read_fastafile(fastafile)
    related_data_set = read_fastafile("OPN1LW_refseq_transcript.fasta")

    cline = ClustalwCommandline(r"C:\Program Files (x86)\ClustalW2\clustalw2",
                                infile=fastafile, outfile=alignment_output)
    cline()
    print(cline)

def convertClustal_toFASTA(clustalalignment):
    records = SeqIO.parse(clustalalignment, "clustal")
    count = SeqIO.write(records, "convertedbasicalignment.fasta", "fasta")

if __name__ == "__main__":
    clustalW_alignment(basic_set, basicset_alignment_output)
    #clustalW_alignment(related_set, relatedset_alignment_output)
    with open(basicset_alignment_output, "r") as align:
        alignment = AlignIO.read(align, "clustal")
    i = 0
    for rec in alignment:
        print(rec, i )
        i = i + 1
    convertClustal_toFASTA("basicaligned.aln")
