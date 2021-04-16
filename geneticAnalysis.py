""" +165 lines from collecting basic_data_set
diff = 90 lines """

from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline

basic_set = "basicset.fasta"
related_set = "OPN1LW_refseq_transcript.fasta"
basicset_alignment_output = "basicaligned.fasta"
relatedset_alignment_output = "relatedaligned.fasta"

def clustalW_alignment(fastafile, alignment_output):
    basic_data_set = list(SeqIO.parse(fastafile, "fasta"))
    related_data_set = list(
        SeqIO.parse("OPN1LW_refseq_transcript.fasta", "fasta"))

    cline = ClustalwCommandline(r"C:\Program Files (x86)\ClustalW2\clustalw2",
                                infile=fastafile, outfile=alignment_output)
    cline()
    print(cline)

if __name__ == "__main__":
    clustalW_alignment(basic_set, basicset_alignment_output)
    clustalW_alignment(related_set, relatedset_alignment_output)