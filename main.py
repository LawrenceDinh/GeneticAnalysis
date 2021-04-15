from Bio import SeqIO
from Bio import Entrez
from Bio.Blast import NCBIWWW
import lxml.etree as etree
from Bio import SearchIO

Entrez.email = "luat.dinh@sjsu.edu"

def readsequences(accession):
    try:
        with Entrez.efetch(db="nucleotide", rettype="gb", retmod="text",
                           id=accession) as request:
            for result in SeqIO.parse(request, "gb"):
                nuc_des = str(result.description)
                nuc_seq = str(result.seq)
                print(len(nuc_seq))
    except:
        print("Request failed")


def runBlast(GI, hits=50):
    blast_results = NCBIWWW.qblast("blastn", "nt", GI, megablast=True,
                                   hitlist_size=hits)
    return blast_results


def writeBlast(blast_results, filename):
    with open(filename, "w") as out_handle:
        out_handle.write(blast_results.read())
    blast_results.close()


def processXML(filename):
    blastresults = SearchIO.read(filename, 'blast-xml')
    accessionDict = {}
    for hit in blastresults:
        accession = getAccession(hit).id
        print(hit.id)

    return accessionDict

def getAccession(hit):
    accession = hit.id.split('|')
    hit.id = accession[3]
    return hit

def hitinfo(hit):
    namesearch = Entrez.efetch(db="nucleotide", id=hit, rettype="gb",
                           retmode="text")
    genbank = SeqIO.read(namesearch, 'genbank')
    scientificname = genbank.annotations['organism']
    return scientificname

def refinedataset(hitarray):
    print()

"""
An HSP has a query sequence (qseq), hit sequence (hseq) and midline--- '|' 
characters showing matches or a space for no match
"""
if __name__ == "__main__":
    """
    execute blast using the OPN1LW gene as the query sequence,
    only ran for the first time to retrieve basic training set.
    Then saves 100 blast hits to XML format. 
    """
    # writeBlast(runBlast("NM_020061.6", 100), "basicblast.xml")
    xmlFile = "basicblast.xml"
    result_handle = open(xmlFile)

    arr = processXML(xmlFile)
    print(arr)
    for e in arr:
        w = hitinfo(e)
        print(w)

