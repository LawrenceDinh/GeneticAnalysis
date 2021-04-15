from Bio import SeqIO
from Bio import Entrez
from Bio.Blast import NCBIWWW
from Bio import SearchIO
from lxml import etree
from Bio.Blast import NCBIXML

Entrez.email = "luat.dinh@sjsu.edu"

def readsequences(accession):
    try:
        with Entrez.efetch(db="nucleotide", rettype="gb", retmod="text",
                           id=accession) as request:
            for result in SeqIO.parse(request, "gb"):
                nuc_des = str(result.description)
                nuc_seq = str(result.seq)
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


def processXML(blastresults):
    #(number, acc, key)
    accessionArray ={}
    i = 1
    for hit in blastresults:
        accession = getAccession(hit).id
        # print(hit.id)
        accessionArray[i] = [accession , hit, accession]
        i = i+1
    return accessionArray

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


def readXML(filename):
    blastresults = SearchIO.read(filename, 'blast-xml')
    return blastresults

def createalignmentarray():
    with open('basicblast.xml') as result_handle:
        blast_record = NCBIXML.read(result_handle)
    E_VALUE_THRESH = 0.04
    alignmentarray = []
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
                alignmentarray.append(hsp.sbjct)
    return alignmentarray

def savespecies(resultsdict):
    with open ("speciesarray.txt", "w") as f:
        for key,val in resultsdict.items():
            f.write(val[4] + "\n")
    f.close()

def addtodict(dictionary, array):
    i = 0
    for key, value in arr.items():
        #speciesname = hitinfo(value[0])
        #print(speciesname)
        value.append(array[i])
        #value.append(speciesname)
        print(value)
        i = i+1

def readspecies(filename):
    with open(filename, "r") as f:
        writer = open("sparray.txt","w")
        i = 1
        for line in f:
            sp = line.strip() + "," + str(i) + "\n"
            #print(sp)
            i = i+1
            writer.write(sp)

def finddiffspecies(filename):
    sp = []
    diffsp = []
    results = []
    with open(filename,"r") as f:
        for line in f:
            sp.append(line.strip().split(","))
    for species in sp:
        print(species[0])
        if (species[0] not in diffsp):
            diffsp.append(species[0])
            results.append(species)
    return results

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
    blastresults = readXML(xmlFile)
    #getAlignment(blastresults)
    arr = processXML(blastresults)
    agn = createalignmentarray()
    addtodict(arr, agn)
    print(arr)
    for key, val in arr.items():
        print(val[1])
    #savespecies(arr)
    readspecies("speciesarray.txt")
    uniquespecies = finddiffspecies("sparray.txt")
    for x in uniquespecies:
        print(x)

