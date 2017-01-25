# this scipt will sort out supplied loci according to their minor allele frequency (MAF)
# used primarily for weighting in Networks analysis (http://www.fluxus-engineering.com/sharenet.htm)


# IMPORTS
# -------
from datetime import datetime


# INPUTS
# ------

# file with mtDNA haplotypes in following form
# name      freq    sequence
# sample1   1       ATTAGA...
# sample2   1       AGTTGA...
INPUT = "mtDNA_haplotypes.txt"

# read the sequences (3rd part of the split line) into 'seq_list' variable
seq_list = []
with open(INPUT, mode="r", encoding="utf-8") as seqIn:
    for line in seqIn:
        line = line.strip().split()
        seq_list.append(line[2])


# provide names of the variable positions/mutations
# if you dont have names, you can use something like: loci = [n for n in range(len(seq_list[0]))]
variable_snps = "0;1;2;3;53;149;150;151;194;198;227;262;315;316;515;724;751;887;957;1009;1440"
loci = variable_snps.split(";")


# FUNCTIONS
# ---------
def minor_allele_frequency(locus, collection):
    """Suggests just 2 alleles per locus
    input: locus [integer]; collection of sequences [something iterable = string, list, tuple...]
    returns: frequency of the minor allele (MAF)"""
    
    snpList = [x[locus] for x in collection]
    maf = snpList.count(snpList[0])/len(snpList)
    #print(maf)
    
    # I dont know which allele is the minor one, so I return one with the lower frequency
    return min(maf, 1 - maf)
    

# SCRIPT
# -------
# cycle through the loci and compute maf for each one

#print(seq_list)
#print("--")
#print(loci)

maf_list = []
    
for n, i in enumerate(loci):
    a = minor_allele_frequency(n, seq_list)
    maf_list.append((int(i), a))


# OUTPUT
# ------
print("=============================")
print("Loci sorted by descending MAF")
print("-----------------------------")
print("INPUT:", INPUT)
print("Timestamp: {:%Y-%m-%d %H:%M:%S}".format(datetime.now()))
print("=============================")
    
# sorted according to the 2nd position (maf) in [(position, maf), (position2, maf2),...]
maf_list.sort(key=lambda x: x[1], reverse=True)

for i in maf_list:
    print(i[0], i[1], sep="\t")