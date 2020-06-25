mutate = "DNAfile.txt"

DNAseq = """ TTTTTCTTATTGCTTCTCCTACTGATTATCATAATGGTTGTC
GTAGTGTCTTCCTCATCGCCTCCCCCACCGACTACCACAACGGCTGCCGGAGGGTATTA
CCATCACCAACAGAATAACAAAAAGGATGACGAAGAGTGTTGCTGGCGTCGCCGA
CGGAGTAGCAGAAGGGGTGGCGGAGGG  """
DNA = DNAseq.replace('\n','')


print(len(DNA))
print(DNAseq)

codons = """ 
		'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
		'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
		'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
		'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',				 
		'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
		'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
		'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
		'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
		'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
		'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
		'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
		'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
		'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
		'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
		'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
		'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
	"""

print(codons)
myDNA = input("Enter your DNA: ")
if len(myDNA) < 3:
    print("You have entered an invalid DNA sequence")
    
else:
    print("checking your DNA codon...")
    C = 'G'
    T = 'A'
    A = 'T'
    G = 'C'
    DNAcodonT = print(myDNA.replace("T","A"))
