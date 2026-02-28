import sys
from Bio import SeqIO #pip install biopython
from Bio import SeqRecord
from Bio.Seq import Seq

sys.argv = ['RigorE_align.py'] #for hardcoding files and such

seqOne = ""
seqTwo = ""
seqOneLen = 0
seqTwoLen = 0
seqCount = 0
file_path = ""
opening = True

if len(sys.argv) > 1:
    file_path = sys.argv[1]
    print(f"File is: {file_path}.")
else:
    print("Missing Filename / input file instead?")
    file_path = input("input file:")
    if file_path == "1": # for debuging
        seqOne = Seq("GAT")
        seqOneLen = len(seqOne)
        seqTwo = Seq("ACCC")
        seqTwoLen = len(seqOne)
        opening = False

while opening:
    try: 
        with open(file_path, "r") as handle:
            opening = False
        for record in SeqIO.parse(handle, "fasta"):
            if count == 0:
                seqOne = record.seq
                seqOneLen = len(record)
            elif count == 1:
                seqTwo = record.seq
                seqOneLen = len(record)
                count += 1        
                print(record)
    except FileNotFoundError:   
        print("Error / input file instead?")
        file_path = input("input file:")
        if file_path == "1": # for debuging
            print("EXAMPLE_ONE")
            seqOne = Seq("GAT")
            seqTwo = Seq("ACCC")
            seqOneLen = len(seqOne)
            seqTwoLen = len(seqTwo)
            opening = False

IDENTITY = 4
TRANSITION = -1 # (a<>g, t<>c)
TRANSVERSION = -2 # (a<>t g<>t c<>a c<>g)
GAP = -10

matrix = [[0 for seqOneDNA in range(seqOneLen)] for seqTwoDNA in range(seqTwoLen)]

pairwise_score = float('-inf')
pairwise_gap = float('-inf')

#base case; 0 * gap should equal 0 for the top-left corner
for i in range(seqOneLen):
    matrix[0][i] = GAP * i
for j in range(seqTwoLen):
    matrix[j][0] = GAP * j
for x in range(len(matrix)):
    print(matrix[x])

'''
#TODO     
for i in range(1, seqOneLen):
  for j in range(1, seqTwoLen):
    if (seqOne[i] == 'A' && seqTwo[j] == 'G') || (seqOne[i] == 'T' && seqTwo[j] == 'C'):
      pairwise_score = seqOne[i-1]+seqTwo[j-1]+TRANSITION
    elif (seqOne[i] == seqTwo[j]):
      pairwise_score = seqOne[i-1]+seqTwo[j-1]+IDENTITY
    else:
      pairwise_score = seqOne[i-1]+seqTwo[j-1]+TRANSVERSION
    pairwise_gap = seqOne[i-1]+seqTwo[j-1]+GAP
    matrix[i][j] = max(pairwise_score, pairwise_gap)
    pairwise_score = float('-inf')
    pairwise_gap = float('-inf')
'''
