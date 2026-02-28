import sys
from Bio import SeqIO #pip install biopython
from Bio.Seq import Seq

#TODO put in an object for better defined scopes
seqOne= ""
seqTwo = ""
seqOneLen = 0
seqTwoLen = 0
file_path = ""
opening = True
traceback = True

def printMatrix(seqOne, seqTwo, matrix):
    print(list(str(seqOne)))
    for x in range(len(matrix)):
        print(f"{seqTwo[x]} {matrix[x]}")
    print(seqOneLen, seqTwoLen)

if len(sys.argv) > 1:
    file_path = sys.argv[1]
    print(f"File is: {file_path}.")
else:
    print("Missing Filename / input file instead?")
    file_path = input("input file:")
    if file_path == "1": # for debuging
        seqOne = Seq(" GAT")
        seqOneLen = len(seqOne)
        seqTwo = Seq(" ACT")
        seqTwoLen = len(seqTwo)
        opening = False

while opening:
    try: 
        with open(file_path, "r") as handle:
            opening = False
            count = 0
            for record in SeqIO.parse(handle, "fasta"):
                if count == 0:
                    seqOne = Seq(" " + str(record.seq)) #hacky solution, but used for the base csaes
                    seqOneLen = len(record)
                elif count == 1:
                    seqTwo = Seq(" " + str(record.seq))
                    seqTwoLen = len(record)
                count += 1        
                print(record)
    except FileNotFoundError:   
        print("Error / input file instead?")
        file_path = input("input file:")
        if file_path == "1": # for debuging
            print("EXAMPLE_ONE")
            seqOne = Seq(" GAT")
            seqTwo = Seq(" ACC")
            seqOneLen = len(seqOne)
            seqTwoLen = len(seqTwo)
            opening = False

#seqOne = str(seqOne)
#seqTwo = str(seqTwo)

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


printMatrix(seqOne, seqTwo, matrix)


#TODO: clean up magic numbers (indexs) use x,y
for i in range(1, seqTwoLen):
  for j in range(1, seqOneLen):
    if ((seqOne[j] == 'A') & (seqTwo[i] == 'G')) | ((seqOne[j] == 'T') & (seqTwo[i] == 'C')):
      pairwise_score = matrix[i-1][j-1]+TRANSITION
    elif (seqOne[j] == seqTwo[i]):
      pairwise_score = matrix[i-1][j-1]+IDENTITY
    else:
      pairwise_score = matrix[i-1][j-1]+TRANSVERSION
    pairwise_gap = max(matrix[i-1][j]+GAP, matrix[i][j-1]+GAP)
    printMatrix(seqOne, seqTwo, matrix)
    print(i,j)
    matrix[i][j] = max(pairwise_score, pairwise_gap)
    pairwise_score = float('-inf')
    pairwise_gap = float('-inf')

printMatrix(seqOne, seqTwo, matrix)

#make the traceback from bottom right

seqOneScore = ""
seqTwoScore = ""
x = seqOneLen -1
y = seqTwoLen -1
trace = []

while traceback: #would rather just use pointers atp
    #if(x == seqOneLen | y = seqTwoLen):
        


    if (x!=0) & (y!=0): #bounds safety
        if (matrix[x-1][y-1] > matrix[x-1][y]) & (matrix[x-1][y-1] > matrix[x][y-1]): #case diagonal
            seqOneScore = seqOne[y] + seqOneScore
            seqTwoScore = seqTwo[x] + seqTwoScore
            x-=1
            y-=1
        elif (matrix[x-1][y] > matrix[x][y-1]): #gap horizontal
            seqOneScore = "-" + seqOneScore
            seqTwoScore = seqTwo[x] + seqTwoScore 
            x-=1
        else: #gap vertical 
            seqOneScore = seqTwo[x] + seqOneScore
            seqTwoScore = "-" + seqTwoScore
            y-=1
    elif (x==0):
        seqOneScore = seqTwo[x]+ seqOneScore 
        seqTwoScore = "-" + seqTwoScore
        y-=1
    else: #only veritcals remain
        seqOneScore = "-" + seqOneScore 
        seqTwoScore = seqTwo[x] + seqTwoScore
        x-=1
    trace.append((x,y))
    if (x == 0) & (y == 0):
        print(trace)
        traceback = False

#Print score and matching sequences
print(seqOneScore, seqTwoScore)
print(matrix[seqOneLen -1][seqTwoLen -1])

#Implement Gotoh's algo in the future