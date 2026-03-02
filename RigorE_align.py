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
    for x in range(seqTwoLen):
        print(seqTwo[x], end="   ")
        for y in range(seqOneLen):
            print(matrix[x][y], end="  ")
        print("")
    print("Sequence One length", seqOneLen, "Sequence Two length", seqTwoLen)

if len(sys.argv) > 1:
    file_path = sys.argv[1]
    print(f"File is: {file_path}.")
else:
    print("Missing Filename / input file instead?")
    file_path = input("input file:")
    if file_path == "1": # for debuging
        print("EXAMPLE_ONE")
        seqOne = Seq(" GCC")
        seqOneLen = len(seqOne+1)
        seqTwo = Seq(" ACT")
        seqTwoLen = len(seqTwo+1)
        opening = False

while opening:
    try: 
        with open(file_path, "r") as handle:
            opening = False
            count = 0
            for record in SeqIO.parse(handle, "fasta"):
                if count == 0:
                    seqOne = Seq(" " + str(record.seq)) #hacky solution, but used for the base csaes
                    seqOneLen = len(record.seq)+1
                elif count == 1:
                    seqTwo = Seq(" " + str(record.seq))
                    seqTwoLen = len(record.seq)+1
                count += 1        
                print(record)
    except FileNotFoundError:   
        print("Error / input file instead?")
        file_path = input("input file:")
        if file_path == "1": # for debuging
            print("EXAMPLE_ONE")
            seqOne = Seq(" GAT")
            seqTwo = Seq(" ACC")
            seqOneLen = len(seqOne+1)
            seqTwoLen = len(seqTwo+1)
            opening = False

IDENTITY = 4
TRANSITION = -1 # (a<>g, t<>c)
TRANSVERSION = -2 # (a<>t g<>t c<>a c<>g)
GAP = -10

matrix = [[0 for seqOneDNA in range(seqOneLen)] for seqTwoDNA in range(seqTwoLen)]

pairwise_score = float('-inf')
pairwise_gap = float('-inf')

#base case; 0 * gap should equal 0 for the top-left corner
for x in range(seqOneLen):
    matrix[0][x] = GAP * x
for y in range(seqTwoLen):
    matrix[y][0] = GAP * y

#TODO: clean up magic numbers (indexs)
for x in range(1, seqTwoLen):
  for y in range(1, seqOneLen):
    if (seqOne[y] == seqTwo[x]):
      pairwise_score = matrix[x-1][y-1]+IDENTITY
    elif ((seqOne[y] == 'A') and (seqTwo[x] == 'G')) or ((seqOne[y] == 'T') and (seqTwo[x] == 'C')
        or (seqOne[y] == 'G') and (seqTwo[x] == 'A')) or ((seqOne[y] == 'C') and (seqTwo[x] == 'T')): #use table atp
      pairwise_score = matrix[x-1][y-1]+TRANSITION
    else:
      pairwise_score = matrix[x-1][y-1]+TRANSVERSION
    pairwise_gap = max(matrix[x-1][y]+GAP, matrix[x][y-1]+GAP)
    if "-m" in sys.argv:
        printMatrix(seqOne, seqTwo, matrix)
    matrix[x][y] = max(pairwise_score, pairwise_gap)
    pairwise_score = float('-inf')
    pairwise_gap = float('-inf')

#make the traceback from bottom right
seqOneScore = ""
seqTwoScore = ""
x = seqTwoLen -1
y = seqOneLen -1
trace = [(x,y)]
while traceback and (seqOneLen -1 or seqTwoLen -1): #would rather just use pointers atp
    if (x!=0) and (y!=0): #bounds safety
        if(seqOne[y] == seqTwo[x]): #case match
            seqOneScore = seqOne[y] + seqOneScore
            seqTwoScore = seqTwo[x] + seqTwoScore
            x-=1
            y-=1
        elif (matrix[x-1][y-1] >= matrix[x-1][y]) and (matrix[x-1][y-1] >= matrix[x][y-1]): #case diagonal
            seqOneScore = seqOne[y] + seqOneScore
            seqTwoScore = seqTwo[x] + seqTwoScore
            x-=1
            y-=1
        elif (matrix[x-1][y] > matrix[x][y-1]): #gap horizontal
            seqOneScore = "-" + seqOneScore
            seqTwoScore = seqTwo[x] + seqTwoScore 
            x-=1
        else: #gap vertical 
            seqOneScore = seqOne[y] + seqOneScore
            seqTwoScore = "-" + seqTwoScore
            y-=1
    elif (x==0): #only horizontals remain
        seqOneScore = seqOne[y]+ seqOneScore 
        seqTwoScore = "-" + seqTwoScore
        y-=1
    else: #only veritcals remain
        seqOneScore = "-" + seqOneScore 
        seqTwoScore = seqTwo[x] + seqTwoScore
        x-=1
    trace.append((x,y))
    if (x == 0) and (y == 0):
        print("Sequence one:" + seqOneScore + "\n Sequence two:" + seqTwoScore)
        traceback = False

#Convert matrix into strings for underlined trace
for x in range(1, seqTwoLen):
  for y in range(1, seqOneLen):
      matrix[x][y] = str(matrix[x][y])
      if (x,y) in trace:
          matrix[x][y] = '\033[4m' + matrix[x][y] +'\033[0m'

#Print score and matching sequences

printMatrix(seqOne, seqTwo, matrix)
print("Trace back: ", trace)
print("Score: ", matrix[seqTwoLen -1][seqOneLen -1])

#Implement Gotoh's algo in the future