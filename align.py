from numpy import *
# sequence alignment script
# it uses linear programming
# written by Jose Flores-Canales 08/2015

# define values for PAM matrix
data_matrix=[2,-2,0,0,-2,0,0,1,-1,-1,-2,-1,-1,-4,1,1,1,-6,-4,0,
6,0,-1,-4,1,-1,-3,2,-2,-3,3,0,-5,0,0,-1,2,-4,-3,
2,2,-4,1,1,0,2,-2,-3,1,-2,-4,-1,1,0,-4,-2,-2,
4,-5,2,3,1,1,-2,-4,0,-3,-6,-1,0,0,-7,-4,-2,
12,-5,-5,-3,-3,-2,-6,-5,-5,-4,-3,0,-2,-8,0,-2,
4,3,-1,3,-2,-2,1,-1,-5,0,-1,-1,-5,-4,-2,
4,0,1,-2,-3,0,-2,-5,-1,0,0,-7,-4,-2,
5,-2,-3,-4,-2,-3,-5,-1,1,0,-7,-5,-1,
7,-2,-2,0,-2,-2,0,-1,-1,-3,0,-2,
5,2,-2,2,1,-2,-1,0,-5,-1,4,
6,-3,4,2,-3,-3,-2,-2,-1,2,
5,0,-5,-1,0,0,-4,-4,-2,
6,0,-2,-2,-1,-4,-2,2,
9,-5,-3,-3,0,7,-1,
6,1,0,-6,-5,-1,
2,1,-3,-3,-1,
3,-5,-3,0,
17,0,-6,
10,-3,
4]

# build PAM matrix (lower triangle)
score=zeros((20,20))
a,b=[],[]
for i in xrange(20):
  for j in xrange(i,20,1):
    a.append(j)
    b.append(i)
indices=(array(a),array(b))
score[indices]=data_matrix

# define map between residue name and matrix index
mapkey={'A':0,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8,'I':9,'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19}

# define score Rab
def Rab(a,b):
  if (a!='-'):
    ind1=mapkey[a]
  if (b!='-'):
    ind2=mapkey[b]
  if (a=='-' and b=='-'):
    return 0 
  elif (a=='-' or b=='-'):
    return 0
  else:
    if(ind1>ind2):
      #print a,ind1,b,ind2,score[ind1][ind2]
      return score[ind1][ind2]
    else:
      #print b,ind2,a,ind1,score[ind2][ind1]
      return score[ind2][ind1]

# define linear score Rab
def Linear_Rab(a,b):
  if (a!='-'):
    ind1=mapkey[a]
  if (b!='-'):
    ind2=mapkey[b]
  if (a=='-' and b=='-'):
    return 0 
  elif (a=='-' or b=='-'):
    return -g
  if (a==b):
    return 1
  else:
    return 0

# initialize global variables q,d
# first element (index 0) of q,d is null
q,d=[],[] 

# get length of sequences m,n
def initialize_seq(name1,name2):
  global q
  global d
  qseq=open(name1,'r')
  dseq=open(name2,'r')
  qline=qseq.readline()
  dline=dseq.readline()
  dseq.close()
  qseq.close()
  m=len([c for c in  qline if c.isalpha()])
  n=len([c for c in  dline if c.isalpha()])
  print "m,n : ", m , n

  [q.append(c) for c in qline if c.isalpha()]
  [d.append(c) for c in dline if c.isalpha()]
  return m,n

# construct H matrix of indexes m+1,n+1
def calculate_H():
  H = [[] for i in xrange(m+1)]
  E = [[] for i in xrange(m+1)] # optimal aligment that ends q,-
  G = [[] for i in xrange(m+1)] # optimal aligment that ends -,d
  F = [[] for i in xrange(m+1)] # optimal aligment that ends with score q,d
  for i in xrange(m+1):
    for j in xrange(n+1):
      # H[0][0], E[0],[0],... =0,0,0,0
      H[i].append(0)
      E[i].append(0)
      F[i].append(0)
      G[i].append(0)

  # initialize aligments that end in their corresponding type
  # F[i][0],G[i][0]=-infinity,-infinity
  # E[0][j],G[0][j]=-infinity,-infinity
  for i in xrange(1,m+1):
    H[i][0]=-(g_open+(i-1)*g_extend)
    E[i][0]=-(g_open+(i-1)*g_extend)
    F[i][0]=-999
    G[i][0]=-999
  for j in xrange(1,n+1):
    H[0][j]=-(g_open+(j-1)*g_extend)
    E[0][j]=-999
    F[0][j]=-(g_open+(i-1)*g_extend)
    G[0][j]=-999
  # calculate H
  for i in xrange(1,m+1,1):
    for j in xrange(1,n+1,1):
      E[i][j]=max(E[i-1][j]-g_extend,F[i-1][j]-g_open,G[i-1][j]-g_open)
      F[i][j]=max(E[i][j-1]-g_open,F[i][j-1]-g_extend,G[i][j-1]-g_open)
      G[i][j]=max(E[i-1][j-1]+Rab(q[i-1],d[j-1]),F[i-1][j-1]+Rab(q[i-1],d[j-1]),G[i-1][j-1]+Rab(q[i-1],d[j-1]))
      H[i][j]=max(E[i][j],F[i][j],G[i][j]) # optimal of three types of aligments
  return H,E,F,G

# print matrix H of dimensions m+1,n+1
def printmatrix(name,H,m,n):
  outfile=open(name,'w')
  for i in xrange(0,m+1,1):
    buf = "print >> outfile, \""
    for j in xrange(0,n+1,1):
      buf += str(H[i][j])+" "
    buf+="\""
    exec buf 
  outfile.close()

# B is a matrix of 2,K dimensions, K is at most m+n+1
def CreateB():
 
  B =[[] for i in xrange(2)]
  for i in xrange(2):
    for j in xrange(m+n+1):
      B[i].append(0)
  return B

# receives matrix B and prints string chains
def printmatrixB(name,copyB):
  outfile=open(name,'w')
  line1,line2="",""
  ln1=copyB[0][::-1]
  ln2=copyB[1][::-1]
  t1=[str(i) for i in ln1 if str(i).isalpha() or str(i)=='-']
  t2=[str(i) for i in ln2 if str(i).isalpha() or str(i)=='-']
  for i,j in zip(t1,t2):
    line1+=i
    line2+=j
  print >> outfile, line1
  print >> outfile, line2
  outfile.close()
  return [t1,
          t2]

# get score of B aligment
def getscore(cB):
  m=len(cB[0])
  n=len(cB[1])
  scoreB=0
  if (m!=n):
    print "Dimensions of B rows are different, getscore()\n"
  for j in range(m):
      scoreB+=Rab(cB[0][j],cB[1][j])
  return scoreB


# number of best aligments
occur=0

# backtracking algorithm
def backtrack(i,j,B,k,mode='m'):
  global occur
  if (i==0):
    while(j>0):
      B[0][k],B[1][k]='-',d[j-1]
      k+=1
      j-=1
    #write B
    name='Bmatrix'+str(occur)+'.dat'
    NewB=printmatrixB(name,B)
    #sc=getscore(NewB)
    print "Printing %s, Score %.1f" % (name, H[m][n])
    occur+=1
  elif (j==0):
    while(i>0):
      B[0][k],B[1][k]=q[i-1],'-'
      k+=1
      i-=1
    #write B
    name='Bmatrix'+str(occur)+'.dat'
    print "Printing %s" % (name)
    NewB=printmatrixB(name,B)
    #sc=getscore(NewB)
    print "Printing %s, Score %.1f:" % (name, H[m][n])
    occur+=1
  else:
    if (mode=='u'):
      EE=E[i][j]-g_extend
      FF=F[i][j]-g_open
      GG=G[i][j]-g_open
      ij_state= max(EE,FF,GG)
      
      if (ij_state==EE):
        B[0][k],B[1][k]=q[i-1],'-'
        backtrack(i-1,j,B,k+1,'u')
      if (ij_state==FF):
        B[0][k],B[1][k]='-',d[j-1]
        backtrack(i,j-1,B,k+1,'l')
      if (ij_state==GG):
        B[0][k],B[1][k]=q[i-1],d[j-1]
        backtrack(i-1,j-1,B,k+1,'m')

    if (mode=='l'):
      EE=E[i][j]-g_open
      FF=F[i][j]-g_extend
      GG=G[i][j]-g_open
      ij_state= max(EE,FF,GG)

      if (ij_state==EE):
        B[0][k],B[1][k]=q[i-1],'-'
        backtrack(i-1,j,B,k+1,'u')
      if (ij_state==FF):
        B[0][k],B[1][k]='-',d[j-1]
        backtrack(i,j-1,B,k+1,'l')
      if (ij_state==GG):
        B[0][k],B[1][k]=q[i-1],d[j-1]
        backtrack(i-1,j-1,B,k+1,'m')

    if (mode=='m'):
      EE=E[i][j]
      FF=F[i][j]
      GG=G[i][j]
      ij_state= max(EE,FF,GG)

      if (ij_state==EE):
        B[0][k],B[1][k]=q[i-1],'-'
        backtrack(i-1,j,B,k+1,'u')
      if (ij_state==FF):
        B[0][k],B[1][k]='-',d[j-1]
        backtrack(i,j-1,B,k+1,'l')
      if (ij_state==GG):
        B[0][k],B[1][k]=q[i-1],d[j-1]
        backtrack(i-1,j-1,B,k+1,'m')
####################################################
##
## Main program
##
###################################################

# define constants for affine
g_open = 1 
g_extend = 0.1
print "g_open: %.1f g_extend: %.1f" % (g_open,g_extend)
# enter name of sequences to align
name1='qseq.txt'
name2='dseq.txt'

m,n=initialize_seq(name1,name2)   
H,E,F,G=calculate_H()
printmatrix('Hmatrix.dat',H,m,n)
B=CreateB()
backtrack(m,n,B,1)
