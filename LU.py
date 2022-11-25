import numpy
U = numpy.zeros([5,5])
L = numpy.zeros([5,5])
M = numpy.array([[2, 8, 3, 9, 0, ],
[9, 8, 4, 8, 4, ],
[7, 8, 7, 0, 5, ],
[8, 3, 9, 2, 4, ],
[5, 9, 3, 5, 2, ],
])
for i in range(5):
    for j in range(5):
        if i<=j:
            t=  M[i,j] - sum(L[i,k]*U[k,j] for k in range(i))
            print(f"U[{i},{j}]::",t)
            U[i,j] =t
        else:
            print(f"M[{i},{j}]={M[i,j]}")
            print('+'.join(f"({L[i,k]} * {U[k,j]}:{k})" for k in range(j)))
            t= (M[i,j] - sum(L[i,k]*U[k,j] for k in range(j)))/U[j,j]
            print(f"L[{i},{j}]::",t,"| U[j,j] =",U[j,j])
            L[i,j] =t
print("L:\n",L)
print("U:\n",U)
print("(L+I)*U:\n",(L+numpy.eye(5))@U)
print(mat:=numpy.linalg.solve(L+numpy.eye(5),numpy.array([1,2,3,4,5]).T ))
print(numpy.linalg.solve(U,mat.T))
