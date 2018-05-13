from sympy.parsing.sympy_parser import parse_expr as parse
from gravipy import *
import sympy as s
import numpy as np
from sympy.matrices import *

#Input of Metric:
print('Please input your line element, keeping in mind general python script, x to the power y is written x**y')
g_input = parse(input('Your line element: ds**2='))      #parse interprets input as function

# Input of Variables:
print('Please input the coordinate system you used in inputting your line element,eg: (t,x,y,z)')
x = parse(input('Your Coordinate System: '))

n = len(x)             #Number of Dimentions
dx = []
g_mn = list([[]for a in range(n)])

for i in range(n):
    dx.append(s.Symbol('d'+str(x[i])))

for i in range(n):
    for j in range(n):
        g1_temp = g_input.diff(dx[i])
        g1_temp = g1_temp.diff(dx[j])
        g1_temp = g1_temp/2
        g_mn[i].append(g1_temp)

print('\n\n\n Metric Calcualted')


#            FINDING INVERSE METRIC
gmn = list(np.array(Matrix(g_mn).inv()))

print('Inverse Metric Calcualted')


#           FINDING DET(g_mn)
g = simplify(Matrix(g_mn).det())
print('Determinant Calcualted')


#               CHRISTOFFEL SYMBOLS:
ch = list([[[]for a in range(n)]for b in range(n)]for b in range(n))
cht = list([[]for a in range(n)]for b in range(n))


for i in range(n):
    for j in range(n):
        for k in range(n):
            for v in range(n):
                ch[i][j][k].append(0.5*gmn[i][v]*(g_mn[k][v].diff(x[j])+g_mn[v][j].diff(x[k])-g_mn[j][k].diff(x[v])))
            cht[i][j].append(simplify(sum(ch[i][j][k])))
print('Christoffel Symbols Calcualted')

#               RIEMANN TENSORS:

R = list([[[[]for a in range(n)]for b in range(n)]for b in range(n)]for c in range(n))
Rp_smv = list([[[]for a in range(n)]for b in range(n)]for b in range(n))


for p in range(n):
    for s in range(n):
        for m in range(n):
            for v in range(n):
                t2arr = []
                t1 = cht[p][v][s].diff(x[m])-cht[p][m][s].diff(x[v])
                for l in range(n):
                    t2 = cht[p][m][l]*cht[l][v][s]-cht[p][v][l]*cht[l][m][s]
                    R[p][s][m][v].append(t1+t2)
                    t2arr.append(t2)
                t2 = sum(t2arr)
                Rp_smv[p][s][m].append(simplify(t1+t2))
print('Riemann Tensors Calcualted')

#                      RICCI TENSOR

R_mn = list([[]for a in range(n)])


for i in range(n):
    for j in range(n):
        p = 0
        for k in range(n):
            p += Rp_smv[k][i][k][j]
#            simplify(p)
        R_mn[i].append(simplify(p))


print('Ricci Tensor Calcualted')

#                   RICCI SCALAR
R = 0
for i in range(n):
    for j in range(n):
        R += g_mn[i][j]*R_mn[i][j]
        # R=simplify(R)                                                     #########################3
print('Ricci Scalar Calcualted')


#                   Einstein TENSOR

Gmn = list([[]for a in range(n)])

for i in range(n):
    for j in range(n):
        G = R_mn[i][j]-0.5*R*g_mn[i][j]
        Gmn[i].append((G))
print('Einstein Tensor Calcualted \n\n\n\n')

print('Would you like the output in Matrix or Latex form?')
m = input('(Input \'Latex\' for latex and \'Matrix\' for matrix: ')
if str(m) == 'latex' or m == 'Latex':
    print('\n\n The metric g_mu,nu is: \n' ,latex(Matrix(g_mn)),'\n\n')
    print('Its Determinant is:  ' ,latex(g),'\n\n')
    print('The Inverse metric: \n' , latex(Matrix(gmn)),'\n\n')
    print('The Christoffel Symbols are:\n ',latex(Matrix(cht[0])),'\n\n')
    for i in range(1,n):
        print(latex(Matrix(cht[i])),'\n')
    print('The Ricci Tensor is:\n ',latex(Matrix(R_mn)),'\n\n')
    print('The Ricci Scalar is: ',latex(R),'\n\n')
    for i in range(n):
        for j in range(n):
            if Gmn[i][j] != 0:
                print('G^{',x[i],x[j],'}=',latex(Gmn[i][j]),'\\')
if str(m) == 'Matrix' or m == 'matrix':
    print('\n\n The metric g_mu,nu is: \n' ,np.array(g_mn),'\n\n')
    print('Its Determinant is:  ' ,g,'\n\n')
    print('The Inverse metric: \n' , np.array(gmn),'\n\n')
    print('The Christoffel Symbols are:\n ',(np.array(cht)),'\n\n')
    print('The Ricci Tensor is:\n ',(np.array(R_mn)),'\n\n')
    print('The Ricci Scalar is: ',R,'\n\n')
    print('The Einstein Tensor is:\n',(np.array(Gmn)),'\n\n')

# Short list of test metrics:
'Spherical Isotropic: '
# -((1-ps/p)/(1+ps/p))**2*c**2*dt**2+(1+ps/p)**4*(dp**2+p**2*(dtheta**2+sin(theta)**2*dphi**2))
'Schwarzschild'
#  -(1-rs/r)*c**2*dt**2+1/(1-rs/r)*dr**2+r**2*((dtheta)**2+sin(theta)**2*dphi**2)
'Cant remember oops'
#-dt**2+s*(dr**2/D+dtheta**2)+(r**2+a**2)*sin(theta)**2*dphi**2+2*m*r/ds*(a*sin(theta)**2*dphi-dt)**2
'Jejjala Schwarzchild:'
# -(1-2*G_N*M/(r*c**2))*c**2*dt**2+(1-2*G_N*M/(r*c**2))**-1*dr**2+r**2*(dtheta**2+sin(theta)**2*dphi**2)
'Kerr Metric'
# -(1-Rs*r/rho**2)*c**2*dt**2+rho**2*dr**2/D+rho**2*dtheta**2+(r**2+alpha**2+Rs*r*alpha**2/rho**2*sin(theta)**2)*sin(theta)**2*dphi**2-2*Rs*r*alpha*sin(theta)**2/rho**2*c*dt*dphi
'Anti de Sitter'
# -(1+r**2/L**2)*c**2*dt**2+(1+r**2/L**2)**-1*dr**2+r**2*(dtheta**2+sin(theta)**2*dphi**2)
