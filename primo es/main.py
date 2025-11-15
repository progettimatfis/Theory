# ALGORITMO PER LA RICERCA DELL'AUTOVALORE MINIMO E AUTOVETTORE
import numpy as np
from funzioni import annealing,main
A,b = main()
'''
A = [[0, 1, 0],
     [1, 0, 1],
     [0, 1, 0]]
'''
dim_mat = len(A) # perchè A quadrata
traccia = 0
for i in range (dim_mat) :
    traccia = traccia + A[i][i]
autoval = traccia * 1 / dim_mat
H = A - autoval*np.eye(dim_mat)
# b = 5

# costruzione p vettore precisione
p = np.zeros(b)
p[0] = -1
for i in range(1,b) :
    p[i]=1/(2**i)

# vogliamo Psx*H*Pdx

Psx = np.zeros((dim_mat * b, dim_mat))
# Riempimento matrice con i pt a sx di H = A-lambdaIn
Psx = np.zeros((dim_mat * b, dim_mat))
for j in range(dim_mat) :
    i = j*b
    k = 0
    for i in range(i,(j+1)*b) :
        Psx[i][j] = p[k]
        k = k+1
print('Psx : \n', Psx)
# Riempimento matrice con i p a dx di H = A-lambdaIn
Pdx = np.zeros((dim_mat, dim_mat * b))
for i in range(dim_mat) :
    j = i*b
    k = 0
    for j in range(j,(i+1)*b) :
        Pdx[i][j] = p[k]
        k = k+1
print('Pdx : \n', Pdx)

vett_reale = annealing(dim_mat, b, p, Psx, H, Pdx)
# normalizzazione
vett_reale = vett_reale/np.linalg.norm(vett_reale)

print('vett_reale \n' ,vett_reale)

xtAx = vett_reale.T @ A @ vett_reale
tol = 1e-5
maxiter = 500
k=0
while  (abs(autoval - xtAx) > tol) and k< maxiter :
    print('autoval : ',autoval)
    print('xtAx : ',xtAx)
    print('normavett : ',np.linalg.norm(vett_reale))
    autoval = xtAx / (np.linalg.norm(vett_reale)**2)
    autovett = vett_reale
    H = A - autoval * np.eye(dim_mat)
    vett_reale = annealing(dim_mat, b, p, Psx, H, Pdx)
    vett_reale = vett_reale/ np.linalg.norm(vett_reale)
    xtAx = vett_reale.T @ A @ vett_reale
    k = k+1
print('auto_vett\n : ',autovett)
print ('auto_val\n : ',autoval)

# confronto con python
[D, V] = np.linalg.eig(A)
print(D)
print(V)
lambda_min = D[0]
print('errore sull\'autovalore\n : ',abs(lambda_min-autoval))
print('errore sull\'autovettore\n : ',np.linalg.norm(autovett-V[:,0]))