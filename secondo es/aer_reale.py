# genera un numero reale in un intervallo tra due numeri reali con un simulatore ideale

from qiskit_aer import AerSimulator
from qiskit import QuantumCircuit
import numpy as np

# input
a = float(input('inserisci l\'estremo sx dell\'intervallo\' (reale) : \n'))
b = float(input('inserisci l\'estremo dx dell\'intervallo\'(reale) : \n'))

# vettore di precisione
num_bit = 4
p = [(1/2)**i for i in range(1,num_bit+1)] # avere 0 è inutile , mi da due volte stessi valori,
# così come -1 mi da due volte stessi valori in modulo quindi lo escludo

# creazione del circuito
qc = QuantumCircuit(num_bit)
for i in range(num_bit) : # i = riga del circuito
    qc.h(i)
qc.save_statevector() # funzione d'onda
qc.measure_all()
print(qc)

simulator = AerSimulator(method='statevector')
# NON facciamo traspilazione

# risultato di una misurazione (non importa fare shot perchè gia conosco la prob dallo statevector)
result = simulator.run(qc, shots=1).result()
state = result.get_statevector()
print('state : \n ',state)
# per ottenere il vettore binario misurato
bit_string = list(result.get_counts().keys())[0] # dict_key -> lista con un solo elemento che è una stringa -> ne prendo il primo elemento che è una stringa
numero_binario = [int(b) for b in bit_string] # stringa -> lista di interi (non posso farlo in place perchè stringa è immutabile)
# per ottenere numero reale nell'intervallo [0,sum(p)]
numero_reale = np.dot(p, numero_binario )
# riportiamo nell'intervallo (a,b)
numero_reale = (numero_reale * 1/sum(p))*(b-a)
numero_reale = numero_reale+a
print(f'numero casuale reale tra {a} e {b}: {numero_reale}')

'''
#VERIFICA CHE IL RAGIONAMENTO FILI

# partiamo da 1 perchè altrimenti p = [0,1/2,...] avrei p * [0,1,....] = 0 + 1/2+ ... e anche p * [1,1,....] = 0 +1/2......
print('p : \n',p)
# prodotto tra p e qualsiasi vettore binario di dim = dim(p)
risultati = {}
prodotto=[]
k=0
for bits in product([0, 1], repeat=num_bit):
    print(bits)
    prod = np.dot(p,bits)
    risultati[''.join(map(str, bits))] = prod
    prodotto.append(prod)
    print(prodotto)
    if (k>1) :
        print(float(prodotto[k])-float(prodotto[k-1]))
    k = k+1

# riportiamo i valori in prodotto tra 0 e 1
somma = sum(p)
for i in range(len(prodotto)) :
    prodotto[i] = float (prodotto[i]*1/somma)
    #
    if (i>1) :
        print(float(prodotto[i])-float(prodotto[i-1]))
    #
# per riportarlo nell'intervallo (a,b)
for i in range(len(prodotto)) :
    prodotto[i] = float (prodotto[i]*1/somma)



# x = int(bits, 2)              # converti binario → decimale
'''