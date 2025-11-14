# genera un numero reale in un intervallo tra due numeri reali con un simulatore reale

from qiskit.visualization import plot_histogram
from qiskit import QuantumCircuit, transpile
import matplotlib.pyplot as plt
import numpy as np
from qiskit_ibm_runtime.fake_provider import FakeFez

# input
a = float(input('inserisci l\'estremo sx dell\'intervallo\' (reale) : \n'))
b = float(input('inserisci l\'estremo dx dell\'intervallo\'(reale) : \n'))

# vettore di precisione
num_bit = 4
p = [(1/2)**i for i in range(1,num_bit+1)]
print(1/sum(p))

# creazione del circuito
qc = QuantumCircuit(num_bit)
for i in range(num_bit) : # i = riga del circuito
    qc.h(i)
qc.measure_all()
print(qc)

simulator = FakeFez()
qc_transpiled = transpile(qc,simulator)
print(qc_transpiled)

shot = 1024
result = simulator.run(qc_transpiled, shots=shot).result()
# per verificare che sia equiprobabile facciamo tanti shots
counts = result.get_counts()
plot_histogram(counts, title='counts')
plt.show()

# risultato della misurazione
result = simulator.run(qc_transpiled, shots=1).result()
# otteniamo il vettore binario misurato:
bit_string = list(result.get_counts().keys())[0] # key -> lista con un solo elemento che è una stringa -> ne prendo il primo elemento che è una stringa
numero_binario = [int(b) for b in bit_string] # stringa -> lista di interi (non posso farlo in place perchè stringa è immutabile)
# per ottenere numero reale nell'intervallo [0,sum(p)]
numero_reale = np.dot(p, numero_binario)
# riportiamo nell'intervallo (a,b)
numero_reale = (numero_reale * 1/sum(p))*(b-a) #da qui otteniamo un numero con tante cifre decimali per via di *1/sum(p)
numero_reale = numero_reale+a
print(f'numero casuale reale tra {a} e {b}: {numero_reale}')
