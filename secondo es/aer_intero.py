# genera un numero intero in un intervallo tra due numeri interi con un simulatore ideale

from qiskit_aer import AerSimulator
from qiskit import QuantumCircuit

# input
a = int(input('inserisci l\'estremo sx dell\'intervallo\' (intero) : \n'))
b = int(input('inserisci l\'estremo dx dell\'intervallo\'(intero) : \n'))
# otteniamo prima un numero casuale in [0,(b-a)] e poi trasliamo in [a,b]
num_bit = (b-a).bit_length()

# creazione del circuito
qc = QuantumCircuit(num_bit)
for i in range(num_bit) : # i = riga del circuito
    qc.h(i)
qc.save_statevector()
qc.measure_all()
print(qc)

simulator = AerSimulator(method='statevector') # NON facciamo traspilazione

# risultato della misurazione
prop = 1

while( prop == 1 ) :
    result = simulator.run(qc, shots=1).result()
    numero_binario = list(result.get_counts().keys())  # lista di un unico elemento che è una stringas
    numero_decimale = 0
    k = num_bit - 1
    for cifra in numero_binario[0]:
        numero_decimale = numero_decimale + (int(cifra)) * 2 ** (k)
        k = k - 1
    prop = numero_decimale > (b-a)
    if (prop == 0) :
        state = result.get_statevector() # cosi evita di farlo tante volte se non becca l intervallo [a,b]
        print('state : \n ', state)
        numero_decimale = a + numero_decimale

print(f'numero casuale tra {a} e {b}: {numero_decimale}')