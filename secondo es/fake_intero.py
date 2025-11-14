# genera un numero intero in un intervallo tra due numeri interi con un simulatore reale

from qiskit import QuantumCircuit, generate_preset_pass_manager, transpile
from qiskit.visualization import plot_histogram
from qiskit_ibm_runtime.fake_provider import FakeFez
import matplotlib as mpl

# input
a = int(input('inserisci l\'estremo sx dell\'intervallo\' (intero) : \n'))
b = int(input('inserisci l\'estremo sx dell\'intervallo\'(intero) : \n'))
# lavoriamo prima su [0,b-a] e poi [a,b]
num_bit = (b-a).bit_length()
print(num_bit)

# Creazione circuito.
qc = QuantumCircuit(num_bit)
for i in range(num_bit) : # i = riga del circuito
    qc.h(i)
qc.measure_all()
print(qc)

simulator = FakeFez()
qc_transpiled = transpile(qc,simulator)
print(qc_transpiled)
# verifichiamo probabilità
shot = 1024
result = simulator.run(qc_transpiled, shots=shot).result()
counts = result.get_counts()
plot_histogram(counts, title='counts')
mpl.pyplot.show()

# risultato della misurazione
prop = 1
while( prop == 1 ) :
    result = simulator.run(qc_transpiled, shots=1).result()
    numero_binario = list(result.get_counts().keys())[0]
    numero_decimale = 0
    k = num_bit - 1
    for cifra in numero_binario:
        numero_decimale = numero_decimale + (int(cifra)) * 2 ** (k)
        k = k - 1
    prop = numero_decimale > (b - a)
    if prop == 0 :
        numero_decimale = a + numero_decimale

print(f'numero casuale intero tra {a} e {b}: {numero_decimale}')

# domanda : dobbiamo ottenere una stima delle probabilità nell intervallo a,b e quindi escluder ogni volta
# in cui il circuito mi da un vettore binario che corrisponde a un intero fuori da [a,b]
# oppure mi basta sapere che dal circuito posso ottenre un qualsiasi vettore binario con la stessa proabbilità (poso stare in un intervallo piu grande di [a,b])

# ricorda int bin python