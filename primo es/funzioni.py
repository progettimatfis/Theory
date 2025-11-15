from pyqubo import Binary
import neal
import numpy as np
import random

def inserisci_dati_manualmente(): # SIMMETRICA
    n = int(input("Inserisci la dimensione della matrice simmetrica A : "))

    print("\nInserisci gli elementi della matrice A:")
    A = []
    for i in range(n):
        riga = list(map(float, input(f"Riga {i + 1}: ").split()))
        if len(riga) != n:
            raise ValueError("Numero di elementi non corretto nella riga.")
        A.append(riga)
    A = np.array(A)

    b = int(input("Inserisci il numero di bit: "))

    return A, b

def genera_dati_random():
    """Genera casualmente A (matrice), b (scalare) e numero di bit."""
    n = random.randint(2, 7)
    b = int(np.random.choice([8, 16, 32, 64]))
    A = np.random.randint(-10, 10, size=(n, n))  # genera numeri interi casuali
    A = (A + A.T) / 2  # rende la matrice simmetrica

    print("Matrice simmetrica A:")
    print(A)

    print("\n✅ Dati generati casualmente:")
    print("Matrice A:\n", A)
    print("Numero di bit:", b)

    return A, b

def main():
    print("== MENU PRINCIPALE ==")
    print("1️⃣  Inserisci dati manualmente")
    print("2️⃣  Genera dati casualmente")

    scelta = input("\nScegli un'opzione (1 o 2): ").strip()
    if scelta == '1':
        A, b = inserisci_dati_manualmente()
        return A,b
    elif scelta == '2':
        A, b = genera_dati_random()
        return A,b
    else:
        print("❌ Scelta non valida. Riprova.")
        return main()

    print("\n=== RISULTATI FINALI ===")
    print("Matrice A:\n", A)
    print("Valore b:", b)

def annealing (dim_mat,b,p, Psx, H, Pdx) :
    dim_vett_bin = dim_mat*b
    Q = Psx @ H @ Pdx
    x = [Binary(f'x{i}') for i in range(dim_vett_bin)]

    # === 3. Costruzione dell’espressione x^T Q x ===
    expr = 0
    for i in range(dim_vett_bin):
        for j in range(dim_vett_bin):
            expr += Q[i, j] * x[i] * x[j]

    # === 4. Compila il modello in QUBO ===
    model = expr.compile()
    # compile() converte l’espressione simbolica (expr) in una struttura che
    # PyQUBO può manipolare per creare un Binary Quadratic Model (BQM).
    bqm = model.to_bqm()
    # to_bqm() converte il modello PyQUBO
    # in un Binary Quadratic Model (BQM) standard.

    # === 5. Risolvi con Simulated Annealing ===
    sampler = neal.SimulatedAnnealingSampler()
    sampleset = sampler.sample(bqm, num_reads=100, sweeps=2000)
    # 100 prove indipendenti

    # === 6. Decodifica i risultati ===
    decoded_samples = model.decode_sampleset(sampleset)  # da pyqubo a binario (come abbiamo def x)
    # cerco il minimo tra i decoded_samples in base al valore s.energy (ottengo un oggetto la cui configurazione è binaria)
    best = min(decoded_samples, key=lambda s: s.energy)
    # best = oggetto che tra tre attributi :
    # best.sample → valori ottimali delle variabili binarie
    # best.energy → valore minimo della funzione obiettivo
    # best.occurrence → quante volte è stata trovata questa soluzione
    # per questo printiamo solo best.sample
    print(" Migliore soluzione trovata:")
    vett_bin_dis = best.sample  # dimensione = nb (è un dizionario) # vettore binario disordinato
    print(vett_bin_dis)

    # for i in range(dim_vett_bin) :
    #    print (vett_bin_dis[f'x{i}'])
    # riordino vett_bin_dis (py_qubo non lo restituisce ordinato)
    vett_bin_ord = np.array([vett_bin_dis[f'x{i}'] for i in range(dim_vett_bin)])  # dimensione = nb # vettore binario ordinato

    # trasformazione binario->reale
    vett_reale = np.zeros(dim_mat)
    for i in range(dim_mat):
        print ('lunghezza vett bin ord : ', len (vett_bin_ord[i * b:(i + 1) * b]) )
        vett_reale[i] = np.dot(p, vett_bin_ord[i * b:(i + 1) * b])
    print('vett_reale \n', vett_reale)
    return vett_reale