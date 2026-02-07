from pyqubo import Binary
import neal
import numpy as np
import random

def inserisci_dati_manualmente(): 
    n = int(input("Inserisci la dimensione della matrice simmetrica A : "))

    print(" Inserisci gli elementi della matrice A:")
    A = []
    for i in range(n):
        riga = list(map(float, input(f"Riga {i + 1}: ").split()))
        if len(riga) != n:
            raise ValueError("Numero di elementi non corretto nella riga.")
        A.append(riga)
        
    b = int(input("Inserisci il numero di bit per il vettore di precisione : "))

    return A, b

def genera_dati_random():
    n = random.randint(2, 10)
    b = np.random.choice([4 ,8, 16, 24])
    A = np.random.uniform(-10, 10, size=(n, n)) 
    A = (A + A.T) / 2  # così A è una matrice simmetrica

    print("\n Dati generati casualmente:")
    print("Matrice A:\n", A)
    print("Numero di bit per il vettore di precisione : ", b)

    return A, b

def menu():
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
        print("Scelta non valida. Riprova.")
        return menu()

def annealing (dim_mat,b,p, Psx, H, Pdx,num_reads,sweeps) :
    # Costruzione dell'espressione simbolica xtQx 
    Q = Psx @ H @ Pdx
    dim_vett_bin = dim_mat*b
    x = [Binary(f'x{i}') for i in range(dim_vett_bin)] # vettore di variabili binarie simboliche
    expr = 0
    # xtQx = sum(i=1:n) Q_ii * x_i^2 + 2 * sum(i=1:n) sum(j=i+1:n) Q_ij * x_i * x_j perchè Q è simmetrica
    for i in range(dim_vett_bin):
        expr += Q[i,i] * x[i]   # (x[i]**2 == x[i] perchè x[i] variabile binaria)
        for j in range(i+1, dim_vett_bin):
            expr += 2 * Q[i,j] * x[i] * x[j]

    # Trasformazione dell'espressione simbolica 'expr' in un modello PyQUBO (Python Quadratic Unconstrained Binary Optimization)
    model = expr.compile()
    # Conversione del modello PyQUBO in un BQM (Binary Quadratic Model), formato richiesto dal sampler successivo
    bqm = model.to_bqm()

    # Risoluzione del problema BQM
    # Creazione di un sampler, oggetto della classe SimulatedAnnealingSampler, che implementa l’algoritmo di Simulated Annealing per la risoluzione di problemi BQM
    sampler = neal.SimulatedAnnealingSampler() 
    sampleset = sampler.sample(bqm, num_reads=num_reads, sweeps=sweeps) 
    # In sampleset sono conservati num_reads soluzioni indipendenti del problema generate dal sampler tramite sweeps iterazioni dell'algoritmo

    # Decodifica dei campioni nel sampleset salvati secondo la rappresentazione interna del BQM in configurazioni binarie leggibili secondo le variabili originali del modello PyQUBO
    decoded_samples = model.decode_sampleset(sampleset) 
    # decoded_samples è una lista di oggetti DecodedSolution in cui ogni oggetto ha per attributi la configurazione binaria delle variabili 
    # e l’energia associata cioè il valore della funzione obiettivo (xt*Q*x) per quella configurazione

    # Ricerca della soluzione migliore
    best = min(decoded_samples, key=lambda s: s.energy) # minimo tra i decoded_samples in base al valore s.energy
    vett_bin_dis = best.sample  # dizionario corrispondente configurazione binaria di best in cui le chiavi sono le variabili binarie del modello PYQUBO 
    vett_bin_ord = np.array([vett_bin_dis[f'x{i}'] for i in range(dim_vett_bin)])  # vettore dei valori di vett_bin_dis corrispondenti alle chiavi in ordine crescente

    # Trasformazione del vettore vett_bin_ord da binario a reale
    vett_reale = np.zeros(dim_mat)
    for i in range(dim_mat):
        vett_reale[i] = np.dot(p, vett_bin_ord[i * b:(i + 1) * b])
    return vett_reale

def ger (A,dim_mat) :
    # Lista dei bordi destri dei dischi di Gershgorin
    gershgorin = []
    for i in range(dim_mat):
        R_i = np.sum(np.abs(A[i, :])) - np.abs(A[i,i]) 
        gershgorin.append(A[i,i] - R_i)

    return max(gershgorin)
