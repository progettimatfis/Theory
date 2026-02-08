from pyqubo import Binary
import neal
import numpy as np
import random

def inserisci_dati_manualmente(): 
    n = int(input("Inserisci la dimensione della matrice simmetrica A : "))
    A = np.zeros((n,n))

    for i in range(n):
        for j in range(n): 
            A[i, j] = float(input(f"Inserisci A[{i},{j}]: "))
            if (i > j) and (A[i, j] != A[j, i]):
                raise ValueError(f"Errore : la matrice deve essere simmetrica!")

    b = int(input("Inserisci il numero di bit per il vettore di precisione : "))

    print("Scegli la stima iniziale:")
    print("1️⃣  Media autovalori")
    print("2️⃣  Limite di Gershgorin")
    scelta_in = input("\nScegli un'opzione (1 o 2): ").strip()

    if  scelta_in == '1':
        autoval = np.trace(A) /n
    elif scelta_in == '2' :  
        autoval = ger(A,n)
    else : 
        raise ValueError(f"La scelta non è valida, deve essere 1 o 2")
    
    return A, b, scelta_in, autoval

def ger (A,dim_mat) :
    # Lista dei bordi destri dei dischi di Gershgorin
    gershgorin = []
    for i in range(dim_mat):
        R_i = np.sum(np.abs(A[i, :])) - np.abs(A[i,i]) 
        gershgorin.append(A[i,i] - R_i)

    return max(gershgorin)

def genera_dati_random():
    n = random.randint(2, 10)
    b = np.random.choice([4 ,8, 16, 24])
    A = np.random.uniform(-10, 10, size=(n, n)) 
    A = (A + A.T) / 2  # così A è una matrice simmetrica
    lista = { 1: lambda: np.trace(A) / len(A), 2: lambda: ger(A, len(A))}
    scelta_in = random.choice([1, 2])
    autoval = lista[scelta_in]()

    return A, b, scelta_in,autoval

def menu():
    print("== MENU PRINCIPALE ==")
    print("1️⃣  Inserisci dati manualmente")
    print("2️⃣  Genera dati casualmente")
    scelta = input("Scegli un'opzione (1 o 2): ").strip()

    if scelta == '1':
            return inserisci_dati_manualmente()
    elif scelta== '2' : 
            return genera_dati_random()
    else:
        raise ValueError("La scelta non è valida, deve essere 1 o 2")

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
    # In sampleset sono conservate num_reads soluzioni indipendenti del problema generate dal sampler tramite sweeps iterazioni dell'algoritmo

    # Decodifica dei campioni nel sampleset salvati secondo la rappresentazione interna del BQM in configurazioni binarie leggibili secondo le variabili originali del modello PyQUBO
    decoded_samples = model.decode_sampleset(sampleset) 
    # decoded_samples è una lista di oggetti DecodedSolution in cui ogni oggetto ha per attributi la configurazione binaria e l’energia associata, cioè il valore della funzione obiettivo (xt*Q*x) per quella configurazione

    # Ricerca della soluzione migliore
    best = min(decoded_samples, key=lambda s: s.energy) # minimo tra i decoded_samples in base al valore s.energy
    vett_bin_dis = best.sample  # dizionario corrispondente configurazione binaria di best in cui le chiavi sono le variabili binarie del modello PYQUBO 
    vett_bin_ord = np.array([vett_bin_dis[f'x{i}'] for i in range(dim_vett_bin)])  # vettore dei valori di vett_bin_dis corrispondenti alle chiavi in ordine crescente

    # Trasformazione del vettore vett_bin_ord da binario a reale
    vett_reale = np.zeros(dim_mat)
    for i in range(dim_mat):
        vett_reale[i] = np.dot(p, vett_bin_ord[i * b:(i + 1) * b])
    return vett_reale