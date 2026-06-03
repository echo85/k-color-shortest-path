using JuMP
using DataStructures 

include("ReadInstance.jl")
using .ReadInstance

mutable struct Label
    node::Int
    dist::Float64
    colors::BitSet
    path::Vector{Int}
end

function riduzione_grafo(V_set, A, d_cost, c_color, s_node, t_node, k_max)
    # Costruzione liste di adiacenza provvisorie per il calcolo delle euristiche
    adj = Dict{Int, Vector{Tuple{Int, Float64, Int}}}(u => [] for u in V_set)
    rev_adj = Dict{Int, Vector{Tuple{Int, Float64}}}(u => [] for u in V_set)
    costs = Float64[]
    
    for arco in A
        u, v = arco[1], arco[2]
        push!(adj[u], (v, d_cost[arco], c_color[arco]))
        push!(rev_adj[v], (u, d_cost[arco]))
        push!(costs, d_cost[arco])
    end
    
    # Euristica CCDA (Upper Bound)
    w_min = isempty(costs) ? 0.0 : minimum(costs)
    w_mean = isempty(costs) ? 0.0 : sum(costs) / length(costs)
    w_max = isempty(costs) ? 0.0 : maximum(costs)
    lambdas = [0.0, w_min/4, w_min/2, w_min, 2*w_min, w_mean/4, w_mean/2, w_mean, w_max, 2*w_max]
    
    UB = Inf
    for λ in lambdas
        dist_pen = Dict{Int, Float64}(u => Inf for u in V_set)
        colors_used = Dict{Int, Set{Int}}(u => Set{Int}() for u in V_set)
        parent = Dict{Int, Tuple{Int, Float64, Int}}()
        
        dist_pen[s_node] = 0.0
        pq_ccda = PriorityQueue{Int, Float64}()
        enqueue!(pq_ccda, s_node, 0.0)
        
        while !isempty(pq_ccda)
            u = dequeue!(pq_ccda)
            if u == t_node break end
            
            for (v, cost, col) in adj[u]
                penalized_cost = cost + (col in colors_used[u] ? 0.0 : λ)
                if dist_pen[u] + penalized_cost < dist_pen[v]
                    dist_pen[v] = dist_pen[u] + penalized_cost
                    colors_used[v] = union(colors_used[u], Set([col]))
                    parent[v] = (u, cost, col)
                    pq_ccda[v] = dist_pen[v]
                end
            end
        end
        
        if dist_pen[t_node] < Inf && length(colors_used[t_node]) <= k_max
            curr = t_node
            true_cost = 0.0
            while curr != s_node
                p, c_cost, _ = parent[curr]
                true_cost += c_cost
                curr = p
            end
            UB = min(UB, true_cost)
        end
    end

    # Dijkstra Forward da s_node (distanze minime rilassate dalla sorgente)
    dist_s = Dict{Int, Float64}(u => Inf for u in V_set)
    dist_s[s_node] = 0.0
    pq_s = PriorityQueue{Int, Float64}()
    enqueue!(pq_s, s_node, 0.0)
    
    while !isempty(pq_s)
        u = dequeue!(pq_s)
        for (v, cost, _) in adj[u]
            if dist_s[u] + cost < dist_s[v]
                dist_s[v] = dist_s[u] + cost
                pq_s[v] = dist_s[v]
            end
        end
    end

    # Dijkstra Backward da t_node (calcola h(i) per la riduzione e per A*)
    h = Dict{Int, Float64}(u => Inf for u in V_set)
    h[t_node] = 0.0
    pq_h = PriorityQueue{Int, Float64}()
    enqueue!(pq_h, t_node, 0.0)
    
    while !isempty(pq_h)
        curr = dequeue!(pq_h)
        for (prev, cost) in rev_adj[curr]
            if h[curr] + cost < h[prev]
                h[prev] = h[curr] + cost
                pq_h[prev] = h[prev]
            end
        end
    end

    # Filtraggio degli archi in base all'Upper Bound
    A_ridotto = Vector{Tuple{Int, Int}}()
    for arco in A
        u, v = arco[1], arco[2]
        if dist_s[u] + d_cost[arco] + h[v] <= UB
            push!(A_ridotto, arco)
        end
    end

    #println("   [Riduzione] UB trovato (CCDA): ", UB == Inf ? "Nessuno" : UB)
    println("   [Riduzione] Archi originali: $(length(A)) -> Archi rimanenti: $(length(A_ridotto))")

    # Restituiamo il grafo ridotto, il vettore delle euristiche h e l'UB trovato
    return A_ridotto, h, UB
end

# Funzione per determinare se l1 domina l2 (Definition 2)
function domina(l1::Label, l2::Label)
    # l1 domina l2 se ha costo minore/uguale, è sottoinsieme di colori
    # e almeno una delle due condizioni è strettamente minore
    if l1.dist <= l2.dist && issubset(l1.colors, l2.colors)
        if l1.dist < l2.dist || length(l1.colors) < length(l2.colors)
            return true
        elseif l1.dist == l2.dist && length(l1.colors) == length(l2.colors)
            # A parità assoluta (stesso costo e stessi identici colori), 
            # consideriamo che la prima domini l'altra per evitare percorsi ridondanti.
            return true
        end
    end
    return false
end

function risolvi(V_set, C_set, A, d_cost, c_color, s_node, t_node, k_max, time_limit::Float64=2000.0)
    
    # APPLICAZIONE DELLA RIDUZIONE IN VIA PRELIMINARE
    A_ridotto, h, UB = riduzione_grafo(V_set, A, d_cost, c_color, s_node, t_node, k_max)
    
    if h[s_node] == Inf
        return :infeasible, -1.0, time() - inizio_calcolo
    end
    
    # COSTRUZIONE ADIACENZE SUL GRAFO RIDOTTO
    adj = Dict{Int, Vector{Tuple{Int, Float64, Int}}}(u => [] for u in V_set)
    for arco in A_ridotto
        u, v = arco[1], arco[2]
        push!(adj[u], (v, d_cost[arco], c_color[arco]))
    end
    inizio_calcolo = time()
    # STRUTTURE PER L'ALGORITMO DP A*
    # D Contiene la lista delle un-dominated label per ogni nodo, viene istanziata all'inizio per ogni nodo del Grafo con una label vuota
    D = Dict{Int, Vector{Label}}(u => Label[] for u in V_set)

    # Istanzia la lista di Label da esplorare, all'inizio è sol L_s (label nodo start)
    L = PriorityQueue{Label, Float64}() 
    
    # Inizializzazione nodo sorgente
    l_s = Label(s_node, 0.0, BitSet(), [s_node]) 
    push!(D[s_node], l_s) # Inserimento nelle non dominate del nodo sorgente

    # segue la regola di ordinamento A* f(i) = d_i + h(i)
    # inserimento nella nostra lista con priorità 0 (distanza da nodo sorgente) + h[s_node] (distanza per arrivare al nodo finale senza vincoli di colore)
    enqueue!(L, l_s, l_s.dist + h[s_node]) 
    
    # Inizializziamo il best_cost con l'UB dell'euristica di Cerrone
    best_cost = UB
    # LOOP PRINCIPALE DP
    while !isempty(L)
        tempo_corrente = time() - inizio_calcolo
        if tempo_corrente > time_limit
            # Interrompiamo e restituiamo lo stato di Time Limit
            return :time_limit, -1.0, tempo_corrente
        end

        current_label = dequeue!(L) # Extract(L), estraee la lable with the smallest value according to A*
        
        if current_label.dist < best_cost 
            i = current_label.node 
            if i == t_node # Arrivato al nodo finale
                best_cost = current_label.dist
                break
            else
                # Per ogni nodo j collegato ad i calcoliamo la nuova label (espansione vicini)
                for (j, w_ij, col_ij) in adj[i] 
                    new_dist = current_label.dist + w_ij # Costo nuova distanza
                    new_colors = copy(current_label.colors)
                    push!(new_colors, col_ij) # Inseriamo il nuovo colore dell'arco
                    
                    new_path = copy(current_label.path)
                    push!(new_path, j)

                    # (Algoritmo 2) AddLabel
                    l_j = Label(j, new_dist, new_colors, new_path)
                    add_label(D,L,l_j, k_max,h)
                end
            end
        end
    end
    
    tempo_totale = time() - inizio_calcolo
    
    if best_cost < Inf
        return :optimal, best_cost, tempo_totale
    else
        return :infeasible, -1.0, tempo_totale
    end
end

function add_label(D,L,l_j,k_max,h)
    if length(l_j.colors) <= k_max # Feasibility check
        is_dominated = false
        labels_to_remove = Int[]
        j = l_j.node
                        
        # Controllo Dominanza nell'insieme D[j]
        for (idx, l_exist) in enumerate(D[j])
            if domina(l_exist, l_j)
                is_dominated = true
                break
            elseif domina(l_j, l_exist)
                push!(labels_to_remove, idx)
            end
        end
        
        # Se la nuova label non è dominata, cancella da D ed L
        if !is_dominated
            for l_rem in labels_to_remove
                filter!(x -> x !== l_rem, D[j])
                if haskey(L, l_rem)
                    delete!(L, l_rem)
                end
            end
            
            push!(D[j], l_j)
            # Aggiunge in coda con priorità definita da A* (f = g + h)
            enqueue!(L, l_j, l_j.dist + h[j])
        end
    end
end

function batch_execution(file_path::String, instance_type::String = "ferrone")
    tempo_tot = 0.0
    risolte   = 0

    if instance_type == "castro" 
        seed_range = 1:60
    else
        seed_range = 27000:27200
    end

    for seed in seed_range
        nome = "$(file_path)$seed.txt"
        !isfile(nome) && continue

        println("Risoluzione di $nome in corso... ")
        if instance_type == "castro"
            V, C, A, d, c, s_node, t_node, k_max = read_castro_instance(nome)
        else
            V, C, A, d, c, s_node, t_node, k_max = read_ferrone_instance(nome)
        end

        status, obj, tempo = risolvi(V, C, A, d, c, s_node, t_node, k_max)

        if status == :optimal
            println("OTTIMO (z = $obj) | Tempo: $(round(tempo, digits=4)) sec")
            tempo_tot += tempo
            risolte   += 1
        elseif status == :time_limit
            println("TIME LIMIT")
        elseif status == :infeasible
            println("INFEASIBLE | Tempo: $(round(tempo, digits=4)) sec")
        else
            println("ERRORE ($status) | Tempo: $(round(tempo, digits=4)) sec")
        end
    end

    if risolte > 0
        println("Istanze ottime trovate : $risolte")
        println("Tempo totale di calcolo: $(round(tempo_tot,   digits=4)) sec")
        println("Average Time           : $(round(tempo_tot/risolte, digits=4)) sec")
    else
        println("\nNessuna istanza risolta con successo (controlla i file o i thread).")
    end
end


#batch_execution("extra/full_ferone-instances/Big/Testing/Grid/Grid_250x500_99700/Grid_250x500_99700_", "ferrone") # time limit  on 27006 
#batch_execution("ferone-instances/R1/Random_75000x750000_112500_", "ferrone")
#batch_execution("ferone-instances/R3/Random_75000x750000_150000_", "ferrone")
#batch_execution("ferone-instances/G1/Grid_100x100_5940_", "ferrone")
#batch_execution("ferone-instances/G3/Grid_100x200_11910_", "ferrone")
#batch_execution("castro-instances/", "castro")
batch_execution("new-instances/", "castro")