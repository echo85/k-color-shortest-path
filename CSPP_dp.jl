using JuMP
using DataStructures 

function read_file(file_path::String)
    lines = readlines(file_path)

    # Lettura inteestazione file
    header = split(lines[1])
    num_nodi = parse(Int, header[1])
    s_node = parse(Int, header[3])
    t_node = parse(Int, header[4])
    k_max = parse(Int, header[2])
    
    # Lettura Gradi di ogni nodo dal file
    gradi = zeros(Int, num_nodi)
    for i in 1:num_nodi
        gradi[i] = parse(Int, split(lines[i+1])[1])
    end

    V = Set{Int}(1:num_nodi)
    C = Set{Int}()
    A_set = Set{Tuple{Int, Int}}()
    d = Dict{Tuple{Int, Int}, Float64}()
    c = Dict{Tuple{Int, Int}, Int}()

    linea_corrente = num_nodi + 2
    for u in 1:num_nodi
        for _ in 1:gradi[u]
            parts = split(lines[linea_corrente])
            v = parse(Int, parts[1])
            costo = parse(Float64, parts[2])
            colore = parse(Int, parts[3])
            
            push!(C, colore)
            push!(A_set, (u, v))
            d[(u, v)] = costo
            c[(u, v)] = colore
            linea_corrente += 1
        end
    end
    A = collect(A_set)
    return V, C, A, d, c, s_node, t_node, k_max
end

mutable struct Label
    node::Int
    dist::Float64
    colors::Set{Int}
    path::Vector{Int}
    valid::Bool # Flag per la Lazy Deletion
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

function risolvi(file_path::String; time_limit::Float64=600.0)

    V_set, C_set, A, d_cost, c_color, s_node, t_node, k_max = read_file(file_path)
    
    # 1. COSTRUZIONE ADIACENZE (Normale e Inversa)
    adj = Dict{Int, Vector{Tuple{Int, Float64, Int}}}(u => [] for u in V_set)
    rev_adj = Dict{Int, Vector{Tuple{Int, Float64}}}(u => [] for u in V_set)
    
    for arco in A
        u, v = arco[1], arco[2]
        push!(adj[u], (v, d_cost[arco], c_color[arco]))
        push!(rev_adj[v], (u, d_cost[arco]))
    end
    
    # 2. CALCOLO EURISTICA h(i) (Dijkstra da t_node sul grafo inverso)
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
    # Se t_node non è raggiungibile da s_node nemmeno nel grafo rilassato (senza colori)
    if h[s_node] == Inf
        return MOI.INFEASIBLE, -1.0, 0.0
    end
    
    inizio_calcolo = time()
    
    # 3. STRUTTURE PER L'ALGORITMO 1
    D = Dict{Int, Vector{Label}}(u => Label[] for u in V_set)
    L = PriorityQueue{Label, Float64}()
    
    # Inizializzazione nodo sorgente
    l_s = Label(s_node, 0.0, Set{Int}(), [s_node], true)
    push!(D[s_node], l_s)
    # Il valore di priorità per A* è f(n) = g(n) + h(n)
    enqueue!(L, l_s, l_s.dist + h[s_node]) 
    
    best_cost = Inf
    # 4. LOOP PRINCIPALE DP
    while !isempty(L)
        tempo_corrente = time() - inizio_calcolo
        if tempo_corrente > time_limit
            # Interrompiamo e restituiamo lo stato di Time Limit
            return MOI.TIME_LIMIT, -1.0, tempo_corrente
        end

        current_label = dequeue!(L)
        
        # Ignora le label invalidate (Lazy Deletion)
        if !current_label.valid
            continue
        end
        
        i = current_label.node
        
        # Pruning nel caso avessimo già trovato una soluzione (con A* puro non serve se cerchiamo 1 path, 
        # ma è utile come sistema di guardia)
        if current_label.dist >= best_cost
            continue
        end
        
        # A* garantisce che la prima volta che estraiamo t_node, è la soluzione ottima
        if i == t_node
            best_cost = current_label.dist
            break
        end
        
        # Espansione dei vicini di i
        for (j, w_ij, col_ij) in adj[i]
            # Evita di generare cicli ritornando su nodi già visitati in questo path
            if j in current_label.path
                continue
            end
            
            new_dist = current_label.dist + w_ij
            new_colors = copy(current_label.colors)
            push!(new_colors, col_ij)
            
            # (Algoritmo 2) Feasibility check: controlla se rispetta k_max
            if length(new_colors) <= k_max
                is_dominated = false
                labels_to_remove = Int[]
                
                new_path = copy(current_label.path)
                push!(new_path, j)
                l_j = Label(j, new_dist, new_colors, new_path, true)
                
                # Controllo Dominanza nell'insieme D[j]
                for (idx, l_exist) in enumerate(D[j])
                    if !l_exist.valid
                        continue
                    end
                    if domina(l_exist, l_j)
                        is_dominated = true
                        break
                    elseif domina(l_j, l_exist)
                        push!(labels_to_remove, idx)
                    end
                end
                
                # Se la nuova label non è dominata, aggiorna strutture
                if !is_dominated
                    # Invalida etichette dominate in D[j] e L
                    for idx in labels_to_remove
                        D[j][idx].valid = false 
                    end
                    
                    push!(D[j], l_j)
                    # Aggiunge in coda con priorità definita da A* (f = g + h)
                    enqueue!(L, l_j, l_j.dist + h[j])
                end
            end
        end
    end
    
    tempo_totale = time() - inizio_calcolo
    
    if best_cost < Inf
        return MOI.OPTIMAL, best_cost, tempo_totale
    else
        return MOI.INFEASIBLE, -1.0, tempo_totale
    end
end

function batch_execution(file_path::String)

    tempo_totale_gruppo = 0.0
    istanze_risolte = 0
    
    for seed in 27000:27205
        nome_file = "$(file_path)_$seed.txt"
        if !isfile(nome_file)
            continue
        end
        
        print("Risoluzione di $nome_file in corso... ")
        
        status, obj, tempo = risolvi(nome_file, time_limit=120.0)
        
        if status == MOI.OPTIMAL
            println("OTTIMO (z = $obj) | Tempo: $(round(tempo, digits=4)) sec")
            tempo_totale_gruppo += tempo
            istanze_risolte += 1
        elseif status == MOI.TIME_LIMIT
            println("TIME LIMIT (> 600s) | L'algoritmo è esploso per esplorazione combinatoria.")
        elseif status == MOI.INFEASIBLE
            println("INFEASIBLE | Tempo: $(round(tempo, digits=4)) sec")
        else
            println("ERRORE ($status) | Tempo: $(round(tempo, digits=4)) sec")
        end
    end
    
    if istanze_risolte > 0
        tempo_medio = tempo_totale_gruppo / istanze_risolte
        println("Istanze ottime trovate : $istanze_risolte")
        println("Tempo totale di calcolo: ", round(tempo_totale_gruppo, digits=4), " sec")
        println("Average Time: ", round(tempo_medio, digits=4), " sec")
    else
        println("\nNessuna istanza risolta con successo (controlla i file).")
    end
end

#batch_execution("k-color-instances/Big/Testing/Grid/Grid_500x1000_399400/Grid_500x1000_399400")
#batch_execution("k-color-instances/Big/Testing/Grid/Grid_250x250_49800/Grid_250x250_49800")
#batch_execution("k-color-instances/Big/Testing/Grid/Grid_100x100_7920/Grid_100x100_7920")
#batch_execution("k-color-instances/Small/Grid/Grid_250x500_9970/Grid_250x500_9970")
#batch_execution("k-color-instances/Small/Grid/Grid_500x1000_19970/Grid_500x1000_19970")
#batch_execution("k-color-instances/Small/Grid/Grid_500x1000_39940/Grid_500x1000_39940")
batch_execution("k-color-instances/Big/Testing/Random/Random_125000x2500000_500000/Random_125000x2500000_500000")