using JuMP
using DataStructures
include("ReadInstance.jl")
using .ReadInstance

mutable struct Label
    node::Int
    dist::Float64
    colors::BitSet
    parent::Union{Nothing, Label}
end

function in_path(n::Int, l::Label)
    curr = l
    while true
        if curr.node == n
            return true
        end
        if curr.parent === nothing
            break
        end
        curr = curr.parent
    end
    return false
end
function riduzione_grafo(V_set, A, d_cost, c_color, s_node, t_node, k_max)
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

    # Dijkstra Forward da s_node (h_b per ricerca backward)
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

    # Dijkstra Backward da t_node (h_f per ricerca forward)
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

    # Filtraggio degli archi
    A_ridotto = Vector{Tuple{Int, Int}}()
    for arco in A
        u, v = arco[1], arco[2]
        if dist_s[u] + d_cost[arco] + h[v] <= UB
            push!(A_ridotto, arco)
        end
    end

    println("   [Riduzione] UB trovato (CCDA): ", UB == Inf ? "Nessuno" : UB)
    println("   [Riduzione] Archi originali: $(length(A)) -> Archi rimanenti: $(length(A_ridotto))")

    return A_ridotto, h, dist_s, UB # h è h_f (distanza da t), dist_s è h_b (distanza da s)
end

function domina(l1::Label, l2::Label)
    # Filtro Cardinale: se l1 ha più colori di l2, non può dominarlo sui colori
    if l1.dist <= l2.dist && length(l1.colors) <= length(l2.colors)
        if issubset(l1.colors, l2.colors)
            return true # l1 è strettamente uguale o migliore sotto tutti i punti di vista
        end
    end
    return false
end

# Funzione per determinare se l1 domina l2
function domina2(l1::Label, l2::Label)
    if l1.dist <= l2.dist && issubset(l1.colors, l2.colors)
        if l1.dist < l2.dist || length(l1.colors) < length(l2.colors)
            return true
        elseif l1.dist == l2.dist && length(l1.colors) == length(l2.colors)
            return true
        end
    end
    return false
end

function add_label(D, L, l_j, k_max, h)
    if length(l_j.colors) <= k_max
        u = l_j.node
        is_dominated = false
        
        labels_to_remove = Label[] 
        
        for l_exist in D[u]
            if domina(l_exist, l_j)
                is_dominated = true
                break
            elseif domina(l_j, l_exist)
                push!(labels_to_remove, l_exist) # Salviamo l'oggetto
            end
        end
        
        if !is_dominated
            for l_rem in labels_to_remove
                filter!(x -> x !== l_rem, D[u])
                if haskey(L, l_rem)
                    delete!(L, l_rem)
                end
            end
            
            push!(D[u], l_j)
            enqueue!(L, l_j, l_j.dist + h[u])
        end
    end
end

function risolvi(V_set, C_set, A, d_cost, c_color, s_node, t_node, k_max, time_limit::Float64=600.0)
    if s_node == t_node
        return :optimal, 0.0, 0.0
    end

    A_ridotto, h_f, h_b, UB = riduzione_grafo(V_set, A, d_cost, c_color, s_node, t_node, k_max)

    if h_f[s_node] == Inf
        return :infeasible, -1.0, 0.0
    end

    adj_f = Dict{Int, Vector{Tuple{Int, Float64, Int}}}(u => [] for u in V_set)
    adj_b = Dict{Int, Vector{Tuple{Int, Float64, Int}}}(u => [] for u in V_set)
    for arco in A_ridotto
        u, v = arco[1], arco[2]
        push!(adj_f[u], (v, d_cost[arco], c_color[arco]))
        push!(adj_b[v], (u, d_cost[arco], c_color[arco]))
    end

    inizio_calcolo = time()
    
    best_cost = Ref(UB)
    H = UB / 2.0
    
    D_f = Dict{Int, Vector{Label}}(u => Label[] for u in V_set)
    D_b = Dict{Int, Vector{Label}}(u => Label[] for u in V_set)
    
    L_f = PriorityQueue{Label, Float64}()
    L_b = PriorityQueue{Label, Float64}()
    
    l_s = Label(s_node, 0.0, BitSet(), nothing)
    l_t = Label(t_node, 0.0, BitSet(), nothing)
    
    
    push!(D_f[s_node], l_s) 
    enqueue!(L_f, l_s, l_s.dist + h_f[s_node])
    
    push!(D_b[t_node], l_t)
    enqueue!(L_b, l_t, l_t.dist + h_b[t_node])

    # Task Forward
    task_f = Threads.@spawn begin
        while !isempty(L_f)
            tempo_corrente = time() - inizio_calcolo
            if tempo_corrente > time_limit 
                return :time_limit
            end
            
            current_label = dequeue!(L_f)
            if current_label.dist < best_cost[]
                
                i = current_label.node
                
                for (j, w_ij, col_ij) in adj_f[i]
                    if in_path(j, current_label)
                        continue
                    end
                    new_dist = current_label.dist + w_ij
                    if new_dist > H
                        continue # Raggiunto l'Half-Way Point, fermati!
                    end
                    if new_dist + h_f[j] >= best_cost[] continue end
                    
                    new_colors = copy(current_label.colors)
                    push!(new_colors, col_ij)
                    
                    l_j = Label(j, new_dist, new_colors, current_label)
                    add_label(D_f, L_f, l_j, k_max, h_f) # Per il Forward
                end
            end
        end
    end

    # Task Backward
    task_b = Threads.@spawn begin
        while !isempty(L_b)
            tempo_corrente = time() - inizio_calcolo
            if tempo_corrente > time_limit 
                return :time_limit
            end
            
            current_label = dequeue!(L_b)
            if current_label.dist < best_cost[]
                i = current_label.node
                
                for (j, w_ij, col_ij) in adj_b[i]
                    if in_path(j, current_label)
                        continue
                    end
                    new_dist = current_label.dist + w_ij
                    if new_dist + h_b[j]  >= best_cost[] continue end
                    if new_dist > H
                        continue # Raggiunto l'Half-Way Point, fermati!
                    end
                    new_colors = copy(current_label.colors)
                    push!(new_colors, col_ij)
                    
                    l_j = Label(j, new_dist, new_colors, current_label)
                    add_label(D_b, L_b, l_j, k_max, h_b) # Per il Backward
                end
            end
        end
    end

    res_f = fetch(task_f)
    res_b = fetch(task_b)
    if res_f == :time_limit || res_b == :time_limit
        return :time_limit, -1.0, time() - inizio_calcolo
    end
    println("Join Phase...")
    for u in V_set
        
        # Ordina le etichette per costo non decrescente
        Fw = sort(D_f[u], by = x -> x.dist)
        Bw = sort(D_b[u], by = x -> x.dist)
        #println("Best cost $best_cost ...")
        #println("Size FW: "*string(length(Fw)))
        #println("Size BW: "*string(length(Bw)))
        #println(length(Fw))
        # Taglio preventivo (come Step 1 del Pareto-Join)
        if isempty(Fw) || isempty(Bw) continue end
        
        total = Fw[1].dist + Bw[1].dist
        #println("Total Cost: $total")
        #println("Best Cost: $best_cost")
        if Fw[1].dist + Bw[1].dist >= best_cost[] continue end
        
        # Esplorazione con potatura (Pruning)
        for i in 1:length(Fw)
            # Se la forward corrente unita alla migliore backward supera l'incumbent,
            # tutte le successive forward lo supereranno. Interrompiamo il ciclo esterno.
            if Fw[i].dist + Bw[1].dist >= best_cost[]
                break 
            end
            
            for j in 1:length(Bw)
                costo_totale = Fw[i].dist + Bw[j].dist
                
                # Se la somma supera l'incumbent, tutte le successive backward (con questa Fw)
                # costeranno di più. Interrompiamo il ciclo interno (Taglio!).
                if costo_totale >= best_cost[]
                    break
                end
                
                # Check di fattibilità (multi-risorsa: unione dei colori)
                if length(union(Fw[i].colors, Bw[j].colors)) <= k_max
                    best_cost[] = costo_totale
                end
            end
        end
    end
    tempo_totale = time() - inizio_calcolo

    if best_cost[] < Inf
        return :optimal, best_cost[], tempo_totale
    else
        return :infeasible, -1.0, tempo_totale
    end
end

function batch_execution(file_path::String, instance_type::String = "ferrone")
    tempo_tot = 0.0
    risolte = 0
    if instance_type == "castro"
        seed_range = 1:60
    else
        seed_range = 27000:27200
    end

    for seed in seed_range
        nome = "$(file_path)$(seed).txt"
        !isfile(nome) && continue

        println("Risoluzione di $nome in corso...")
        if instance_type == "castro"
            V, C, A, d, c, s_node, t_node, k_max = read_castro_instance(nome)
        else
            V, C, A, d, c, s_node, t_node, k_max = read_ferrone_instance(nome)
        end

        status, obj, tempo = risolvi(V, C, A, d, c, s_node, t_node, k_max)

        if status == :optimal
            println("OTTIMO (z = $obj) | Tempo: $(round(tempo, digits=4)) sec")
            tempo_tot += tempo
            risolte += 1
        elseif status == :time_limit
            println("TIME LIMIT")
        elseif status == :infeasible
            println("INFEASIBLE | Tempo: $(round(tempo, digits=4)) sec")
        else
            println("ERRORE ($status) | Tempo: $(round(tempo, digits=4)) sec")
        end
    end

    if risolte > 0
        println("Istanze ottime trovate: $risolte")
        println("Tempo totale di calcolo: $(round(tempo_tot, digits=4)) sec")
        println("Average Time           : $(round(tempo_tot/risolte, digits=4)) sec")
    else
        println("\nNessuna istanza risolta con successo (controlla i file o i thread).")
    end
end
#batch_execution("extra/full-ferone-instances/Small/Grid/Grid_500x1000_39940/Grid_500x1000_39940_") #G6 good only with 27002, 205 second
#batch_execution("extra/full-ferone-instances/Small/Grid/Grid_100x100_396/Grid_100x100_396_") #G1 works with k/2 colors nad H/2 colors
#batch_execution("extra/full-ferone-instances/Small/Grid/Grid_100x100_792/Grid_100x100_792_") # works well
#batch_execution("extra/full-ferone-instances/Small/Grid/Grid_100x200_794/Grid_100x200_794_") #G2 doesn't work with both solution
#batch_execution("extra/full-ferone-instances/Small/Grid/Grid_100x200_1588/Grid_100x200_1588_") #G2 doesn't work with both solution
#batch_execution("extra/full-ferone-instances/Big/Testing/Grid/Grid_100x100_5940/Grid_100x100_5940_", "ferrone") # works
#batch_execution("extra/full-ferone-instances/Big/Testing/Grid/Grid_250x500_99700/Grid_250x500_99700_", "ferrone") # G4 workds 212 second on 27004

#batch_execution("ferone-instances/R1/Random_75000x750000_112500_", "ferrone")
#batch_execution("ferone-instances/R3/Random_75000x750000_150000_", "ferrone")
#batch_execution("ferone-instances/G1/Grid_100x100_5940_", "ferrone")
#batch_execution("ferone-instances/G3/Grid_100x200_11910_", "ferrone")
#batch_execution("castro-instances/", "castro")
batch_execution("new-instances/", "castro")