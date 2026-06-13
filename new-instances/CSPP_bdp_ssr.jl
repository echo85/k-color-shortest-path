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
    rev_adj = Dict{Int, Vector{Tuple{Int, Float64, Int}}}(u => [] for u in V_set)
    costs = Float64[]
    for arco in A
        u, v = arco[1], arco[2]
        push!(adj[u], (v, d_cost[arco], c_color[arco]))
        push!(rev_adj[v], (u, d_cost[arco], c_color[arco]))
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
    # 1. Calcolo Frequenze dei Colori sugli archi rimasti
    color_freq = Dict{Int, Int}()
    for arco in A_ridotto
        c = c_color[arco]
        color_freq[c] = get(color_freq, c, 0) + 1
    end
    # Array dei colori ordinati per frequenza decrescente
    sorted_colors = sort(collect(keys(color_freq)), by = x -> color_freq[x], rev=true)

    # 2. Calcolo E_min (Salti minimi verso t_node)
    E_min = Dict{Int, Int}(u => typemax(Int) for u in V_set)
    E_min[t_node] = 0
    pq_e = PriorityQueue{Int, Int}()
    enqueue!(pq_e, t_node, 0)

    while !isempty(pq_e)
        curr = dequeue!(pq_e)
        for (prev, _, _) in rev_adj[curr] # Assumendo che rev_adj abbia (nodo, costo, colore)
            if E_min[curr] + 1 < E_min[prev]
                E_min[prev] = E_min[curr] + 1
                pq_e[prev] = E_min[prev]
            end
        end
    end

    print("[Rid CCDA]", UB == Inf ? "N\\A" : UB)
    print(" Archi:$(length(A))->$(length(A_ridotto)) ")

    return A_ridotto, h, dist_s, UB, E_min, color_freq, sorted_colors

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
function calcola_min_nuovi_colori(E_min_j, new_colors, color_freq, sorted_colors)
    if E_min_j == typemax(Int) || E_min_j == 0
        return 0
    end
    
    # Quanti archi possiamo coprire al massimo "gratis" con i colori già presi?
    max_cover_gratis = 0
    for c in new_colors
        max_cover_gratis += get(color_freq, c, 0)
    end
    
    E_rimanenti = E_min_j - max_cover_gratis
    if E_rimanenti <= 0
        return 0 # I colori attuali potrebbero bastare
    end
    
    # Se avanzano salti, dobbiamo prendere nuovi colori (simuliamo il caso migliore prendendo i più frequenti)
    nuovi_colori_necessari = 0
    archi_coperti = 0
    
    for c in sorted_colors
        if !(c in new_colors)
            archi_coperti += color_freq[c]
            nuovi_colori_necessari += 1
            if archi_coperti >= E_rimanenti
                break
            end
        end
    end
    
    return nuovi_colori_necessari
end
function extend(j, d_j, new_colors, l_i, k_max, h, best_cost, H, E_min, color_freq, sorted_colors)
    # Calcolo SSR: Colori minimi obbligatori
    C_new_min = calcola_min_nuovi_colori(E_min[j], new_colors, color_freq, sorted_colors)
    
    # Pruning di Ammissibilità (Feasibility Bound)
    if length(new_colors) + C_new_min > k_max
        return nothing 
    end

    if d_j + h[j] >= best_cost 
        return nothing 
    end
    
    if d_j > H 
        return nothing 
    end

    if length(new_colors) <= k_max
        l_j = Label(j, d_j, new_colors, l_i)
        return l_j
    end
    
    return nothing
end
function add_label(D, Gamma_bar, l_j::Nothing)
    return nothing
end

function add_label(D, Gamma_bar, l_j::Label)
    u = l_j.node
    is_dominated = false
    labels_to_remove = Label[] 
    
    for l_exist in D[u]
        if domina(l_exist, l_j)
            is_dominated = true
            break
        elseif domina(l_j, l_exist)
            push!(labels_to_remove, l_exist)
        end
    end
    
    if !is_dominated
        for l_rem in labels_to_remove
            filter!(x -> x !== l_rem, D[u])
            # Rimuoviamo l'etichetta dominata anche da quelle da estendere
            filter!(x -> x !== l_rem, Gamma_bar[u]) 
        end
        
        push!(D[u], l_j)
        push!(Gamma_bar[u], l_j)
    end
end

function risolvi(V_set, C_set, A, d_cost, c_color, s_node, t_node, k_max, time_limit::Float64=600.0)
    if s_node == t_node
        return :optimal, 0.0, 0.0
    end

    A_ridotto, h_f, h_b, UB, E_min, color_freq, sorted_colors = riduzione_grafo(V_set, A, d_cost, c_color, s_node, t_node, k_max)

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
    
    best_cost = UB
    H = UB / 2.0 
    D_f = Dict{Int, Vector{Label}}(u => Label[] for u in V_set)
    D_b = Dict{Int, Vector{Label}}(u => Label[] for u in V_set)
    
    Gamma_bar_f = Dict{Int, Vector{Label}}(u => Label[] for u in V_set)
    Gamma_bar_b = Dict{Int, Vector{Label}}(u => Label[] for u in V_set)
    E_queue = Int[s_node, t_node] # Vettore usato come coda
    
    l_s = Label(s_node, 0.0, BitSet(), nothing)
    l_t = Label(t_node, 0.0, BitSet(), nothing)
    
    push!(D_f[s_node], l_s)
    push!(Gamma_bar_f[s_node], l_s)
    
    push!(D_b[t_node], l_t)
    push!(Gamma_bar_b[t_node], l_t)

    while !isempty(E_queue)
        tempo_corrente = time() - inizio_calcolo
        if tempo_corrente > time_limit 
            return :time_limit, -1.0, tempo_corrente
        end
        
        i = popfirst!(E_queue)
        
        if !isempty(Gamma_bar_f[i])
            sort!(Gamma_bar_f[i], by = x -> length(x.colors))
        end
        # Forward extension
        for l_i in copy(Gamma_bar_f[i])
            for (j, w_ij, col_ij) in adj_f[i]    
                new_colors = copy(l_i.colors)
                push!(new_colors, col_ij)                     
                l_j = extend(j, l_i.dist + w_ij, new_colors, l_i, k_max, h_f, best_cost, H, E_min, color_freq, sorted_colors)
                add_label(D_f, Gamma_bar_f, l_j)
                if !isempty(Gamma_bar_f) push!(E_queue, j) end
            end
        end
        empty!(Gamma_bar_f[i])
        
        # Backward extension
        if !isempty(Gamma_bar_b[i])
            sort!(Gamma_bar_b[i], by = x -> length(x.colors))
        end
        for l_i in copy(Gamma_bar_b[i])
            for (j, w_ij, col_ij) in adj_b[i]              
                new_colors = copy(l_i.colors)
                push!(new_colors, col_ij)
                
                l_j = extend(j, l_i.dist + w_ij, new_colors, l_i, k_max, h_b, best_cost, H, E_min, color_freq, sorted_colors)
                add_label(D_b, Gamma_bar_b, l_j)
                if !isempty(Gamma_bar_b) push!(E_queue, j) end
            end
        end
        empty!(Gamma_bar_b[i])
    end

    # JOIN
    for u in V_set
        Fw = sort(D_f[u], by = x -> x.dist)
        Bw = sort(D_b[u], by = x -> x.dist)
        
        if isempty(Fw) || isempty(Bw) continue end
        if Fw[1].dist + Bw[1].dist >= best_cost continue end
        
        for i in 1:length(Fw)
            if Fw[i].dist + Bw[1].dist >= best_cost break end
            
            for j in 1:length(Bw)
                costo_totale = Fw[i].dist + Bw[j].dist
                if costo_totale >= best_cost break end
                
                if length(union(Fw[i].colors, Bw[j].colors)) <= k_max
                    best_cost = costo_totale
                end
            end
        end
    end
    
    tempo_totale = time() - inizio_calcolo

    if best_cost[] < Inf
        return :optimal, best_cost, tempo_totale
    else
        return :infeasible, -1.0, tempo_totale
    end
end

function batch_execution(file_path::String, instance_type::String = "ferrone")
    tempo_tot = 0.0
    risolte = 0
    if instance_type == "castro"
        seed_range = 27:40
    else
        seed_range = 27000:27200
    end

    for seed in seed_range
        nome = "$(file_path)$(seed).txt"
        !isfile(nome) && continue

        print("Solving $nome...")
        if instance_type == "castro"
            V, C, A, d, c, s_node, t_node, k_max = read_castro_instance(nome)
        else
            V, C, A, d, c, s_node, t_node, k_max = read_ferrone_instance(nome)
        end

        status, obj, tempo = risolvi(V, C, A, d, c, s_node, t_node, k_max)

        if status == :optimal
            println("OTTIMO (z = $obj) $(round(tempo, digits=4))s")
            tempo_tot += tempo
            risolte += 1
        elseif status == :time_limit
            println("TIME LIMIT")
        elseif status == :infeasible
            println("INFEASIBLE $(round(tempo, digits=4))s")
        else
            println("ERRORE ($status) $(round(tempo, digits=4))s")
        end
    end

    # if risolte > 0
    #     println("Istanze ottime trovate: $risolte")
    #     println("Tempo totale di calcolo: $(round(tempo_tot, digits=4)) sec")
    #     println("Average Time           : $(round(tempo_tot/risolte, digits=4)) sec")
    # else
    #     println("\nNessuna istanza risolta con successo (controlla i file o i thread).")
    # end
end
inizio_calcolo = time()

batch_execution("ferone-instances/R1/Random_75000x750000_112500_", "ferrone")
batch_execution("ferone-instances/R3/Random_75000x750000_150000_", "ferrone")
batch_execution("ferone-instances/R9/Random_125000x2500000_500000_", "ferrone")
batch_execution("ferone-instances/G1/Grid_100x100_5940_", "ferrone")
batch_execution("ferone-instances/G3/Grid_100x200_11910_", "ferrone")
batch_execution("ferone-instances/G4/Grid_250x500_99700_", "ferrone") 
batch_execution("ferone-instances/G6/Grid_500x1000_39940_", "ferrone")
#batch_execution("castro-instances/", "castro")
#batch_execution("new-instances/", "castro")
tempo_totale = time() - inizio_calcolo
println("Total time: $(round(tempo_totale,   digits=4)) sec")

