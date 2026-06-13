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

    print("[Rid CCDA]", UB == Inf ? "N\\A" : UB)
    print(" Archi:$(length(A))->$(length(A_ridotto)) ")

    return A_ridotto, h, dist_s, UB 
end

function domina(l1::Label, l2::Label)
    if l1.dist <= l2.dist && length(l1.colors) <= length(l2.colors)
        if issubset(l1.colors, l2.colors)
            return true 
        end
    end
    return false
end

function extend(j, d_j, new_colors, l_i, k_max, h, best_cost, H)
    if d_j + h[j] >= best_cost return nothing end
    if d_j > H return nothing end

    if length(new_colors) <= k_max
        l_j = Label(j, d_j, new_colors, l_i)
        return l_j
    end
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
    
    best_cost = UB
    H = UB / 2.0 
    D_f = Dict{Int, Vector{Label}}(u => Label[] for u in V_set)
    D_b = Dict{Int, Vector{Label}}(u => Label[] for u in V_set)
    
    Gamma_bar_f = Dict{Int, Vector{Label}}(u => Label[] for u in V_set)
    Gamma_bar_b = Dict{Int, Vector{Label}}(u => Label[] for u in V_set)
    E_queue = Int[s_node, t_node] 
    
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
                l_j = extend(j, l_i.dist + w_ij, new_colors, l_i, k_max, h_f, best_cost, H)
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
                
                l_j = extend(j, l_i.dist + w_ij, new_colors, l_i, k_max, h_b, best_cost, H)
                add_label(D_b, Gamma_bar_b, l_j)
                if !isempty(Gamma_bar_b) push!(E_queue, j) end
            end
        end
        empty!(Gamma_bar_b[i])
    end
    inizio_join = time()
    # JOIN
    for arco in A_ridotto
        i, j = arco[1], arco[2]
        cost_ij = d_cost[arco]
        col_ij = c_color[arco]
        
        if isempty(D_f[i]) || isempty(D_b[j]) continue end
        Fw = sort(D_f[i], by = x -> x.dist)
        Bw = sort(D_b[j], by = x -> x.dist)
        
        if Fw[1].dist + cost_ij + Bw[1].dist >= best_cost continue end
        
        for l_fw in Fw
            if l_fw.dist + cost_ij + Bw[1].dist >= best_cost break end
            
            for l_bw in Bw
                costo_totale = l_fw.dist + cost_ij + l_bw.dist
                if costo_totale >= best_cost break end # Pruning rapido
                
                # Check Feasibility (Colori)
                merged_colors = union(l_fw.colors, l_bw.colors)
                push!(merged_colors, col_ij)
                
                if length(merged_colors) <= k_max
                    best_cost = costo_totale
                end
            end
        end
    end
    tempo_join = time() - inizio_join
    tempo_totale = time() - inizio_calcolo

    if best_cost[] < Inf
        return :optimal, best_cost, tempo_totale, tempo_join
    else
        return :infeasible, -1.0, tempo_totale, tempo_join
    end
end

function batch_execution(file_path::String, instance_type::String = "ferrone")
    tempo_tot = 0.0
    risolte = 0
    if instance_type == "castro"
        seed_range = 1:40
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

        status, obj, tempo, tempo_join = risolvi(V, C, A, d, c, s_node, t_node, k_max)

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
end

batch_execution("ferone-instances/R1/Random_75000x750000_112500_", "ferrone")
batch_execution("ferone-instances/R3/Random_75000x750000_150000_", "ferrone")
batch_execution("ferone-instances/R9/Random_125000x2500000_500000_", "ferrone")
batch_execution("ferone-instances/G1/Grid_100x100_5940_", "ferrone")
batch_execution("ferone-instances/G3/Grid_100x200_11910_", "ferrone")
batch_execution("ferone-instances/G4/Grid_250x500_99700_", "ferrone") 
batch_execution("ferone-instances/G6/Grid_500x1000_39940_", "ferrone")
batch_execution("new-instances/", "castro")
batch_execution("castro-instances/", "castro")