using JuMP
using CPLEX
using DataStructures
include("ReadInstance.jl")
using .ReadInstance

mutable struct Label
    node::Int
    dist::Float64
    colors::BitSet
    parent::Union{Nothing, Label}
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

    return A_ridotto, h, dist_s, UB # h è h_f (distanza da t), dist_s è h_b (distanza da s)
end

function risolvi(V, C, A, d, c, s_node, t_node, k_max)
    A_r, h_f, h_b, UB = riduzione_grafo(V, A, d, c, s_node, t_node, k_max)
    
    # Use A_r (reduced arc set) instead of the full A
    in_arcs  = Dict(u => Tuple{Int,Int}[] for u in V)
    out_arcs = Dict(u => Tuple{Int,Int}[] for u in V)
    for arco in A_r
        push!(out_arcs[arco[1]], arco)
        push!(in_arcs[arco[2]],  arco)
    end

    rhs = Dict{Int, Int}(u => 0 for u in V)
    rhs[s_node] = -1
    rhs[t_node] =  1

    modello = Model(CPLEX.Optimizer)
    if UB < Inf
        set_optimizer_attribute(modello, "CPX_PARAM_CUTUP", UB)
    end
    set_attribute(modello, "CPXPARAM_MIP_Strategy_NodeSelect", 1)
    set_silent(modello)

    @variable(modello, x[a in A_r], Bin)   # indexed over A_r
    @variable(modello, y[h in C],   Bin)

    @objective(modello, Min, sum(d[a] * x[a] for a in A_r))

    @constraint(modello, Flow_Balance[u in V],
        sum(x[a] for a in in_arcs[u]) - sum(x[a] for a in out_arcs[u]) == rhs[u]
    )

    @constraint(modello, MaxColor,     sum(y[h] for h in C) <= k_max)
    @constraint(modello, Activation[a in A_r], x[a] <= y[c[a]])  # indexed over A_r

    ottimizza_inizio = time()
    set_time_limit_sec(modello, 600.0) 
    optimize!(modello)
    tempo_totale = time() - ottimizza_inizio

    status  = termination_status(modello)
    obj_val = (status == MOI.OPTIMAL) ? objective_value(modello) : -1.0

    return status, obj_val, tempo_totale
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

        print("Solving $nome.. ")
        if instance_type == "castro"
            V, C, A, d, c, s_node, t_node, k_max = read_castro_instance(nome)
        else
            V, C, A, d, c, s_node, t_node, k_max = read_ferrone_instance(nome)
        end

        status, obj, tempo = risolvi(V, C, A, d, c, s_node, t_node, k_max)

        if status == MOI.OPTIMAL
            println("OTTIMO (z = $obj) | Tempo: $(round(tempo, digits=4)) sec")
            tempo_tot += tempo
            risolte   += 1
        elseif status == MOI.TIME_LIMIT
            println("TIME LIMIT")
        elseif status == MOI.INFEASIBLE
            println("INFEASIBLE | Tempo: $(round(tempo, digits=4)) sec")
        else
            println("ERRORE ($status) | Tempo: $(round(tempo, digits=4)) sec")
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