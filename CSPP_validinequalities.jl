using JuMP
using CPLEX
using DataStructures
include("ReadInstance.jl")
using .ReadInstance

function calcola_raggiungibilita(V_set, A)
    n = maximum(V_set)
    T = zeros(Bool, n, n)
    adj = Dict{Int, Vector{Int}}(u => [] for u in V_set)
    for (u, v) in A
        push!(adj[u], v)
    end

    for start_node in V_set
        queue = Queue{Int}()
        enqueue!(queue, start_node)
        visited = falses(n)
        visited[start_node] = true
        T[start_node, start_node] = true
        
        while !isempty(queue)
            curr = dequeue!(queue)
            for neighbor in adj[curr]
                if !visited[neighbor]
                    visited[neighbor] = true
                    T[start_node, neighbor] = true
                    enqueue!(queue, neighbor)
                end
            end
        end
    end
    return T
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

# MODELLO PRINCIPALE (FFP + Valid Inequalities)
function risolvi(V_set, C_set, A, d_cost, c_color, s_node, t_node, k_max)
    # Calcolo Matrice di Raggiungibilità
    A_r, h_f, h_b, UB = riduzione_grafo(V_set, A, d_cost, c_color, s_node, t_node, k_max)
    T = calcola_raggiungibilita(V_set, A_r)
    
    model = Model(CPLEX.Optimizer)
    set_silent(model)
    if UB < Inf
       set_optimizer_attribute(model, "CPX_PARAM_CUTUP", UB)
    end
    # Variabili
    @variable(model, x[A_r], Bin)         # x_uv = 1 se l'arco è usato
    @variable(model, y[C_set], Bin)     # y_h = 1 se il colore è usato
    
    @objective(model, Min, sum(d_cost[a] * x[a] for a in A_r))
    for u in V_set
        archi_in  = [a for a in A_r if a[2] == u]
        archi_out = [a for a in A_r if a[1] == u]
        
        if u == s_node
            @constraint(model, sum(x[a] for a in archi_out) - sum(x[a] for a in archi_in) == 1)
        elseif u == t_node
            @constraint(model, sum(x[a] for a in archi_out) - sum(x[a] for a in archi_in) == -1)
        else
            @constraint(model, sum(x[a] for a in archi_out) - sum(x[a] for a in archi_in) == 0)
        end
    end

    for a in A_r
        @constraint(model, x[a] <= y[c_color[a]])
    end
    
    # Vincolo: Massimo k_max colori (4)
    @constraint(model, sum(y[c] for c in C_set) <= k_max)


    # AGGIUNTA DELLE VALID INEQUALITIES
    # --- TAGLIO 19 (Cut-Set Randomizzati) ---
    tagli_19_aggiunti = 0
    for _ in 1:60
        S = Set{Int}([s_node])
        for v in V_set
            if v != s_node && v != t_node && rand() > 0.5
                push!(S, v)
            end
        end
        # Verifichiamo la condizione [V\S, S] = vuoto (nessun arco che torna indietro)
        archi_back = [a for a in A_r if a[1] ∉ S && a[2] ∈ S]
        if isempty(archi_back)
            archi_cut = [a for a in A_r if a[1] ∈ S && a[2] ∉ S]
            for h in C_set
                archi_h = [a for a in archi_cut if c_color[a] == h]
                if !isempty(archi_h)
                    @constraint(model, sum(x[a] for a in archi_h) <= y[h])
                    tagli_19_aggiunti += 1
                end
            end
        end
    end

    # --- TAGLIO 18 (Semplificato A Priori) ---
    tagli_18_aggiunti = 0
    for h in C_set
        archi_h = [a for a in A_r if c_color[a] == h]
        for i in 1:length(archi_h)
            for j in (i+1):length(archi_h)
                e1, e2 = archi_h[i], archi_h[j]
                # Se la testa di e1 non può raggiungere la coda di e2 e viceversa
                if !T[e1[2], e2[1]] && !T[e2[2], e1[1]]
                    @constraint(model, x[e1] + x[e2] <= y[h])
                    tagli_18_aggiunti += 1
                end
            end
        end
    end

    # --- TAGLIO 17 (Nodi non raggiungibili limitati) ---
    tagli_17_aggiunti = 0
    for v in V_set
        if v == s_node || v == t_node; continue; end
        
        out_colors = Set(c_color[a] for a in A_r if a[1] == v)
        for h in out_colors
            W_h_v = Tuple{Int,Int}[]
            for u in [a[1] for a in A_r if a[2] == v]
                if !T[v, u] # u in \bar{R}(v)
                    for j in [a[2] for a in A_r if a[1] == u && a[2] != v]
                        if !T[j, v] && c_color[(u,j)] == h
                            push!(W_h_v, (u,j))
                        end
                    end
                end
            end
            
            if !isempty(W_h_v)
                archi_out_v = [a for a in A_r if a[1] == v && c_color[a] == h]
                @constraint(model, sum(x[a] for a in W_h_v) + sum(x[a] for a in archi_out_v) <= y[h])
                tagli_17_aggiunti += 1
            end
        end
    end

    set_time_limit_sec(model, 600.0) 
    optimize!(model)
    
    if termination_status(model) == MOI.OPTIMAL
        cost = objective_value(model)
        
        colori_usati = Int[]
        for c in C_set
            if value(y[c]) > 0.5  # Essendo binaria, se è circa 1 il colore è attivo
                push!(colori_usati, c)
            end
        end
        archi_attivi = [a for a in A_r if value(x[a]) > 0.5]
        
        cammino_nodi = [s_node]
        cammino_archi = Tuple{Int, Int}[]
        corrente = s_node
        
        while corrente != t_node
            prossimo_arco = filter(a -> a[1] == corrente, archi_attivi)
            if isempty(prossimo_arco)
                println("ATTENZIONE: Errore nella ricostruzione del cammino (struttura corrotta).")
                break
            end
            
            arco_scelto = prossimo_arco[1]
            push!(cammino_archi, arco_scelto)
            corrente = arco_scelto[2] # Muoviti al nodo successivo
            push!(cammino_nodi, corrente)
        end
        
        return :optimal, cost, solve_time(model)
    else
        println("Nessuna soluzione ottima trovata. Status: ", termination_status(model))
        return termination_status(model), Inf, solve_time(model)
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

        print("Solving $nome..")
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

end

batch_execution("ferone-instances/R1/Random_75000x750000_112500_", "ferrone")
batch_execution("ferone-instances/R3/Random_75000x750000_150000_", "ferrone")
batch_execution("ferone-instances/R9/Random_125000x2500000_500000_", "ferrone")
batch_execution("ferone-instances/G1/Grid_100x100_5940_", "ferrone")
batch_execution("ferone-instances/G3/Grid_100x200_11910_", "ferrone")
batch_execution("ferone-instances/G4/Grid_250x500_99700_", "ferrone") 
batch_execution("ferone-instances/G6/Grid_500x1000_39940_", "ferrone")
batch_execution("castro-instances/", "castro")
batch_execution("new-instances/", "castro")