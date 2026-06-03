using JuMP
using GLPK
using DataStructures
include("ReadInstance.jl")
using .ReadInstance

# FUNZIONE DI SUPPORTO: CALCOLO RAGGIUNGIBILITA' (BFS)
# Restituisce una matrice booleana T dove T[u, v] = true se v è raggiungibile da u
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

# MODELLO PRINCIPALE (FFP + Valid Inequalities)
function risolvi(V_set, C_set, A, d_cost, c_color, s_node, t_node, k_max)
    # Calcolo Matrice di Raggiungibilità
    T = calcola_raggiungibilita(V_set, A)
    
    # Inizializza il modello con GLPK (o CPLEX.Optimizer)
    model = Model(GLPK.Optimizer)
    
    # Variabili
    @variable(model, x[A], Bin)         # x_uv = 1 se l'arco è usato
    @variable(model, y[C_set], Bin)     # y_h = 1 se il colore è usato
    
    # Funzione Obiettivo: Minimizzare il costo
    @objective(model, Min, sum(d_cost[a] * x[a] for a in A))
    
    # Vincoli di Conservazione del Flusso (2)
    for u in V_set
        archi_in  = [a for a in A if a[2] == u]
        archi_out = [a for a in A if a[1] == u]
        
        if u == s_node
            @constraint(model, sum(x[a] for a in archi_out) - sum(x[a] for a in archi_in) == 1)
        elseif u == t_node
            @constraint(model, sum(x[a] for a in archi_out) - sum(x[a] for a in archi_in) == -1)
        else
            @constraint(model, sum(x[a] for a in archi_out) - sum(x[a] for a in archi_in) == 0)
        end
    end
    
    # Vincolo: Se usi l'arco, devi usare il suo colore (3)
    for a in A
        @constraint(model, x[a] <= y[c_color[a]])
    end
    
    # Vincolo: Massimo k_max colori (4)
    @constraint(model, sum(y[c] for c in C_set) <= k_max)


    # AGGIUNTA DELLE VALID INEQUALITIES
    println("Generazione Tagli in corso...")

    # --- TAGLIO 19 (Cut-Set Randomizzati) ---
    # Generiamo 60 partizioni randomizzate per creare tagli S / V\S
    tagli_19_aggiunti = 0
    for _ in 1:60
        S = Set{Int}([s_node])
        for v in V_set
            if v != s_node && v != t_node && rand() > 0.5
                push!(S, v)
            end
        end
        # Verifichiamo la condizione [V\S, S] = vuoto (nessun arco che torna indietro)
        archi_back = [a for a in A if a[1] ∉ S && a[2] ∈ S]
        if isempty(archi_back)
            archi_cut = [a for a in A if a[1] ∈ S && a[2] ∉ S]
            for h in C_set
                archi_h = [a for a in archi_cut if c_color[a] == h]
                if !isempty(archi_h)
                    @constraint(model, sum(x[a] for a in archi_h) <= y[h])
                    tagli_19_aggiunti += 1
                end
            end
        end
    end
    println("Aggiunti $tagli_19_aggiunti tagli di tipo (19)")


    # --- TAGLIO 18 (Semplificato A Priori) ---
    # Impedisce l'attivazione simultanea di archi dello stesso colore non consecutivi 
    # che sono mutuamente irraggiungibili.
    tagli_18_aggiunti = 0
    for h in C_set
        archi_h = [a for a in A if c_color[a] == h]
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
    println("Aggiunti $tagli_18_aggiunti tagli di tipo (18) statici")


    # --- TAGLIO 17 (Nodi non raggiungibili limitati) ---
    tagli_17_aggiunti = 0
    for v in V_set
        if v == s_node || v == t_node; continue; end
        
        # Per ogni colore uscente da v
        out_colors = Set(c_color[a] for a in A if a[1] == v)
        for h in out_colors
            W_h_v = Tuple{Int,Int}[]
            # Troviamo gli archi uscenti dal vicinato in ingresso (u) verso nodi non raggiungibili da v
            for u in [a[1] for a in A if a[2] == v]
                if !T[v, u] # u in \bar{R}(v)
                    for j in [a[2] for a in A if a[1] == u && a[2] != v]
                        if !T[j, v] && c_color[(u,j)] == h
                            push!(W_h_v, (u,j))
                        end
                    end
                end
            end
            
            # (Euristica: prendiamo il primo elemento di W per formare il set W_h(v))
            if !isempty(W_h_v)
                archi_out_v = [a for a in A if a[1] == v && c_color[a] == h]
                @constraint(model, sum(x[a] for a in W_h_v) + sum(x[a] for a in archi_out_v) <= y[h])
                tagli_17_aggiunti += 1
            end
        end
    end
    println("Aggiunti $tagli_17_aggiunti tagli di tipo (17)")


    println("Avvio ottimizzazione...")
    set_time_limit_sec(model, 1800.0) # Time limit di 30 minuti
    optimize!(model)
    
    if termination_status(model) == MOI.OPTIMAL
        cost = objective_value(model)
        println("OTTIMO TROVATO! Costo: $cost")
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

#batch_execution("instances/R1/Random_75000x750000_112500_", "ferrone")
#batch_execution("instances/R3/Random_75000x750000_150000_", "ferrone")
#batch_execution("instances/G1/Grid_100x100_5940_", "ferrone")
#batch_execution("instances/G3/Grid_100x200_11910_", "ferrone")
#batch_execution("castro-instances/", "castro")
batch_execution("new-instances/", "castro")