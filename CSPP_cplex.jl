using JuMP
using CPLEX

include("ReadInstance.jl")
using .ReadInstance

function risolvi(V, C, A, d, c, s_node, t_node, k_max)
    # Popolamento dizionare contenente archi entrata N^+ e e archi uscita N^-
    in_arcs = Dict(u => Tuple{Int,Int}[] for u in V)
    out_arcs = Dict(u => Tuple{Int,Int}[] for u in V)
    for arco in A
        push!(out_arcs[arco[1]], arco)
        push!(in_arcs[arco[2]], arco)
    end
    
    # Vincolo di flusso (Righ Hand Side)
    rhs = Dict{Int, Int}(u => 0 for u in V)
    rhs[s_node] = -1
    rhs[t_node] = 1
    
    modello = Model(CPLEX.Optimizer)
    set_attribute(modello, "CPXPARAM_MIP_Strategy_NodeSelect", 1) # 1 = Depth-First | 2 = Best-Bound
    set_silent(modello) 
    
    @variable(modello, x[a in A], Bin)
    @variable(modello, y[h in C], Bin)
    
    @objective(modello, Min, sum(d[a] * x[a] for a in A))
    
    @constraint(modello, Flow_Balance[u in V],
        sum(x[a] for a in in_arcs[u]) - sum(x[a] for a in out_arcs[u]) == rhs[u]
    )
    
    @constraint(modello, MaxColor, sum(y[h] for h in C) <= k_max)
    @constraint(modello, Activation[a in A], x[a] <= y[c[a]])
    
    ottimizza_inizio = time()
    optimize!(modello)
    tempo_totale = time() - ottimizza_inizio
    
    status = termination_status(modello)
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

        print("Risoluzione di $nome in corso... ")
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

    if risolte > 0
        println("Istanze ottime trovate : $risolte")
        println("Tempo totale di calcolo: $(round(tempo_tot,   digits=4)) sec")
        println("Average Time           : $(round(tempo_tot/risolte, digits=4)) sec")
    else
        println("\nNessuna istanza risolta con successo (controlla i file o i thread).")
    end
end

#batch_execution("ferrone-instances/Small/Grid/Grid_100x100_396/Grid_100x100_396_", seed_range=27000:27009)
batch_execution("castro-instances/", "castro")