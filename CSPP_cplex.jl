using JuMP
using CPLEX

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

function risolvi_cplex(file_path::String)
    
    V, C, A, d, c, s_node, t_node, k_max = read_file(file_path)
    
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

function esegui_esperimento(file_path::String)

    tempo_totale_gruppo = 0.0
    istanze_risolte = 0
    
    for seed in 27000:27009
        nome_file = "$(file_path)_$seed.txt"
        if !isfile(nome_file)
            continue
        end
        
        print("Risoluzione di $nome_file in corso... ")
        
        status, obj, tempo = risolvi_cplex(nome_file)
        
        if status == MOI.OPTIMAL
            println("OTTIMO (z = $obj) | Tempo: $(round(tempo, digits=4)) sec")
            tempo_totale_gruppo += tempo
            istanze_risolte += 1
        elseif status == MOI.INFEASIBLE
            println("INFEASIBLE | Tempo: $(round(tempo, digits=4)) sec")
        else
            println("ERRORE ($status) | Tempo: $(round(tempo, digits=4)) sec")
        end
    end
    
    if istanze_risolte > 0
        tempo_medio = tempo_totale_gruppo / istanze_risolte
        println("\n==================================================")
        println("              RISULTATI FINALI                     ")
        println("==================================================")
        println("Istanze ottime trovate : $istanze_risolte")
        println("Tempo totale di calcolo: ", round(tempo_totale_gruppo, digits=4), " sec")
        println("TEMPO MEDIO: ", round(tempo_medio, digits=4), " sec")
        println("==================================================")
    else
        println("\nNessuna istanza risolta con successo (controlla i file).")
    end
end

esegui_esperimento("k-color-instances/Small/Grid/Grid_100x100_396/Grid_100x100_396")