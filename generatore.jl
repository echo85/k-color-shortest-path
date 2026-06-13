function genera_labirinto_frazionario(W::Int, L::Int, k_max::Int)
    # W = larghezza griglia (densità)
    # L = profondità griglia (deve essere maggiore di k_max * 2 per l'inganno totale)
    
    A = Vector{Tuple{Int,Int}}()
    d_cost = Dict{Tuple{Int,Int}, Int}()
    c_color = Dict{Tuple{Int,Int}, Int}()

    s_node = 1
    t_node = 2
    num_nodes = 2

    colore_bypass = 999 
    nodo_prec = s_node
    
    for i in 1:(L+1)
        num_nodes += 1
        nuovo_nodo = num_nodes
        arco = (nodo_prec, nuovo_nodo)
        push!(A, arco)
        d_cost[arco] = 500  # Costo alto, ma fattibile
        c_color[arco] = colore_bypass
        nodo_prec = nuovo_nodo
    end
    arco_fin = (nodo_prec, t_node)
    push!(A, arco_fin)
    d_cost[arco_fin] = 500
    c_color[arco_fin] = colore_bypass

    nodi_griglia = Dict{Int, Vector{Int}}()
    for layer in 1:L
        nodi_griglia[layer] = Int[]
        for w in 1:W
            num_nodes += 1
            push!(nodi_griglia[layer], num_nodes)
        end
    end

    # usiamo esattamente k_max * 2 colori in totale
    tot_colori_trappola = k_max * 2 
    
    # Crea una "finestra scorrevole" di k_max colori per ogni livello
    function colori_livello(lvl::Int)
        start_idx = (lvl - 1) % tot_colori_trappola
        return [ ((start_idx + i) % tot_colori_trappola) + 1 for i in 0:(k_max-1) ]
    end

    # Archi da S al Layer 1
    col_L1 = colori_livello(1)
    for (i, v) in enumerate(nodi_griglia[1])
        arco = (s_node, v)
        push!(A, arco)
        d_cost[arco] = 10
        c_color[arco] = col_L1[(i % k_max) + 1] # Distribuisce equamente i colori
    end

    # Archi densi tra i layer intermedi (tutti verso tutti)
    for layer in 1:L-1
        col_Lnext = colori_livello(layer + 1)
        for u in nodi_griglia[layer]
            for (i, v) in enumerate(nodi_griglia[layer+1])
                arco = (u, v)
                push!(A, arco)
                d_cost[arco] = 10
                c_color[arco] = col_Lnext[(i % k_max) + 1]
            end
        end
    end

    # Archi dal Layer L a T
    for u in nodi_griglia[L]
        arco = (u, t_node)
        push!(A, arco)
        d_cost[arco] = 10
        c_color[arco] = 1 
    end

    num_colors = max(colore_bypass, tot_colori_trappola)
    num_arcs = length(A)

    return num_nodes, num_arcs, s_node, t_node, k_max, num_colors, A, d_cost, c_color
end

function scrivi_istanza_txt(nome_file::String, n, m, s, t, k, num_colors, A, d_cost, c_color)
    open(nome_file, "w") do f
        # Formato esatto richiesto: n m s t k |C|
        write(f, "$n $m $s $t $k $num_colors\n")
        
        for arco in A
            u, v = arco[1], arco[2]
            cost = d_cost[arco]
            color = c_color[arco]
            write(f, "$u $v $cost $color\n")
        end
    end
    println("$nome_file | Nodi: $n, Archi: $m, k_max: $k, Colori totali trappola: $(k*2)")
end

n, m, s, t, k, c, A, d, col = genera_labirinto_frazionario(70, 35, 30)
scrivi_istanza_txt("new-instances/12.txt", n, m, s, t, k, c, A, d, col)

n, m, s, t, k, c, A, d, col = genera_labirinto_frazionario(80, 45, 40)
scrivi_istanza_txt("new-instances/13.txt", n, m, s, t, k, c, A, d, col)

n, m, s, t, k, c, A, d, col = genera_labirinto_frazionario(90, 55, 50)
scrivi_istanza_txt("new-instances/14.txt", n, m, s, t, k, c, A, d, col)