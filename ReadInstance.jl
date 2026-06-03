module ReadInstance

export read_ferrone_instance, read_castro_instance

function read_ferrone_instance(file_path::String)
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

function read_castro_instance(file_path::String)
    lines = readlines(file_path)

    # Lettura intestazione file (riga 1)
    header = split(lines[1])
    num_nodi  = parse(Int, header[1])
    num_archi = parse(Int, header[2]) # Novità: il numero di archi è al secondo posto
    s_node    = parse(Int, header[3])
    t_node    = parse(Int, header[4])
    k_max     = parse(Int, header[5]) # k_max si è spostato in quinta posizione
    # header[6] è probabilmente il numero totale di colori, ma non serve salvarlo
    
    # Inizializzazione insiemi e dizionari
    V = Set{Int}(1:num_nodi)
    C = Set{Int}()
    A_set = Set{Tuple{Int, Int}}()
    d = Dict{Tuple{Int, Int}, Float64}()
    c = Dict{Tuple{Int, Int}, Int}()

    # Lettura degli archi (partendo dalla riga 2)
    # Cicliamo esattamente per il numero di archi letto dall'intestazione
    for i in 2:(num_archi + 1)
        parts = split(lines[i])
        
        u      = parse(Int, parts[1])
        v      = parse(Int, parts[2])
        costo  = parse(Float64, parts[3])
        colore = parse(Int, parts[4])
        
        push!(C, colore)
        push!(A_set, (u, v))
        d[(u, v)] = costo
        c[(u, v)] = colore
    end
    
    A = collect(A_set)
    
    return V, C, A, d, c, s_node, t_node, k_max
end

end