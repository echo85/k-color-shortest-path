using DataStructures
using Base.Threads

function read_file(file_path::String)
    lines = readlines(file_path)

    header    = split(lines[1])
    num_nodi  = parse(Int,     header[1])
    k_max     = parse(Int,     header[2])
    s_node    = parse(Int,     header[3])
    t_node    = parse(Int,     header[4])

    gradi = [parse(Int, split(lines[i+1])[1]) for i in 1:num_nodi]

    V   = Set{Int}(1:num_nodi)
    C   = Set{Int}()
    A_s = Set{Tuple{Int,Int}}()
    d   = Dict{Tuple{Int,Int}, Float64}()
    c   = Dict{Tuple{Int,Int}, Int}()

    ptr = num_nodi + 2
    for u in 1:num_nodi
        for _ in 1:gradi[u]
            parts  = split(lines[ptr])
            v      = parse(Int,     parts[1])
            costo  = parse(Float64, parts[2])
            colore = parse(Int,     parts[3])
            push!(C, colore);  push!(A_s, (u, v))
            d[(u,v)] = costo;  c[(u,v)] = colore
            ptr += 1
        end
    end
    return V, C, collect(A_s), d, c, s_node, t_node, k_max
end

# ─────────────────────────────────────────────────────────────
#  STRUTTURA LABEL
# ─────────────────────────────────────────────────────────────
mutable struct Label
    node   :: Int
    dist   :: Float64
    colors :: Set{Int}
    path   :: Vector{Int}   # forward: [s,…,v]  |  backward: [t,…,v]
    valid  :: Bool
end

# ─────────────────────────────────────────────────────────────
#  DOMINANZA  (Definition 2 – invariata)
# ─────────────────────────────────────────────────────────────
function domina(l1::Label, l2::Label)::Bool
    if l1.dist <= l2.dist && issubset(l1.colors, l2.colors)
        return l1.dist < l2.dist ||
               length(l1.colors) < length(l2.colors) ||
               l1.colors == l2.colors          # stesso costo + stessi colori
    end
    return false
end

# ─────────────────────────────────────────────────────────────
#  COMBINAZIONE DI DUE HALF-PATH
#
#  l_fwd.path = [s, a, …, v]      (estremo = nodo di incontro v)
#  l_bwd.path = [t, b, …, v]      (estremo = nodo di incontro v)
#
#  Restituisce il costo totale se:
#    (1) |C_fwd ∪ C_bwd| ≤ k_max
#    (2) i nodi interni dei due half-path sono disgiunti (no cicli)
#  altrimenti Inf.
# ─────────────────────────────────────────────────────────────
function combine_cost(l_fwd::Label, l_bwd::Label, k_max::Int)::Float64
    length(union(l_fwd.colors, l_bwd.colors)) > k_max && return Inf

    # internal del forward: path[1:end-1]  (s … penultimo)
    # internal del backward: path[1:end-1] (t … penultimo)
    # Il nodo di incontro v = path[end] in entrambi: escluso dal check.
    internal_fwd = Set(@view l_fwd.path[1:end-1])
    for v in @view l_bwd.path[1:end-1]
        v in internal_fwd && return Inf
    end

    return l_fwd.dist + l_bwd.dist
end

# ─────────────────────────────────────────────────────────────
#  AGGIUNTA LABEL CON DOMINANZA + TENTATIVO DI COMBINAZIONE
#
#  Tutta la sezione critica è sotto nlock (lock del nodo j).
#  Gerarchia dei lock sempre rispettata: nlock → bc_lock  (no deadlock).
# ─────────────────────────────────────────────────────────────
function try_add_label!(
    new_lbl  :: Label,
    D_same   :: Vector{Label},          # etichette stessa direzione al nodo j
    D_opp    :: Vector{Label},          # etichette direzione opposta al nodo j
    queue    :: PriorityQueue{Label,Float64},
    h_val    :: Float64,                # euristica per questa direzione in j
    k_max    :: Int,
    best_cost:: Ref{Float64},
    bc_lock  :: ReentrantLock,
    nlock    :: ReentrantLock           # lock specifico del nodo j
)
    lock(nlock) do
        # ── 1. Controllo dominanza ──
        dominated = false
        to_remove = Int[]
        for (i, l) in enumerate(D_same)
            !l.valid && continue
            if domina(l, new_lbl)
                dominated = true; break
            elseif domina(new_lbl, l)
                push!(to_remove, i)
            end
        end
        dominated && return

        # ── 2. Invalida etichette dominate e inserisce la nuova ──
        for i in to_remove; D_same[i].valid = false; end
        push!(D_same, new_lbl)
        enqueue!(queue, new_lbl, new_lbl.dist + h_val)

        # ── 3. Tentativo di combinazione con etichette opposte ──
        #       (lock su bc_lock acquisito DENTRO nlock → ordine coerente)
        for l_opp in D_opp
            !l_opp.valid && continue
            c = combine_cost(new_lbl, l_opp, k_max)
            if c < Inf
                lock(bc_lock) do
                    if c < best_cost[]; best_cost[] = c; end
                end
            end
        end
    end
end

# ─────────────────────────────────────────────────────────────
#  WORKER GENERICO  (forward e backward usano la stessa funzione)
#
#  Parametri direzionali:
#    fwd_adj  → adj        (forward)  o  back_adj   (backward)
#    D_self   → D_fwd      (forward)  o  D_bwd      (backward)
#    D_other  → D_bwd      (forward)  o  D_fwd      (backward)
#    h        → h_fwd      (forward)  o  h_bwd      (backward)
#    source   → s_node     (forward)  o  t_node     (backward)
#    target   → t_node     (forward)  o  s_node     (backward)
# ─────────────────────────────────────────────────────────────
function run_search!(
    fwd_adj    :: Dict{Int,Vector{Tuple{Int,Float64,Int}}},
    D_self     :: Dict{Int,Vector{Label}},
    D_other    :: Dict{Int,Vector{Label}},
    node_locks :: Dict{Int,ReentrantLock},
    h          :: Dict{Int,Float64},
    source     :: Int,
    target     :: Int,
    k_max      :: Int,
    best_cost  :: Ref{Float64},
    bc_lock    :: ReentrantLock,
    start_time :: Float64,
    time_limit :: Float64
)::Symbol

    L  = PriorityQueue{Label,Float64}()
    l0 = Label(source, 0.0, Set{Int}(), [source], true)

    lock(node_locks[source]) do
        push!(D_self[source], l0)
    end
    enqueue!(L, l0, h[source])

    while !isempty(L)
        # ── Verifica time limit ──
        time() - start_time > time_limit && return :time_limit

        curr = dequeue!(L)
        !curr.valid && continue

        # ── Pruning A*: f(label) ≥ best già trovato ──
        curr.dist + h[curr.node] >= best_cost[] && continue

        # ── Raggiunto il target direttamente (path completo mono-direzione) ──
        if curr.node == target
            lock(bc_lock) do
                if curr.dist < best_cost[]; best_cost[] = curr.dist; end
            end
            continue
        end

        # ── Espansione ──
        for (j, w, col) in fwd_adj[curr.node]
            j in curr.path && continue          # evita cicli nel path

            new_dist   = curr.dist + w
            new_colors = copy(curr.colors); push!(new_colors, col)
            length(new_colors) > k_max                    && continue
            new_dist + h[j] >= best_cost[]                && continue

            try_add_label!(
                Label(j, new_dist, new_colors, vcat(curr.path, j), true),
                D_self[j], D_other[j],
                L, h[j],
                k_max, best_cost, bc_lock, node_locks[j]
            )
        end
    end

    return :finished
end

function risolvi(file_path::String; time_limit::Float64=600.0)

    V_set, C_set, A, d_cost, c_color, s_node, t_node, k_max =
        read_file(file_path)

    # ── Costruzione grafi di adiacenza ──
    adj      = Dict{Int,Vector{Tuple{Int,Float64,Int}}}(u => [] for u in V_set)
    back_adj = Dict{Int,Vector{Tuple{Int,Float64,Int}}}(u => [] for u in V_set)
    rev_nc   = Dict{Int,Vector{Tuple{Int,Float64}}}(u => [] for u in V_set)

    for (u, v) in A
        w, col = d_cost[(u,v)], c_color[(u,v)]
        push!(adj[u],      (v, w, col))
        push!(back_adj[v], (u, w, col))   # grafo inverso con colori (backward)
        push!(rev_nc[v],   (u, w))         # grafo inverso senza colori (Dijkstra)
    end

    # ── Euristica forward: h_fwd[v] = dist(v → t) ──
    #    Dijkstra da t sul grafo inverso (senza colori)
    h_fwd = Dict{Int,Float64}(u => Inf for u in V_set)
    h_fwd[t_node] = 0.0
    pq = PriorityQueue{Int,Float64}(); enqueue!(pq, t_node, 0.0)
    while !isempty(pq)
        u = dequeue!(pq)
        for (v, w) in rev_nc[u]
            d = h_fwd[u] + w
            if d < h_fwd[v]; h_fwd[v] = d; pq[v] = d; end
        end
    end
    h_fwd[s_node] == Inf && return :infeasible, -1.0, 0.0

    # ── Euristica backward: h_bwd[v] = dist(s → v) ──
    #    Dijkstra da s sul grafo forward (senza colori)
    h_bwd = Dict{Int,Float64}(u => Inf for u in V_set)
    h_bwd[s_node] = 0.0
    pq2 = PriorityQueue{Int,Float64}(); enqueue!(pq2, s_node, 0.0)
    while !isempty(pq2)
        u = dequeue!(pq2)
        for (v, w, _) in adj[u]
            d = h_bwd[u] + w
            if d < h_bwd[v]; h_bwd[v] = d; pq2[v] = d; end
        end
    end

    start_time = time()
    D_fwd      = Dict{Int,Vector{Label}}(u => Label[] for u in V_set)
    D_bwd      = Dict{Int,Vector{Label}}(u => Label[] for u in V_set)
    node_locks = Dict{Int,ReentrantLock}(u => ReentrantLock() for u in V_set)
    best_cost  = Ref(Inf)
    bc_lock    = ReentrantLock()

    t_fwd = Threads.@spawn run_search!(
        adj, D_fwd, D_bwd, node_locks, h_fwd,
        s_node, t_node, k_max, best_cost, bc_lock, start_time, time_limit
    )
    t_bwd = Threads.@spawn run_search!(
        back_adj, D_bwd, D_fwd, node_locks, h_bwd,
        t_node, s_node, k_max, best_cost, bc_lock, start_time, time_limit
    )

    s_fwd = fetch(t_fwd)
    s_bwd = fetch(t_bwd)
    elapsed = time() - start_time

    if s_fwd == :time_limit || s_bwd == :time_limit
        return :time_limit, -1.0, elapsed
    elseif best_cost[] < Inf
        return :optimal, best_cost[], elapsed
    else
        return :infeasible, -1.0, elapsed
    end
end

function batch_execution(file_path::String)
    tempo_tot = 0.0
    risolte   = 0

    for seed in 27000:27200
        nome = "$(file_path)_$seed.txt"
        !isfile(nome) && continue

        print("Risoluzione di $nome in corso... ")
        status, obj, tempo = risolvi(nome, time_limit=120.0)

        if status == :optimal
            println("OTTIMO (z = $obj) | Tempo: $(round(tempo, digits=4)) sec")
            tempo_tot += tempo
            risolte   += 1
        elseif status == :time_limit
            println("TIME LIMIT (> 120s)")
        elseif status == :infeasible
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

#batch_execution("k-color-instances/Big/Testing/Grid/Grid_500x1000_399400/Grid_500x1000_399400")
#batch_execution("k-color-instances/Big/Testing/Grid/Grid_250x250_49800/Grid_250x250_49800")
#batch_execution("k-color-instances/Big/Testing/Grid/Grid_100x100_7920/Grid_100x100_7920")
#batch_execution("k-color-instances/Small/Grid/Grid_250x500_9970/Grid_250x500_9970")
#batch_execution("k-color-instances/Small/Grid/Grid_500x1000_19970/Grid_500x1000_19970")
#batch_execution("k-color-instances/Small/Grid/Grid_500x1000_39940/Grid_500x1000_39940")
batch_execution("k-color-instances/Big/Testing/Random/Random_125000x2500000_500000/Random_125000x2500000_500000")