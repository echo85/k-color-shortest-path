
set V; # Nodi Grafo
set A within {V,V};
set C; # Colori
param k; # Vincolo massimo colori
param d {A};	# Costi Archi
param s symbolic in V;
param t symbolic in V;
param c {A} symbolic in C; # assegnamo ogni arco ad un colore

var x {A} binary;
var y {C} binary;

minimize z: sum{(u, v) in A} d[u,v] * x[u,v];

subject to Flow_Balance {u in V}:
   sum {(v, u) in A} x[v, u] - sum {(u, v) in A} x[u, v] = 
      (if u == s then -1 else if u == t then 1 else 0);
    
subject MaxColor:
	sum{h in C} y[h] <= k;
	
subject Activation {(u,v) in A}:
	x[u,v] <= y[c[u,v]]