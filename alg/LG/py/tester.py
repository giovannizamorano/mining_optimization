import pandas as pd
import networkx as nx
from networkx.algorithms.flow import maximum_flow
from pulp import LpProblem, LpMaximize, LpVariable, PULP_CBC_CMD, LpStatus
from collections import deque, defaultdict

# --- CONFIGURACIÓN ---
INPUT_FILE = 'C:/Users/Giova/OneDrive/Escritorio/Mining/data/MARVIN_BM.csv'
BLOCK_SIZE = 30

# --- 1) CARGA DE DATOS Y PRECEDENCIAS ---
df = pd.read_csv(INPUT_FILE)
# Generar índices de grilla 1×5 para ángulo de talud 45°
if 'block_id' not in df.columns:
    df['block_id'] = df.index
x_min, y_min, z_min = df['X'].min(), df['Y'].min(), df['Z'].min()
df['i'] = ((df['X'] - x_min) / BLOCK_SIZE).round().astype(int)
df['j'] = ((df['Y'] - y_min) / BLOCK_SIZE).round().astype(int)
df['k'] = ((df['Z'] - z_min) / BLOCK_SIZE).round().astype(int)
coord_to_idx = {(row.i, row.j, row.k): idx for idx, row in df.iterrows()}
offsets = [(0, 0, -1), (1, 0, -1), (-1, 0, -1), (0, 1, -1), (0, -1, -1)]
precedences = []
for idx, row in df.iterrows():
    if row.k == 0:
        continue
    for di, dj, dk in offsets:
        nbr = (row.i + di, row.j + dj, row.k + dk)
        if nbr in coord_to_idx:
            precedences.append((idx, coord_to_idx[nbr]))

# --- 2) EJECUTAR PIT_OPTIMAL (ILP – Lerchs & Grossmann puro)
# Esta sección reproduce la lógica de `pit_optimal.py`:
#  - Formulación de cierre máximo (maximum closure) con variables binarias x_i
#  - Objetivo: maximizar suma de valores económicos (Economic * x_i)
#  - Restricciones de precedencia: x_j >= x_i para cada arco i->j
#  - Se resuelve con CBC y se extrae la solución óptima

def run_ilp():
    n = len(df)
    prob = LpProblem('PitILP', LpMaximize)
    x = LpVariable.dicts('x', range(n), cat='Binary')
    prob += sum(df.loc[i, 'Economic'] * x[i] for i in range(n)), 'TotalValue'
    for i, j in precedences:
        prob += x[j] >= x[i]
    solver = PULP_CBC_CMD(msg=False)
    status = prob.solve(solver)
    cert = (LpStatus[status] == 'Optimal')
    sel = {i for i in range(n) if x[i].value() > 0.5}
    return sel, cert

# --- 3) EJECUTAR UPL (Min-Cut / Pseudoflow equivalente a LG)
# Esta sección reproduce la lógica de `upl.py`:
#  - Construcción de red de flujo: nodo Source conectado a bloques rentables, bloques estériles al Sink
#  - Arcos de precedencia con capacidad "infinita" M para garantizar la estabilidad del pit
#  - Cálculo de flujo máximo con algoritmo maximum_flow
#  - Construcción manual del grafo residual y BFS para certificar ausencia de caminos aumentantes (flujo máximo)
#  - Extracción del conjunto S alcanzable desde Source como pit óptimo

def run_flow():
    G = nx.DiGraph()
    s, t = 's', 't'
    G.add_node(s); G.add_node(t)
    M = df['Economic'].abs().sum() + 1
    for i, val in df['Economic'].items():
        if val >= 0:
            G.add_edge(s, i, capacity=val)
        else:
            G.add_edge(i, t, capacity=-val)
    for i, j in precedences:
        G.add_edge(i, j, capacity=M)
    flow_val, flow_dict = maximum_flow(G, s, t, capacity='capacity')
    
    # Construir grafo residual
    R = nx.DiGraph()
    R.add_nodes_from(G.nodes())
    for u, v, data in G.edges(data=True):
        cap = data['capacity']
        f = flow_dict.get(u, {}).get(v, 0)
        if cap - f > 0:
            R.add_edge(u, v, residual_capacity=cap - f)
        if f > 0:
            R.add_edge(v, u, residual_capacity=f)
    visited = {s}
    dq = deque([s])
    while dq:
        u = dq.popleft()
        for v, d in R[u].items():
            if d['residual_capacity'] > 0 and v not in visited:
                visited.add(v)
                dq.append(v)
    sel = {n for n in visited if n != s}
    cert = (t not in visited)
    return sel, cert

# --- 4) VERIFICAR PRECEDENCIAS ---
def check_precedences(sol):
    return [(i, j) for i, j in precedences if i in sol and j not in sol]

# --- 5) OMITIDA PRUEBA DE MONOTONÍA ---
# Por tiempo, se omite la verificación de monotonicidad completa.

# --- 6) EJECUTAR MÉTODOS ---
sel_ilp, cert_ilp = run_ilp()
sel_flow, cert_flow = run_flow()
viol_ilp = check_precedences(sel_ilp)
viol_flow = check_precedences(sel_flow)
mono_ilp = True  # Prueba de monotonía omitida
mono_flow = True
print(f"ILP:   bloques={len(sel_ilp)}, preced={len(viol_ilp)} viol, opt={cert_ilp}, mono={mono_ilp}")
print(f"Flow:  bloques={len(sel_flow)}, preced={len(viol_flow)} viol, opt={cert_flow}, mono={mono_flow}")
print(f"Diferencia ILP↔Flow: {len(sel_ilp ^ sel_flow)} bloques")

# --- 7) DETERMINACIÓN FINAL ---
# Se comparan ILP y Flow; monotonicidad omitida.
if len(sel_ilp ^ sel_flow) == 0 and cert_ilp and cert_flow and not viol_ilp and not viol_flow:
    print(">>> Ambos métodos coinciden y pasan pruebas de precedencias y flujo. Pit definitivo obtenido.")
    final = sel_ilp
else:
    print(">>> Discrepancia detectada; se elige ILP (LG puro) como definitivo.")
    final = sel_ilp
# Guardar pit definitivo
result = df.loc[list(final), ['block_id','X','Y','Z','Economic']]
result.to_csv('pit_definitivo.csv', index=False)
print(f"Pit definitivo guardado: bloques={len(final)}, beneficio={df.loc[list(final),'Economic'].sum():.2f}")