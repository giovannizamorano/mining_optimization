# pit_optimization_fixed_full.py

# --- PARÁMETROS ---

import pandas as pd
import math
from collections import defaultdict
import pulp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# --------------------------- CONFIGURACIÓN ---------------------------
BLOCK_SIZE = 30
INPUT_FILE = "C:/Users/Giova/OneDrive/Escritorio/Mining/data/MARVIN_BM.csv"
SLOPE_ANGLE_DEGREES = 45    # Ángulo de talud en grados
# ---------------------------------------------------------------------

# --- 1) CARGAR DATOS ---
df = pd.read_csv(INPUT_FILE)
df["block_id"] = df.index
# Valores económicos ya por bloque (float)

# --- 2) INDICES DE GRILLA (i,j,k) ---
# Facilita generación de precedencias 1x5 (45°)
x_min, y_min, z_min = df["X"].min(), df["Y"].min(), df["Z"].min()
df["i"] = ((df["X"] - x_min) / BLOCK_SIZE).round().astype(int)
df["j"] = ((df["Y"] - y_min) / BLOCK_SIZE).round().astype(int)
df["k"] = ((df["Z"] - z_min) / BLOCK_SIZE).round().astype(int)

# Mapeo de coordenadas discretas a fila de DataFrame
coord_to_idx = { (r.i, r.j, r.k): idx for idx, r in df.iterrows() }

# --- 3) GENERAR PRECEDENCIAS 1x5 ---
#offsets = [ ( 0,  0, -1), ( 1,  0, -1), (-1,  0, -1), ( 0,  1, -1), ( 0, -1, -1) ]
# --- 3) GENERAR PRECEDENCIAS 1x9 ---
offsets = [
    (dx, dy, -1)
    for dx in (-1, 0, 1)
    for dy in (-1, 0, 1)
]

precedences = []  # lista de tuplas (i -> j)
for idx, row in df.iterrows():
    if row.k == 0:
        continue
    for di, dj, dk in offsets:
        nbr = (row.i + di, row.j + dj, row.k + dk)
        if nbr in coord_to_idx:
            precedences.append((idx, coord_to_idx[nbr]))

# --- 4) FORMULAR ILP (Cierre máximo = LG puro) ---
num_blocks = len(df)
prob = pulp.LpProblem("Pit_Optimal_LG", pulp.LpMaximize)
# Variables binarias x_i: 1 si el bloque i entra al pit, 0 si no
x = pulp.LpVariable.dicts('x', range(num_blocks), cat='Binary')

# Objetivo: maximizar suma de valores por bloque
prob += pulp.lpSum(df.loc[i, 'Economic'] * x[i] for i in range(num_blocks)), "Total_Value"

# Restricciones de precedencia: si incluimos i, debemos incluir j
for i, j in precedences:
    prob += x[j] >= x[i], f"prec_{i}_requires_{j}"

# --- 5) RESOLVER ---
solver = pulp.PULP_CBC_CMD(msg=True)
status = prob.solve(solver)
if pulp.LpStatus[prob.status] != 'Optimal':
    raise RuntimeError(f"Solver no encontró óptimo: {pulp.LpStatus[prob.status]}")

# --- 6) EXTRAER Y MOSTRAR RESULTADOS ---
selected = [i for i in range(num_blocks) if pulp.value(x[i]) > 0.5]
benefit  = df.loc[selected, 'Economic'].sum()
print('--- RESULTADOS DEFINITIVOS (LG Puro) ---')
print(f'Bloques en pit óptimo: {len(selected)}')
print(f'Beneficio neto total: {benefit:.2f}\n')
print('Primeros 20 bloques (block_id, X, Y, Z):')
print(df.loc[selected, ['block_id','X','Y','Z']].head(20).to_string(index=False))

# --- 7) VISUALIZACIÓN 3D ---
fig = plt.figure(figsize=(10,8))
ax  = fig.add_subplot(111, projection='3d')
coords = df.loc[selected, ['X','Y','Z']]
ax.scatter(coords['X'], coords['Y'], coords['Z'], s=10, alpha=0.7)
ax.set_title(f'Pit Óptimo 3D - LG ({len(selected)} bloques)')
ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
plt.show()
