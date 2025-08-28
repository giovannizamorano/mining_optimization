#!/usr/bin/env python3
import pandas as pd
import pulp

# 1) CARGA DE DATOS
df = pd.read_csv('C:/Users/Giova/OneDrive/Escritorio/Mining/data/MARVIN_BM.csv')
n = len(df)

# 2) GENERACIÓN DE PRECEDENCIAS (patrón 1×5 Talud 45°)
BLOCK = 30
x_min, y_min, z_min = df['X'].min(), df['Y'].min(), df['Z'].min()
df['i'] = ((df['X'] - x_min) / BLOCK).round().astype(int)
df['j'] = ((df['Y'] - y_min) / BLOCK).round().astype(int)
df['k'] = ((df['Z'] - z_min) / BLOCK).round().astype(int)
coord_to_idx = {(r.i, r.j, r.k): idx for idx, r in df.iterrows()}
offs = [(0,0,-1),(1,0,-1),(-1,0,-1),(0,1,-1),(0,-1,-1)]
precedences = []
for idx, r in df.iterrows():
    if r.k == 0: continue
    for di,dj,dk in offs:
        nb = (r.i+di, r.j+dj, r.k+dk)
        if nb in coord_to_idx:
            precedences.append((idx, coord_to_idx[nb]))

# 3) FORMULACIÓN ILP: CIERRE MÁXIMO (Lerchs & Grossmann puro)
prob = pulp.LpProblem('UltimatePitLimit', pulp.LpMaximize)
x = pulp.LpVariable.dicts('x', range(n), cat='Binary')

# Objetivo
prob += pulp.lpSum(df.loc[i, 'Economic'] * x[i] for i in range(n)), 'TotalBenefit'
# Precedencias
for i, j in precedences:
    prob += x[j] >= x[i]

# 4) RESOLVER Y EXTRAER SOLUCIÓN
solver = pulp.PULP_CBC_CMD(msg=True)  # CBC devuelve status 'Optimal' al garantizado
status = prob.solve(solver)
if pulp.LpStatus[status] != 'Optimal':
    raise RuntimeError(f'Solver no halló óptimo: {pulp.LpStatus[status]}')

sel = [i for i in range(n) if pulp.value(x[i]) > 0.5]
benefit = df.loc[sel, 'Economic'].sum()

# 5) RESULTADOS
print('=== PIT ÓPTIMO (ILP – LG Puro) ===')
print(f'Bloques en pit: {len(sel)}')
print(f'Beneficio total: {benefit:.2f}\n')

print('Primeros 20 bloques (id, X, Y, Z):')
print(df.loc[sel, ['X','Y','Z']].head(20).to_string(index=False))

# 6) GUARDAR
df.loc[sel, ['X','Y','Z','Economic']].to_csv('pit_definitivo.csv', index=False)
print('\nPit definitivo guardado en pit_definitivo.csv')
