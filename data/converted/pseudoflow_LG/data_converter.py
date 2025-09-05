#!/usr/bin/env python3
"""
Convierte MARVIN_BM.csv en nodes.csv y arcs.csv con talud 1V:1H (45°)
y guarda todo en la carpeta 'output'.
"""

import pandas as pd
from pathlib import Path
import numpy as np
import sys

# ----------------------------------------------------------------------
# 1. Rutas de entrada / salida  ----------------------------------------
SRC  = Path(r"C:/Users/Giova/OneDrive/Escritorio/Mining/data/MARVIN_BM.csv")                 # <-- ajusta si es necesario
OUT  = Path("output")                              # carpeta destino
DST_N = OUT / "nodes.csv"
DST_A = OUT / "arcs.csv"
OUT.mkdir(exist_ok=True)

if not SRC.exists():
    sys.exit(f"No se encuentra el CSV de origen: {SRC}")

# ----------------------------------------------------------------------
# 2. Leer modelo de bloques -------------------------------------------
df = pd.read_csv(SRC)

# ----------------------------------------------------------------------
# 3. Deducir dimensiones de bloque (dx,dy,dz) --------------------------
def cell_size(coord_series):
    diffs = coord_series.sort_values().diff().dropna()
    diffs = diffs[diffs > 0]          # descartar ceros
    if diffs.empty:
        sys.exit("Error: sólo existe un nivel en alguna coordenada; "
                 "no se puede deducir tamaño de bloque.")
    return diffs.min()

dx = cell_size(df["X"])
dy = cell_size(df["Y"])
dz = cell_size(df["Z"])

print(f"Tamaño de bloque detectado: {dx} × {dy} × {dz}")

# ----------------------------------------------------------------------
# 4. Asignar ID y calcular tonelaje ------------------------------------
df = df.sort_values(["Z", "Y", "X"], ascending=[False, True, True]).reset_index(drop=True)
df["id"] = df.index
VOLUMEN = dx * dy * dz
df["tonelaje"] = df["Density"] * VOLUMEN

nodes = df[["id","X","Y","Z","Economic","tonelaje","CU","AU"]]
nodes.columns = ["id","x","y","z","valor","tonelaje","gradeCu","gradeAu"]
nodes.to_csv(DST_N, index=False)

# ----------------------------------------------------------------------
# 5. Generar arcs.csv (talud 45°) --------------------------------------
index_lookup = {(r.X, r.Y, r.Z): int(r.id) for r in df.itertuples()}
xs = sorted(df["X"].unique())
ys = sorted(df["Y"].unique())
zs = sorted(df["Z"].unique(), reverse=True)

arcs = []
for r in df.itertuples():
    x0, y0, z0, id0 = r.X, r.Y, r.Z, r.id
    for z1 in zs:
        if z1 <= z0:
            break
        k = int(round((z1 - z0) / dz))
        max_off = k * dz                       # 45°
        for x1 in xs:
            if abs(x1 - x0) > max_off:
                continue
            for y1 in ys:
                if abs(y1 - y0) > max_off:
                    continue
                id1 = index_lookup.get((x1, y1, z1))
                if id1 is not None:
                    arcs.append((id1, id0))    # id1 encima de id0

pd.DataFrame(arcs, columns=["from", "to"]).to_csv(DST_A, index=False)

print(f"Creado {DST_N}  ({len(nodes)} bloques)")
print(f"Creado {DST_A}  ({len(arcs)} arcos)")
