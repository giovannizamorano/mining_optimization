#!/usr/bin/env python3
"""
Lerchs–Grossmann (LG) – Optimized & Visual 3‑D
=============================================
* Precendence graph built with exact slope check (KD‑tree).
* Topological strong‑block extraction → zero violations guaranteed.
* Fast convergence (`NO_PROGRESS_LIMIT = 10_000`).
* Reports Gross Benefit, NPV, Strip Ratio.
* Dual visualization:
  1. 3‑D scatter of pit blocks (Plotly Express).
  2. Smoothed Digital Terrain Model (DTM) via `scipy.interpolate.griddata`.

Dependencies: numpy, scipy, pandas, plotly
"""
import math
import time
from collections import defaultdict, deque, Counter
from typing import Dict, List, Set, Tuple

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.interpolate import griddata
import plotly.express as px
import plotly.graph_objects as go

# ----- CONFIGURATION ---------------------------------------------------------
CSV_PATH = "C:/Users/Giova/OneDrive/Escritorio/Mining/data/MARVIN_BM.csv"
VERBOSE = True
PRINT_EVERY = 50
MAX_VISITS_PER_NODE = 15
NO_PROGRESS_LIMIT = 10_000
BLOCK_SIZE = 30          # meters
SLOPE_ANGLE_DEGREES = 45
DISCOUNT_RATE = 0.08     # 8%
ROOT = -1

# -----------------------------------------------------------------------------
# Tree structure for LG merge/prune
# -----------------------------------------------------------------------------
class Tree:
    def __init__(self, values: Dict[int, float]):
        self.values = values.copy()
        self.parent = {n: ROOT for n in values}
        self.children = defaultdict(set)
        for n in values:
            self.children[ROOT].add(n)
        self.weight = values.copy()
        self._update_all_weights()

    def _update_all_weights(self):
        for n in self.values:
            self.weight[n] = self.values[n]
        for n in sorted(self.values, key=lambda k: self.parent[k], reverse=True):
            p = self.parent[n]
            if p != ROOT:
                self.weight[p] += self.weight[n]

    def _update_up(self, node: int):
        while node != ROOT:
            new_w = self.values[node] + sum(self.weight[c] for c in self.children[node])
            if new_w == self.weight[node]:
                break
            self.weight[node] = new_w
            node = self.parent[node]

    def strong(self, node: int) -> bool:
        return self.weight[node] >= 0

    def merge(self, s: int, w: int) -> bool:
        prev_w = self.weight[s]
        old_p = self.parent[w]
        if old_p == s:
            return False
        self.children[old_p].discard(w)
        self.parent[w] = s
        self.children[s].add(w)
        self._update_up(w)
        if VERBOSE and self.weight[s] < prev_w:
            print(f"[⚠] Merge reduced weight of {s}: {prev_w:.3f} → {self.weight[s]:.3f}")
        return self.weight[s] != prev_w

    def prune(self, root: int):
        stack = deque([root])
        while stack:
            u = stack.pop()
            for v in list(self.children[u]):
                if self.weight[v] >= 0:
                    self.children[u].discard(v)
                    self.parent[v] = ROOT
                    self.children[ROOT].add(v)
                    self._update_up(u)
                else:
                    stack.append(v)

    def strong_blocks(self, edges: Dict[int, List[int]]) -> List[int]:
        indeg = {n: len(edges[n]) for n in self.values}
        depend = defaultdict(list)
        for n, preds in edges.items():
            for p in preds:
                depend[p].append(n)
        q = deque([n for n in self.values if indeg[n] == 0 and self.strong(n)])
        pit = []
        while q:
            n = q.popleft()
            pit.append(n)
            for m in depend[n]:
                indeg[m] -= 1
                if indeg[m] == 0 and self.strong(m):
                    q.append(m)
        return pit

# -----------------------------------------------------------------------------
# Build precedence graph
# -----------------------------------------------------------------------------
def build_ip_edges(df: pd.DataFrame) -> Dict[int, List[int]]:
    coords = df[['X','Y','Z']].values
    slope_rad = math.radians(SLOPE_ANGLE_DEGREES)
    max_xy = math.tan(slope_rad)*BLOCK_SIZE
    tree = cKDTree(coords)
    edges = defaultdict(list)
    for i,(x,y,z) in enumerate(coords):
        idxs = tree.query_ball_point([x,y,z],r=max_xy)
        for j in idxs:
            x2,y2,z2 = coords[j]
            if z2 < z:
                dz = z-z2
                dxy = math.hypot(x-x2,y-y2)
                slope=math.atan2(dz,dxy) if dxy else math.pi/2
                if slope<=slope_rad:
                    edges[i].append(j)
    return edges

# -----------------------------------------------------------------------------
# Utility metrics
# -----------------------------------------------------------------------------
def prec_viol(pit:Set[int],edges:Dict[int,List[int]])->List[Tuple[int,int]]:
    return [(u,v) for u in pit for v in edges[u] if v not in pit]

def calc_npv(df_pit:pd.DataFrame)->float:
    z0=df_pit.Z.min()
    df2=df_pit.copy()
    df2['t']=(df2.Z-z0)/BLOCK_SIZE
    df2['disc']=1/(1+DISCOUNT_RATE)**df2.t
    return float((df2.Economic*df2.disc).sum())

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def run():
    t0=time.time()
    df=pd.read_csv(CSV_PATH).reset_index().rename(columns={'index':'block_id'})
    # volume & tonnage
    dx,dy,dz=[df[c].sort_values().diff().replace(0,np.nan).dropna().mode()[0] for c in ['X','Y','Z']]
    vol=dx*dy*dz
    df['Ton']=df.Density.replace(0,np.nan).fillna(1.0)*vol

    edges=build_ip_edges(df)
    print(f"Precedences: {sum(len(v) for v in edges.values()):,} arcs, avg {np.mean([len(v) for v in edges.values()]):.2f}")

    tree=Tree(dict(zip(df.block_id,df.Economic)))
    q=deque(tree.values.keys()); visits=Counter(); merges=iters=0; stagn=0
    print("Running LG…")
    while q:
        s=q.popleft(); iters+=1; visits[s]+=1
        if visits[s]>MAX_VISITS_PER_NODE or not tree.strong(s): continue
        prog=False
        for w in edges[s]:
            if tree.strong(w): continue
            if tree.merge(s,w): tree.prune(s); q.append(s); merges+=1; prog=True; break
        stagn=0 if prog else stagn+1
        if stagn>=NO_PROGRESS_LIMIT: print(f"Stagnant {NO_PROGRESS_LIMIT} iters, stop"); break

    pit=set(tree.strong_blocks(edges))
    # fix violations
    added=0
    while True:
        viol=prec_viol(pit,edges)
        if not viol: break
        for _,v in viol: pit.add(v)
        added+=len(viol)
    # metrics
    df_pit=df[df.block_id.isin(pit)]
    ore=df_pit[df_pit.Economic>=0]; waste=df_pit[df_pit.Economic<0]
    raw= df_pit.Economic.sum(); npv=calc_npv(df_pit)
    strip = waste.Ton.sum()/ore.Ton.sum() if not ore.empty else float('inf')
    print(f"Iters {iters}, merges {merges}, viol remain {len(prec_viol(pit,edges))}")
    print(f"Blocks {len(pit):,}, benefit {raw:,.0f}, NPV {npv:,.0f}, strip {strip:.2f}, added {added}")
    print(f"Time: {time.time()-t0:.2f}s")

    # visuals
    df_vis = df_pit.copy()
    fig1=px.scatter_3d(df_vis, x='X', y='Y', z='Z', color='Economic', title='Pit 3D Scatter')
    fig1.show()
        # 3-D surface showing pit depth on a flat plane
    # Prepare grid of block centers
    Xi, Yi = np.meshgrid(sorted(df_vis.X.unique()), sorted(df_vis.Y.unique()))
    # Get depth values: relative depth from max surface
    Z_grid = griddata(df_vis[['X','Y']].values, df_vis.Z.values,
                      (Xi, Yi), method='nearest')
    maxZ = df_vis.Z.max()
    depth = maxZ - Z_grid  # deeper = larger depth
    # Create flat surface at elevation just above maxZ
    Z_surf = np.full_like(depth, maxZ + BLOCK_SIZE * 0.1)
    fig2 = go.Figure(data=[
        go.Surface(
            x=Xi, y=Yi, z=Z_surf,
            surfacecolor=depth,
            colorscale='Viridis',
            colorbar=dict(title='Depth'),
            showscale=True
        )
    ])
    fig2.update_layout(
        title='Pit Depth Contour on Flat Plane',
        scene=dict(
            xaxis_title='X', yaxis_title='Y', zaxis_title='Elevation'
        )
    )
    fig2.show()

if __name__=='__main__': run()
