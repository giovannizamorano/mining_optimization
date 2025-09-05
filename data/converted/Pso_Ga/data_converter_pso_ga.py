import csv
import random
import numpy as np
import json

# Este script genera la representación de soluciones (listas de individuos/partículas) en formato JSON
# Input: MARVIN_BM.csv con columnas id,x,y,z,value,tonnage
# Output: JSON con campos:
#   - blocks: lista de diccionarios con datos de cada bloque
#   - precedences: dict id → lista de ids predecesores para pendiente de 45°
#   - ga_population: lista de soluciones GA (permutaciones)
#   - pso_population: lista de soluciones PSO (posición continua + velocidad)

import sys
CSV_BLOCKS = sys.argv[1] if len(sys.argv) > 1 else 'MARVIN_BM.csv'
GA_SIZE = 500   # tamaño de población GA
PSO_SIZE = 500  # tamaño de población PSO

# Carga bloques con coordenadas y datos
def load_blocks(path):
    """
    Lee MARVIN_BM.csv con columnas:
      X, Y, Z, Rock Type, Density, Economic, CU, AU, Perfin
    Asigna 'id' como el índice de fila (1-based),
    mapea 'Economic' → value, 'Density' → tonnage,
    e incluye los grades CU y AU si los necesitas.
    """
    blocks = []
    with open(path, newline='') as f:
        reader = csv.DictReader(f)
        for idx, row in enumerate(reader, start=1):
            blocks.append({
                'id': idx,
                'x': float(row['X']),
                'y': float(row['Y']),
                'z': float(row['Z']),
                'value': float(row['Economic']),
                'tonnage': float(row['Density']),
                # si quieres usar calidades:
                'cu': float(row['CU']),
                'au': float(row['AU'])
            })
    return blocks

# Genera precedencias para un ángulo de talud de 45°
def gen_precedences(blocks):
    precedences = {b['id']: [] for b in blocks}
    for bi in blocks:
        for bj in blocks:
            dz = bj['z'] - bi['z']
            if dz > 0:
                dx = abs(bj['x'] - bi['x'])
                dy = abs(bj['y'] - bi['y'])
                # pendiente de 45°: distancia horizontal <= diferencia vertical
                if (dx**2 + dy**2)**0.5 <= dz:
                    precedences[bi['id']].append(bj['id'])
    return precedences

# Genera población inicial para GA: permutaciones de IDs
def gen_ga_population(blocks, pop_size):
    ids = [b['id'] for b in blocks]
    return [{'perm': random.sample(ids, len(ids))} for _ in range(pop_size)]

# Genera población inicial para PSO: posiciones y velocidades
def gen_pso_population(blocks, pop_size):
    dim = len(blocks)
    pop = []
    for _ in range(pop_size):
        position = np.random.rand(dim).tolist()
        velocity = [0.0] * dim
        pop.append({'position': position, 'velocity': velocity})
    return pop

if __name__ == '__main__':
    blocks = load_blocks(CSV_BLOCKS)
    precedences = gen_precedences(blocks)
    ga_population = gen_ga_population(blocks, GA_SIZE)
    pso_population = gen_pso_population(blocks, PSO_SIZE)

    # Empaquetar todo en un dict y volcar a JSON
    data = {
        'blocks': blocks,
        'precedences': precedences,
        'ga_population': ga_population,
        'pso_population': pso_population
    }
    with open('data.json', 'w') as f:
        json.dump(data, f, indent=4)
    print("Datos generados y guardados en data.json")
