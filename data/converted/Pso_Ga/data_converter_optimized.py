import csv
import random
import numpy as np
import json
from collections import defaultdict
import sys

# Versión optimizada del convertidor de datos
CSV_BLOCKS = sys.argv[1] if len(sys.argv) > 1 else 'MARVIN_BM.csv'

def load_blocks(path):
    """Lee el CSV y carga los bloques"""
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
                'cu_grade': float(row['CU']),  # cambié el nombre para coincidir con C++
                'au_grade': float(row['AU'])   # cambié el nombre para coincidir con C++
            })
    return blocks

def gen_precedences_optimized(blocks):
    """Genera precedencias de forma más eficiente usando grids espaciales"""
    print(f"Generando precedencias para {len(blocks)} bloques...")
    
    # Organizar bloques por coordenadas z (niveles)
    z_levels = defaultdict(list)
    for block in blocks:
        z_levels[block['z']].append(block)
    
    # Ordenar niveles de z
    sorted_z = sorted(z_levels.keys())
    precedences = {b['id']: [] for b in blocks}
    
    # Para cada nivel, verificar precedencias con niveles superiores
    for i, z_lower in enumerate(sorted_z[:-1]):  # No procesar el último nivel
        blocks_lower = z_levels[z_lower]
        
        # Solo verificar algunos niveles superiores (optimización)
        for j in range(i + 1, min(i + 5, len(sorted_z))):  # Máximo 4 niveles arriba
            z_upper = sorted_z[j]
            blocks_upper = z_levels[z_upper]
            dz = z_upper - z_lower
            
            if dz <= 0:
                continue
                
            # Para cada bloque en el nivel inferior
            for b_lower in blocks_lower:
                # Verificar precedencias con bloques del nivel superior
                for b_upper in blocks_upper:
                    dx = abs(b_upper['x'] - b_lower['x'])
                    dy = abs(b_upper['y'] - b_lower['y'])
                    
                    # Pendiente de 45°: distancia horizontal <= diferencia vertical
                    if (dx**2 + dy**2)**0.5 <= dz:
                        precedences[b_lower['id']].append(b_upper['id'])
    
    # Estadísticas
    total_prec = sum(len(precs) for precs in precedences.values())
    print(f"Total de precedencias generadas: {total_prec}")
    
    return precedences

if __name__ == '__main__':
    print("Cargando bloques...")
    blocks = load_blocks(CSV_BLOCKS)
    print(f"Cargados {len(blocks)} bloques")
    
    print("Generando precedencias...")
    precedences = gen_precedences_optimized(blocks)
    
    # Solo generar datos básicos para el algoritmo genético
    data = {
        'blocks': blocks,
        'precedences': precedences
    }
    
    print("Guardando archivo JSON...")
    with open('data.json', 'w') as f:
        json.dump(data, f, indent=2)  # Menos indentación para archivo más pequeño
    
    print("¡Datos generados y guardados en data.json!")
    print(f"Archivo contiene {len(blocks)} bloques y {sum(len(precs) for precs in precedences.values())} precedencias")
