#pragma once

#include <vector>
#include "block.hpp"    // define BlockArray y ArcAdj

// Decodifica un vector continuo x ∈ [0,1]ⁿ⁻¹ a una permutación de índices [1..n–1]
// ordenando de mayor a menor valor de x[i].
std::vector<int> perm_from_continuous(const std::vector<double>& x);

// Dada una permutación de extracción y el grafo de precedencias,
// marca el pit (talud 45°) y devuelve el NPV total:
//   ∑ (value_por_tonelada × tonnage) de los bloques incluidos.
double pit_npv(const std::vector<int>& perm,
               const BlockArray& B,
               const ArcAdj& pred);
