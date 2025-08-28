#include "fitness.hpp"
#include <algorithm>
#include <numeric>
#include <functional>

// 1) Convierte la posición continua x en permutación de bloques [1..n–1]
std::vector<int> perm_from_continuous(const std::vector<double>& x) {
    int n = static_cast<int>(x.size());
    std::vector<int> idx(n);
    std::iota(idx.begin(), idx.end(), 1);  // [1,2,...,n]
    std::stable_sort(idx.begin(), idx.end(),
        [&](int a, int b){ return x[a-1] > x[b-1]; });
    return idx;
}

// 2) Calcula el NPV total sumando (value * tonnage) de cada bloque extraído
double pit_npv(const std::vector<int>& perm,
               const BlockArray& B,
               const ArcAdj& pred) {
    const int n = static_cast<int>(B.size());
    std::vector<bool> in_pit(n, false);

    // Lambda recursiva: marca v y todos sus predecesores
    std::function<void(int)> mark = [&](int v) {
        if (in_pit[v]) return;
        in_pit[v] = true;
        auto it = pred.find(v);
        if (it != pred.end()) {
            for (int u : it->second) {
                mark(u);
            }
        }
    };

    // Recorre la permutación: extrae solo si value>0
    for (int v : perm) {
        if (v > 0 && v < n && B[v].value > 0.0) {
            mark(v);
        }
    }

    // Suma el flujo neto total: value_por_tonelada * toneladas
    double npv = 0.0;
    for (int v = 1; v < n; ++v) {
        if (in_pit[v]) {
            npv += B[v].value * B[v].tonnage;
        }
    }
    return npv;
}
