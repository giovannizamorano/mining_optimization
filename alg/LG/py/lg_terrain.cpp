// lg2_approach2.cpp
// Implementación 100% fiel al Approach 2 (All-Root Connection) del capítulo 5.7.5 de Hustrulid et al.

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <tuple>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>

const int BLOCK = 30;  // Tamaño de celda para discretizar coordenadas

namespace std {
    // Especialización de hash para tuple<int,int,int>
    template <>
    struct hash<tuple<int,int,int>> {
        size_t operator()(tuple<int,int,int> const& t) const noexcept {
            size_t h1 = hash<int>{}(get<0>(t));
            size_t h2 = hash<int>{}(get<1>(t));
            size_t h3 = hash<int>{}(get<2>(t));
            return h1 ^ (h2 << 1) ^ (h3 << 2);
        }
    };
}

struct Block {
    int id;
    int i, j, k;           // Coordenadas discretizadas
    double value;          // Valor económico
    std::vector<int> neigh; // IDs de bloques predecesores (k-1)
};

int main(int argc, char** argv) {
    using namespace std;
    // --- RUTA HARDCODEADA AL CSV ---
    const string inputPath =
        R"(C:\Users\Giova\OneDrive\Escritorio\Mining\data\MARVIN_BM.csv)";
    ifstream file(inputPath);
        if (!file.is_open()) {
        cerr << "Error abriendo archivo " << inputPath << "\n";
        return 1;
    }
    string line;
    getline(file, line); // Se salta la cabecera

    vector<tuple<double,double,double,double>> rawData;
    double x_min = numeric_limits<double>::infinity();
    double y_min = numeric_limits<double>::infinity();
    double z_min = numeric_limits<double>::infinity();

    while (getline(file, line)) {
        if (line.empty()) continue;
        stringstream ss(line);
        string token;
        double x, y, z, value;
        // Se asume: X,Y,Z,<otro>,Value
        getline(ss, token, ','); x = stod(token);
        getline(ss, token, ','); y = stod(token);
        getline(ss, token, ','); z = stod(token);
        getline(ss, token, ',');            // se descarta un campo extra
        if (!getline(ss, token, ',')) continue;
        value = stod(token);

        rawData.emplace_back(x, y, z, value);
        x_min = min(x_min, x);
        y_min = min(y_min, y);
        z_min = min(z_min, z);
    }
    file.close();

    // --- Discretización y asignación de IDs ---
    vector<Block> blocks;
    blocks.reserve(rawData.size());
    unordered_map<tuple<int,int,int>,int> coordToId;
    int nextId = 0;
    for (auto& tup : rawData) {
        double x, y, z, value;
        tie(x, y, z, value) = tup;
        int ii = int(round((x - x_min) / BLOCK));
        int jj = int(round((y - y_min) / BLOCK));
        int kk = int(round((z - z_min) / BLOCK));
        blocks.push_back({ nextId, ii, jj, kk, value, {} });
        coordToId[{ii, jj, kk}] = nextId;
        ++nextId;
    }

    int N = blocks.size();
    // --- Construcción de listas de vecinos (patrón 1–5) ---
    for (auto& b : blocks) {
        int ii = b.i, jj = b.j, kk = b.k;
        vector<tuple<int,int,int>> preds = {
            {ii,   jj,   kk-1},
            {ii-1, jj,   kk-1},
            {ii+1, jj,   kk-1},
            {ii,   jj-1, kk-1},
            {ii,   jj+1, kk-1}
        };
        for (auto& c : preds) {
            auto it = coordToId.find(c);
            if (it != coordToId.end()) {
                b.neigh.push_back(it->second);
            }
        }
    }

    // --- Approach 2: All-Root Connection ---
    int source = N;  // nodo raíz ficticio
    vector<int> parent(N+1, source);
    parent[source] = -1;

    vector<vector<int>> children(N+1);
    auto build_children = [&]() {
        for (auto& ch : children) ch.clear();
        for (int u = 0; u <= N; ++u) {
            int p = parent[u];
            if (p >= 0) children[p].push_back(u);
        }
    };

    vector<double> subtreeW(N+1, 0.0);
    function<double(int)> dfs = [&](int u) {
        double sum = (u < N ? blocks[u].value : 0.0);
        for (int v : children[u]) {
            sum += dfs(v);
        }
        subtreeW[u] = sum;
        return sum;
    };

    // Pasos 2–8: iterar reconectando **un solo** bloque cada vez
    while (true) {
        build_children();
        dfs(source);

        // Paso 2: Y₀ = hijos directos del root con peso positivo
        vector<int> Y;
        for (int v : children[source]) {
            if (subtreeW[v] > 0.0) {
                Y.push_back(v);
            }
        }
        if (Y.empty()) break;

        // Paso 3: elegir el mejor candidato u→v
        double bestW = 0.0;
        int bestU = -1, bestV = -1;
        for (int u : Y) {
            for (int v : blocks[u].neigh) {
                if (parent[v] == source && subtreeW[v] > bestW) {
                    bestW = subtreeW[v];
                    bestU = u;
                    bestV = v;
                }
            }
        }
        // Si no hay candidato con peso > 0, terminamos
        if (bestV < 0) break;

        // Paso 3b: reconectar **solo** bestV bajo bestU
        parent[bestV] = bestU;
        // y volvemos a re-normalizar en la próxima iteración
    }

    // Paso 9: extracción del cierre máximo
    build_children();
    dfs(source);

    vector<int> closure;
    function<void(int)> collect = [&](int u) {
        for (int v : children[u]) {
            if (subtreeW[v] > 0.0) {
                closure.push_back(v);
                collect(v);
            }
        }
    };
    collect(source);

    vector<bool> inPit(N, false);
    for (int id : closure) {
        if (id < N) inPit[id] = true;
    }

    // --- Salida de resultados ---
    double totalProfit = 0.0;
    int totalBlocks = 0;
    ofstream out("pit_blocks.csv");
    out << "X,Y,Z,Value\n";
    for (const auto& b : blocks) {
        if (inPit[b.id]) {
            out << b.i << "," << b.j << "," << b.k << "," << b.value << "\n";
            ++totalBlocks;
            totalProfit += b.value;
        }
    }
    out.close();

    cout << "Bloques en pit: " << totalBlocks << "\n";
    cout << "Beneficio total: " << totalProfit << "\n";
    cout << "Exportado a pit_blocks.csv\n";

    return 0;
}
