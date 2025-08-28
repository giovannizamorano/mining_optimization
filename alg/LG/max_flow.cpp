#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <unordered_map>
#include <tuple>
#include <limits>
#include <set>
#include <cmath> // Para round

// ====================== CONFIGURACIÓN Y ESTRUCTURAS ===========================

const int INF = std::numeric_limits<int>::max();
const int BLOCK = 30; // Tamaño de discretización espacial

// --- Etapa 1: Definición de tipos de arcos según Lerchs-Grossmann ---
enum ArcType { PLUS, MINUS, STRONG };

// --- Etapa 1: Representación de aristas etiquetadas del árbol ---
struct TreeEdge {
    int from, to;
    ArcType type;
    int capacity;
};

// --- Etapa 1: Representación de cada nodo/bloque en el árbol ---
struct TreeNode {
    int id;
    double value;
    std::vector<TreeEdge> outEdges;
    bool reachableFromRoot = false;
};

// ======================== ESTRUCTURA DE GRAFO DE FLUJO =========================

struct Edge {
    int to, rev;
    int capacity;
};

// --- Etapa 4: Clase para grafo de flujo, base del modelo de cierre máximo ---
class Graph {
public:
    Graph(int n) : adj(n) {}

    void addEdge(int u, int v, int capacity) {
        adj[u].push_back({v, static_cast<int>(adj[v].size()), capacity});
        adj[v].push_back({u, static_cast<int>(adj[u].size()) - 1, 0});
    }

    bool bfs(int s, int t, std::vector<int>& level, std::vector<int>& parent, std::vector<int>& edgeIndex) {
        std::fill(level.begin(), level.end(), -1);
        level[s] = 0;
        std::queue<int> q;
        q.push(s);
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (int i = 0; i < adj[u].size(); ++i) {
                const Edge& e = adj[u][i];
                if (e.capacity > 0 && level[e.to] == -1) {
                    level[e.to] = level[u] + 1;
                    parent[e.to] = u;
                    edgeIndex[e.to] = i;
                    if (e.to == t) return true;
                    q.push(e.to);
                }
            }
        }
        return false;
    }

    int maxFlow(int s, int t) {
        int flow = 0;
        std::vector<int> level(adj.size()), parent(adj.size()), edgeIndex(adj.size());
        while (bfs(s, t, level, parent, edgeIndex)) {
            int augFlow = INF;
            for (int v = t; v != s; v = parent[v]) {
                int u = parent[v];
                augFlow = std::min(augFlow, adj[u][edgeIndex[v]].capacity);
            }
            for (int v = t; v != s; v = parent[v]) {
                int u = parent[v];
                Edge& e = adj[u][edgeIndex[v]];
                e.capacity -= augFlow;
                adj[v][e.rev].capacity += augFlow;
            }
            flow += augFlow;
        }
        return flow;
    }

    std::vector<bool> reachableFromSource(int s) {
        std::vector<bool> visited(adj.size(), false);
        std::queue<int> q;
        q.push(s);
        visited[s] = true;
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (const Edge& e : adj[u]) {
                if (e.capacity > 0 && !visited[e.to]) {
                    visited[e.to] = true;
                    q.push(e.to);
                }
            }
        }
        return visited;
    }

private:
    std::vector<std::vector<Edge>> adj;
};

// ======================== BLOQUES Y COORDENADAS DISCRETIZADAS ===================

struct Block {
    int id;
    int i, j, k; // Coordenadas en grilla
    double value; // Valor económico
};

namespace std {
    template<>
    struct hash<std::tuple<int, int, int>> {
        size_t operator()(const std::tuple<int, int, int>& t) const {
            auto h1 = hash<int>{}(std::get<0>(t));
            auto h2 = hash<int>{}(std::get<1>(t));
            auto h3 = hash<int>{}(std::get<2>(t));
            return h1 ^ (h2 << 1) ^ (h3 << 2);
        }
    };
}

// ======================== INICIO FUNCIÓN PRINCIPAL ==============================

int main() {
    std::ifstream file("C:/Users/Giova/OneDrive/Escritorio/Mining/data/MARVIN_BM.csv");
    std::string line;
    std::getline(file, line); // Saltar encabezado

    // --- Etapa 2: Lectura de datos y obtención de mínimos ---
    std::vector<Block> blocks;
    std::unordered_map<std::tuple<int, int, int>, int> coordToId;
    int id = 0;
    double x_min = 1e9, y_min = 1e9, z_min = 1e9;
    std::vector<std::tuple<int, int, int, double>> rawData;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        double x, y, z, value;
        std::getline(ss, token, ','); x = std::stod(token);
        std::getline(ss, token, ','); y = std::stod(token);
        std::getline(ss, token, ','); z = std::stod(token);
        std::getline(ss, token, ','); std::getline(ss, token, ',');
        std::getline(ss, token, ','); value = std::stod(token);

        rawData.emplace_back(x, y, z, value);
        x_min = std::min(x_min, x);
        y_min = std::min(y_min, y);
        z_min = std::min(z_min, z);
    }

    // --- Etapa 3: Discretización espacial y asignación de IDs ---
    for (const auto& [x, y, z, value] : rawData) {
        int i = static_cast<int>(round((x - x_min) / BLOCK));
        int j = static_cast<int>(round((y - y_min) / BLOCK));
        int k = static_cast<int>(round((z - z_min) / BLOCK));
        blocks.push_back({id, i, j, k, value});
        coordToId[{i, j, k}] = id++;
    }

    int n = blocks.size();
    int source = n, sink = n + 1;
    Graph g(n + 2);

    // ====================== INICIO MODELO LERCHS-GROSSMANN ===========================

    // --- Etapa 4: Construcción explícita del árbol y grafo dirigido con etiquetas ---
    std::vector<TreeNode> tree(n);
    for (const auto& b : blocks) {
        tree[b.id] = {b.id, b.value};
        int scaled = static_cast<int>(round(b.value * 100));

        // Arcos desde raíz (source) o hacia sink según signo del valor
        if (b.value >= 0) {
            g.addEdge(source, b.id, scaled);
            tree[b.id].outEdges.push_back({source, b.id, PLUS, scaled});
        } else {
            g.addEdge(b.id, sink, -scaled);
            tree[b.id].outEdges.push_back({b.id, sink, MINUS, -scaled});
        }

        // Arcos de precedencia 1–5 (STRONG)
        std::vector<std::tuple<int, int, int>> preds = {
            {b.i, b.j, b.k - 1}, {b.i - 1, b.j, b.k - 1},
            {b.i + 1, b.j, b.k - 1}, {b.i, b.j - 1, b.k - 1},
            {b.i, b.j + 1, b.k - 1}
        };

        for (const auto& coord : preds) {
            if (coordToId.count(coord)) {
                int predId = coordToId[coord];
                g.addEdge(b.id, predId, INF);
                tree[b.id].outEdges.push_back({b.id, predId, STRONG, INF});
            }
        }
    }

    // --- Etapa 5: Resolución del cierre máximo mediante flujo ---
    int maxflow = g.maxFlow(source, sink);
    std::vector<bool> inPit = g.reachableFromSource(source);

    // ====================== FIN MODELO LERCHS-GROSSMANN ==============================

    // --- Etapa 6: Salida y exportación del pit óptimo ---
    double totalProfit = 0;
    int totalBlocks = 0;
    std::ofstream out("pit_blocks.csv");
    out << "X,Y,Z,Value\n";

    for (const auto& b : blocks) {
        if (inPit[b.id]) {
            totalProfit += b.value;
            totalBlocks++;
            out << b.i << "," << b.j << "," << b.k << "," << b.value << "\n";
        }
    }
    out.close();

    std::cout << "Bloques en pit: " << totalBlocks << std::endl;
    std::cout << "Beneficio total: " << totalProfit << std::endl;
    std::cout << "Exportado a pit_blocks.csv" << std::endl;

    return 0;
}
