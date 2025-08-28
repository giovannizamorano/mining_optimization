/*  pseudoflow.cpp  — versión íntegra con
 *  · Capacidades enteras (Cap = long long) y escala SCALE = 1000
 *  · “Pit” = lado T (no alcanzable desde s)  →  respeta teoría de máxima‑clausura
 *  · Controles de depuración:
 *        – ids de arcos sin bloque    (missingU / missingV)
 *        – valores de INF y suma de pesos positivos
 *        – precedencias violadas en el pit resultante
 */

 #include <algorithm>
 #include <chrono>
 #include <cmath>
 #include <fstream>
 #include <iomanip>
 #include <iostream>
 #include <limits>
 #include <queue>
 #include <sstream>
 #include <string>
 #include <unordered_map>
 #include <vector>
 
 using namespace std;
 
 /*-----------------------------  Tipos y constantes  --------------------------*/
 using Cap = long long;
 constexpr int SCALE = 1000;              // 3 decimales → 1 u.m. = 0.001 “pesos”
 
 struct Edge { int to, rev; Cap cap; };
 struct Node { vector<Edge> adj; int label = 0; Cap excess = 0; };
 
 /*------------------------------  Clase HPF  ----------------------------------*/
 class HPF {
     vector<Node> G;
     int n, s, t, maxLabel = 0, activeCnt = 0;
     vector<vector<int>> bucket;
 
     void addEdgeInt(int u, int v, Cap c) {
         G[u].adj.push_back({v, (int)G[v].adj.size(), c});
         G[v].adj.push_back({u, (int)G[u].adj.size() - 1, 0});
     }
 
     void addActive(int v) {
         if (v == s || v == t || G[v].excess <= 0) return;
         bucket[G[v].label].push_back(v);
         ++activeCnt;
         maxLabel = max(maxLabel, G[v].label);
     }
 
     void push(int u, Edge &e) {
         Cap df = min(G[u].excess, e.cap);
         e.cap -= df;
         G[e.to].adj[e.rev].cap += df;
         G[u].excess -= df;
         G[e.to].excess += df;
         addActive(e.to);
     }
 
     void relabel(int u) {
         int d = numeric_limits<int>::max();
         for (auto &e : G[u].adj)
             if (e.cap > 0) d = min(d, G[e.to].label);
 
         G[u].label = (d == numeric_limits<int>::max()) ?
                      2 * n - 1 : min(d + 1, 2 * n - 1);
         addActive(u);
     }
 
 public:
     HPF(int N, int S, int T) : G(N), n(N), s(S), t(T), bucket(2 * N) {}
 
     void addEdge(int u, int v, Cap c) { addEdgeInt(u, v, c); }
 
     /** Devuelve pit (1 = bloque en el pit) = lado T del corte mínimo */
     vector<char> run() {
         G[s].label = n;          // pre‑flujo
         for (auto &e : G[s].adj) {
             if (e.cap > 0) {
                 Cap df = e.cap;
                 e.cap = 0;
                 G[e.to].adj[e.rev].cap += df;
                 G[e.to].excess += df;
                 addActive(e.to);
             }
         }
 
         while (activeCnt) {
             while (maxLabel >= 0 && bucket[maxLabel].empty()) --maxLabel;
             int u = bucket[maxLabel].back();
             bucket[maxLabel].pop_back();
             --activeCnt;
 
             while (G[u].excess > 0) {
                 bool pushed = false;
                 for (auto &e : G[u].adj) {
                     if (e.cap > 0 && G[u].label == G[e.to].label + 1) {
                         push(u, e);
                         pushed = true;
                         if (G[u].excess == 0) break;
                     }
                 }
                 if (!pushed) { relabel(u); break; }
             }
         }
 
         /* BFS: nodos alcanzables desde s  →  lado S */
         vector<char> side(n + 2, 0);
         queue<int> q; side[s] = 1; q.push(s);
         while (!q.empty()) {
             int u = q.front(); q.pop();
             for (auto &e : G[u].adj)
                 if (e.cap > 0 && !side[e.to]) {
                     side[e.to] = 1; q.push(e.to);
                 }
         }
 
         /* PIT = lado S (alcanzable desde s) */
         vector<char> pit(n);
         for (int i = 0; i < n; ++i) pit[i] = !side[i];
         return pit;
     }
 };
 
 /*-----------------------------------  main  ----------------------------------*/
 int main() {
     ios::sync_with_stdio(false); cin.tie(nullptr);
 
     const string BASE = "C:/Users/Giova/OneDrive/Escritorio/Mining/output/";
     const string NODES = BASE + "nodes.csv";
     const string ARCS  = BASE + "arcs.csv";
     const string OUT_DIR = BASE + "pseudoflujo/";
 
     ifstream fn(NODES), fa(ARCS);
     if (!fn || !fa) { cerr << "No se pueden abrir los archivos.\n"; return 1; }
 
     /* Leer bloques */
     string line; getline(fn, line);             // salta cabecera
     struct Block { int id; Cap valInt; double ton, gCu, gAu; };
     vector<Block> blocks;
     unordered_map<int,int> idx;
 
     while (getline(fn, line)) {
         string s; stringstream ss(line); Block b;
         getline(ss, s, ','); b.id = stoi(s);
         for (int i = 0; i < 3; ++i) getline(ss, s, ',');      // X, Y, Z
         getline(ss, s, ','); b.valInt = (Cap) llround(stod(s) * SCALE);
         getline(ss, s, ','); b.ton = stod(s);
         getline(ss, s, ','); b.gCu = stod(s);
         getline(ss, s, ','); b.gAu = stod(s);
         idx[b.id] = (int)blocks.size();
         blocks.push_back(b);
     }
     fn.close();
 
     /* Construir red */
     int n = blocks.size(), S = n, T = n + 1;
     HPF pf(n + 2, S, T);
 
     for (int i = 0; i < n; ++i) {
         if (blocks[i].valInt > 0) pf.addEdge(S, i, blocks[i].valInt);
         else if (blocks[i].valInt < 0) pf.addEdge(i, T, -blocks[i].valInt);
     }
 
     /* Capacidad infinita */
     long double big = 0;
     for (auto &b : blocks) if (b.valInt > 0) big += b.valInt;
     Cap INF = static_cast<Cap>(big) + 1;
     cerr << "Suma de pesos positivos = " << (INF - 1)
          << "  |  INF usado = "           << INF << '\n';
 
     /* Leer aristas y contar ids perdidos */
     size_t missingU = 0, missingV = 0;
     getline(fa, line);                          // salta cabecera
     while (getline(fa, line)) {
         if (line.empty()) continue;
         auto p = line.find(',');
         int idU = stoi(line.substr(0, p));
         int idV = stoi(line.substr(p + 1));
         if (!idx.count(idU)) { ++missingU; continue; }
         if (!idx.count(idV)) { ++missingV; continue; }
         int u = idx[idU];
         int v = idx[idV];
         pf.addEdge(u, v, INF);
     }
     fa.close();
     cerr << "Arcos con u sin nodo: " << missingU
          << " | v sin nodo: "        << missingV << '\n';
 
     /* Ejecutar */
     auto t0 = chrono::steady_clock::now();
     vector<char> pit = pf.run();
     double secs = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
 
     /* Calcular métricas */
     Cap benefitInt = 0; double oreT = 0, wasteT = 0, cuT = 0, auT = 0;
     for (int i = 0; i < n; ++i) if (pit[i]) {
         benefitInt += blocks[i].valInt;
         if (blocks[i].valInt >= 0) {
             oreT   += blocks[i].ton;
             cuT    += blocks[i].ton * blocks[i].gCu;
             auT    += blocks[i].ton * blocks[i].gAu;
         } else  wasteT += blocks[i].ton;
     }
     double strip = oreT > 0 ? wasteT / oreT : 0;
     double gCu   = oreT > 0 ? cuT / oreT : 0;
     double gAu   = oreT > 0 ? auT / oreT : 0;
 
     /* Comprobar precedencias violadas */
     size_t bad = 0; ifstream fa2(ARCS); getline(fa2, line);
     while (getline(fa2, line)) {
         auto p = line.find(',');
         int idU = stoi(line.substr(0,p));
         int idV = stoi(line.substr(p+1));
         if (!idx.count(idU) || !idx.count(idV)) continue;
         int u = idx[idU]; int v = idx[idV];
         bool inPitU = pit[u], inPitV = pit[v];
         if (inPitV && !inPitU) ++bad;
     }
     cout << "Precedencias violadas: " << bad << '\n';
 
     /* Guardar resultados */
     ofstream ids(OUT_DIR + "pit_ids_pf.csv");
     for (int i = 0; i < n; ++i) if (pit[i]) ids << blocks[i].id << '\n';
     ids.close();
 
     ofstream fout(OUT_DIR + "pf_resultados.txt", ios::app);
     auto now = chrono::system_clock::to_time_t(chrono::system_clock::now());
     fout   << "----- " << ctime(&now)
            << fixed << setprecision(3)
            << "Beneficio: " << benefitInt / double(SCALE) << "\n"
            << "Rec Cu: "   << cuT << " | Au: " << auT << "\n"
            << "Strip: "    << strip << "\n"
            << "Ley: "      << gCu * 100 << "% Cu, "
                            << gAu << " g/t Au\n"
            << "Tiempo: "   << secs << " s\n\n";
     fout.close();
 
     cout  << fixed << setprecision(3)
           << "Beneficio: " << benefitInt / double(SCALE) << "\n"
           << "Rec Cu: "   << cuT << " | Au: " << auT << "\n"
           << "Strip: "    << strip << "\n"
           << "Ley: "      << gCu * 100 << "% Cu, "
                           << gAu << " g/t Au\n"
           << "Tiempo: "   << secs << " s\n"
           << "Bloques en pit: "
           << count(pit.begin(), pit.end(), 1) << '\n';
 }
 