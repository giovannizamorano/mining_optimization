// lg2_approach2.cpp - Implementacion integra del Approach 2 (cap. 5)
// g++ -O3 -std=c++17 lg2.cpp -o lg2
// Para depuracion: g++ -DDEBUG -O3 -std=c++17 -Wall -Wextra lg2.cpp -o lg2
// ---------------------------------------------------------------

#ifdef DEBUG
  // Macro DBG original para flexibilidad, pero usaremos una especifica para el resumen
  #define DBG_PRINT(x) do { std::cerr << "[DBG] " << x << std::endl; } while(0)
  // Macro para el resumen de una linea
  #define DBG_SUMMARY(x) do { if (DEBUG_SUMMARY_ENABLED) { std::cerr << "[SUM] " << x << std::endl; } } while(0)
#else
  #define DBG_PRINT(x) do {} while(0)
  #define DBG_SUMMARY(x) do {} while(0)
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <string>
#include <queue>
#include <stack>
#include <iomanip> // Para std::fixed y std::setprecision

using namespace std;

/* ---------- clave 3-D + hash ---------- */
struct Key3{ double x,y,z; bool operator==(const Key3&o)const noexcept{
    return x==o.x&&y==o.y&&z==o.z;}};
struct Key3Hash{
    size_t operator()(const Key3&k)const noexcept{
        size_t h=0x9e3779b97f4a7c15ULL; hash<double> hd;
        auto mix=[&](double v){h^=hd(v)+0x9e3779b9+(h<<6)+(h>>2);};
        mix(k.x); mix(k.y); mix(k.z); return h;}};

/* ---------- estructuras ---------- */
struct Node{ double w; double x,y,z; };
struct Arc { int u,v; bool plus=false,strong=false; };

// --- Variables Globales para Control de Depuracion ---
#ifdef DEBUG
bool DEBUG_SUMMARY_ENABLED = true; // Habilita/deshabilita los resumenes DBG_SUMMARY
bool DEBUG_DETAILED_IN_WINDOW = true; // Habilita DBG_PRINT detallados dentro de la ventana
const int DEBUG_TARGET_ITERATION = 0; // Poner a 0 para desactivar, o a la iteracion especifica
const int DEBUG_WINDOW = 0;           // Ventana alrededor del TARGET (ej. T = TARGET +/- WINDOW)
#else
// En modo no-DEBUG, estas no importan porque las macros DBG no hacen nada
bool DEBUG_SUMMARY_ENABLED = false;
bool DEBUG_DETAILED_IN_WINDOW = false;
const int DEBUG_TARGET_ITERATION = 0;
const int DEBUG_WINDOW = 0;
#endif


/* ---------- lectura CSV ---------- */
vector<Node> readCSV(const string&file,
                     unordered_map<Key3,int,Key3Hash>&idx,char sep = ','){
    ifstream f(file); if(!f){cerr<<"No se pudo abrir el archivo: " << file << "\n";exit(1);}
    string hdr; getline(f,hdr);
    vector<Node>N(1);
    string ln, tmp;
    int line_count = 1;
    while(getline(f, ln)){
        line_count++;
        if(ln.empty() || ln.find_first_not_of(" \t\n\v\f\r") == std::string::npos) continue;
        stringstream ss(ln);
        vector<string> tok;
        while(getline(ss, tmp, sep))
            tok.push_back(tmp);
        
        if (tok.size() < 6) {
            std::cerr << "Error en CSV (linea " << line_count << "): Columnas insuficientes. Esperadas >=6, encontradas " << tok.size() << ". Linea: " << ln << std::endl;
            continue;
        }
        try {
            double X = stod(tok[0]);
            double Y = stod(tok[1]);
            double Z = stod(tok[2]);
            double w = stod(tok[5]);
            int id = N.size();
            N.push_back({ w, X, Y, Z });
            idx[{X,Y,Z}] = id;
        } catch (const std::exception& e) {
            std::cerr << "Error procesando CSV (linea " << line_count << "): " << ln << " | Error: " << e.what() << std::endl;
        }
    }
    DBG_PRINT("readCSV: Leidos " << N.size()-1 << " bloques desde " << file);
    return N;
}

/* ---------- diferencia minima positiva ---------- */
template<class T>
double minDiff(std::vector<T> v){
    if (v.size() < 2) return std::numeric_limits<double>::infinity();
    std::sort(v.begin(), v.end());
    double d = std::numeric_limits<double>::infinity();
    bool found_diff = false;
    for(size_t i = 1; i < v.size(); ++i){
        double t = v[i] - v[i-1];
        if (t > 1e-9) {
             if (t < d) d = t;
             found_diff = true;
        }
    }
    return found_diff ? d : std::numeric_limits<double>::infinity();
}

/* ---------- aristas 1-5 ---------- */
vector<pair<int,int>> arcs1_5(
        const vector<Node>&N,
        const unordered_map<Key3,int,Key3Hash>&idx)
{
    DBG_PRINT("arcs1_5: Iniciando generacion de arcos de precedencia geometrica (feas).");
    vector<double>xs,ys,zs;
    for(size_t i=1;i<N.size();++i){
        xs.push_back(N[i].x); ys.push_back(N[i].y); zs.push_back(N[i].z);}
    
    double dx=minDiff(xs); 
    double dy=minDiff(ys); 
    double dz=minDiff(zs);

    double default_spacing = 1.0; // Un valor por defecto si no se pueden calcular las diferencias
    if (!std::isfinite(dx)) { DBG_PRINT("arcs1_5: dx es infinito/NaN. Usando fallback."); dx = (std::isfinite(dy) ? dy : (std::isfinite(dz) ? dz : default_spacing)); }
    if (!std::isfinite(dy)) { DBG_PRINT("arcs1_5: dy es infinito/NaN. Usando fallback."); dy = (std::isfinite(dx) ? dx : (std::isfinite(dz) ? dz : default_spacing)); }
    if (!std::isfinite(dz)) { DBG_PRINT("arcs1_5: dz es infinito/NaN. Usando fallback."); dz = default_spacing; }


    vector<pair<int,int>>E;
    for(size_t id_actual=1; id_actual<N.size(); ++id_actual){
        const auto&bloque_actual=N[id_actual]; 
        double z_predecesor = bloque_actual.z - dz; 

        Key3 predecesores_coords[5]={
            {bloque_actual.x, bloque_actual.y, z_predecesor},            
            {bloque_actual.x - dx, bloque_actual.y, z_predecesor},      
            {bloque_actual.x + dx, bloque_actual.y, z_predecesor},      
            {bloque_actual.x, bloque_actual.y - dy, z_predecesor},      
            {bloque_actual.x, bloque_actual.y + dy, z_predecesor}       
        };
        for(const auto& coord_pred : predecesores_coords){
            auto it=idx.find(coord_pred);
            if(it!=idx.end()){ 
                E.push_back({(int)id_actual, it->second}); 
            }
        }
    }
    DBG_PRINT("arcs1_5: Generados " << E.size() << " arcos feas.");
    return E;
}

/* ---------- etiquetar (Tabla 5.15) ---------- */
void label(const vector<Node>& N,
           vector<Arc>&        A,
           vector<double>&  acc,
           int current_step, bool detailed_debug_active) // Pasar bandera de depuracion detallada
{
    if (detailed_debug_active) DBG_PRINT("label (T=" << current_step << "): Entrando. A.size()=" << A.size() << ", N.size()=" << N.size());
    
    int n = N.size()-1; 
    if (n < 0) { 
        acc.clear();
        if (detailed_debug_active) DBG_PRINT("label (T=" << current_step << "): No hay nodos de datos.");
        return;
    }
    acc.assign(n + 1, 0.0);

    if (detailed_debug_active) DBG_PRINT("label (T=" << current_step << "): Poblando adj. n=" << n);
    vector<vector<int>> adj(n + 1);
    for(const auto& a_arc : A){ 
        if (a_arc.u < 0 || a_arc.u > n || a_arc.v < 0 || a_arc.v > n) {
            std::cerr << "[ERROR FATAL LABEL] (T=" << current_step
                      << ") Arco con indices fuera de rango para adj! Arc: ("
                      << a_arc.u << "," << a_arc.v << ") donde n=" << n << std::endl;
            continue; 
        }
        adj[a_arc.u].push_back(a_arc.v);
        adj[a_arc.v].push_back(a_arc.u);
    }
    if (detailed_debug_active) DBG_PRINT("label (T=" << current_step << "): Termino de poblar adj.");

    vector<int> depth(n + 1, -1);
    queue<int> q_bfs;
    if (n >= 0) { 
       depth[0] = 0;
       q_bfs.push(0);
    }
    while(!q_bfs.empty()){
        int u=q_bfs.front(); q_bfs.pop();
        if (u < 0 || u > n) continue; 
        for(int v: adj[u]) {
            if (v < 0 || v > n) continue; 
            if(depth[v]==-1){ depth[v]=depth[u]+1; q_bfs.push(v); }
        }
    }
    if (detailed_debug_active) DBG_PRINT("label (T=" << current_step << "): BFS para depth completado.");

    // --- DFS ITERATIVO ---
    if (detailed_debug_active) DBG_PRINT("label (T=" << current_step << "): Iniciando DFS ITERATIVO para acc.");
    int max_dfs_simulated_depth = 0; 
    std::fill(acc.begin(), acc.end(), 0.0); 
    std::stack<pair<int, int>> dfs_iter_stack; 
    std::vector<int> post_order_nodes;   
    post_order_nodes.reserve(n + 1);     
    std::vector<int> dfs_parent_map(n + 1, -1); 
    if (n >= 0) dfs_iter_stack.push({0, 0}); 
    int current_simulated_depth = 0; 
    while(!dfs_iter_stack.empty()){
        int u = dfs_iter_stack.top().first;
        int &child_idx = dfs_iter_stack.top().second; 
        if (child_idx == 0 && u>=0 && (size_t)u < N.size()) { 
            current_simulated_depth++;
            if (current_simulated_depth > max_dfs_simulated_depth) max_dfs_simulated_depth = current_simulated_depth;
        }
        bool processed_a_child_this_iteration = false;
        if(u>=0 && (size_t)u < adj.size()){
            while(child_idx < adj[u].size()){
                int v = adj[u][child_idx];
                child_idx++; 
                if (v == dfs_parent_map[u]) continue;
                if(v>=0 && (size_t)v < depth.size() && u>=0 && (size_t)u < depth.size()){
                    if (depth[v] == depth[u] + 1) { 
                        dfs_parent_map[v] = u;      
                        dfs_iter_stack.push({v, 0});   
                        processed_a_child_this_iteration = true;
                        break; 
                    }
                } else if (detailed_debug_active) {
                     std::cerr << "[ADVERTENCIA DFS-ITER] (T=" << current_step << ") Indice v=" << v << " o u=" << u << " fuera de rango para depth." << std::endl;
                }
            }
        }
        if(!processed_a_child_this_iteration){
            post_order_nodes.push_back(u);
            current_simulated_depth--; 
            dfs_iter_stack.pop(); 
        }
    }
    for (int u_node : post_order_nodes) {
        if ((size_t)u_node >= N.size() || u_node < 0) continue; 
        acc[u_node] = N[u_node].w; 
        if ((size_t)u_node >= adj.size()) continue; 
        for (int v_neighbor : adj[u_node]) {
            if ((size_t)v_neighbor >= N.size() || v_neighbor < 0) continue; 
            if (dfs_parent_map[v_neighbor] == u_node) { 
                if ((size_t)v_neighbor < acc.size()) acc[u_node] += acc[v_neighbor];
            }
        }
    }
    if (detailed_debug_active) DBG_PRINT("label (T=" << current_step << "): DFS (iterativo) para acc completado. Max prof simulada: " << max_dfs_simulated_depth);
    // --- FIN DFS ITERATIVO ---

    for(auto& a_arc:A){
        if (a_arc.u < 0 || a_arc.u > n || a_arc.v < 0 || a_arc.v > n ||
            (size_t)a_arc.v >= depth.size() || (size_t)a_arc.u >= depth.size() ||
            (size_t)a_arc.v >= acc.size() || (size_t)a_arc.u >= acc.size() ) { 
            std::cerr << "[ERROR FATAL LABEL - PLUS/STRONG] (T=" << current_step 
                      << ") Indices fuera de rango para depth/acc. Arc: ("
                      << a_arc.u << "," << a_arc.v << ")" << std::endl;
            continue;
        }
        a_arc.plus   = (depth[a_arc.v] > depth[a_arc.u]);
        double branch_node_idx = a_arc.plus ? a_arc.v : a_arc.u;
        // Asegurar que branch_node_idx es valido para acc
        double branch_value = (branch_node_idx >=0 && (size_t)branch_node_idx < acc.size()) ? acc[branch_node_idx] : 0.0;
        if (!(branch_node_idx >=0 && (size_t)branch_node_idx < acc.size()) && detailed_debug_active) {
             std::cerr << "[ADVERTENCIA LABEL - PLUS/STRONG] (T=" << current_step << ") Indice branch_node_idx=" << branch_node_idx 
                       << " invalido para acc. Usando branch_value=0. Arc: (" << a_arc.u << "," << a_arc.v << ")" << std::endl;
        }
        a_arc.strong = ( (a_arc.plus && branch_value>=0) || (!a_arc.plus && branch_value<=0) );
    }
    if (detailed_debug_active) DBG_PRINT("label (T=" << current_step << "): Asignacion plus/strong completada.");
    // No hay DBG_PRINT aqui fuera de la ventana para reducir ruido. El DBG_SUMMARY en main informara.
}

/* ---------- normalizar (tree-cut) ---------- */
int normalize(const vector<Node>& N, // Devuelve el numero de loops internos
               vector<Arc>&        A,
               vector<double>&  acc,
               int current_step, bool detailed_debug_active)
{
    if (detailed_debug_active) DBG_PRINT("normalize (T=" << current_step << "): Iniciando. A.size()=" << A.size());
    
    bool again = true;
    int normalize_loops = 0;
    std::unordered_set<int> nodes_connected_to_root_v; 
    if (A.size() > 0) { // Solo poblar si A no esta vacia
        for(const auto& arc_iter : A) if(arc_iter.u == 0) nodes_connected_to_root_v.insert(arc_iter.v);
    }


    while (again) {
        again = false;
        normalize_loops++;
        if (normalize_loops > std::max(200, (int)A.size() / 10) && A.size() > 0 ) { // Limite mas dinamico
            std::cerr << "[ADVERTENCIA NORMALIZE] (T=" << current_step << ") Bucle de normalizacion (" << normalize_loops << " loops) parece excesivo, rompiendo. A.size()=" << A.size() << std::endl;
            break;
        }

        label(N, A, acc, current_step, detailed_debug_active); // Pasar la bandera
        if (detailed_debug_active && normalize_loops == 1) DBG_PRINT("normalize (T=" << current_step << "): Label completado dentro de normalize loop " << normalize_loops);

        for (size_t i = 0; i < A.size(); ) { 
            Arc current_A_arc = A[i]; // Copia para inspeccion
            bool erased_this_iteration = false;

            if (current_A_arc.strong && !current_A_arc.plus) {
                if(detailed_debug_active) DBG_PRINT("normalize (T=" << current_step << "): Accion 1 (strong-minus) en arco (" << current_A_arc.u << "," << current_A_arc.v << ")");
                int q_node = current_A_arc.u;
                A.erase(A.begin() + i); 
                // No es necesario quitar de nodes_connected_to_root_v porque el arco borrado no era (0,X)
                if(nodes_connected_to_root_v.find(q_node) == nodes_connected_to_root_v.end()){
                    A.push_back({0,q_node,true,false});
                    nodes_connected_to_root_v.insert(q_node); // Anadir al set
                }
                again = true; erased_this_iteration = true; break; 
            }

            if (current_A_arc.strong && current_A_arc.plus && current_A_arc.u != 0) {
                 if(detailed_debug_active) DBG_PRINT("normalize (T=" << current_step << "): Accion 2 (strong-plus no root) en arco (" << current_A_arc.u << "," << current_A_arc.v << ")");
                int v_node = current_A_arc.v;
                A.erase(A.begin() + i); 
                if(nodes_connected_to_root_v.find(v_node) == nodes_connected_to_root_v.end()){
                    A.push_back({0,v_node,true,false});
                    nodes_connected_to_root_v.insert(v_node); // Anadir al set
                }
                again = true; erased_this_iteration = true; break; 
            }

            if (!erased_this_iteration) i++; 
        } 
        // No es necesario reconstruir el set aqui si 'again' es falso y salimos.
        // Si 'again' es true, el bucle continua, A se re-etiqueta.
        // El set nodes_connected_to_root_v podria desincronizarse si un arco (0,X) se borra por otra logica (no ocurre aqui).
        // Es mas seguro reconstruirlo si 'again' es true y A ha cambiado.
        if (again) {
             nodes_connected_to_root_v.clear();
             for(const auto& arc_iter : A) if(arc_iter.u == 0) nodes_connected_to_root_v.insert(arc_iter.v);
        }
    } 
    if (detailed_debug_active) DBG_PRINT("normalize (T=" << current_step << "): Terminado. A.size()=" << A.size() << ", loops=" << normalize_loops);
    return normalize_loops;
}

/* ---------- *** MOD 2: construir Y completo (sub-rama) ---------- */
vector<char> buildY(const vector<Arc>&A, const vector<Node>&N, int current_step, bool detailed_debug_active)
{
    if (detailed_debug_active) DBG_PRINT("buildY (T=" << current_step << "): Entrando. A.size()=" << A.size() << ", N.size()=" << N.size());
    
    int n = N.size()-1;
    if (n < 0) {
        if (detailed_debug_active) DBG_PRINT("buildY (T=" << current_step << "): No hay nodos de datos.");
        return vector<char>();
    }

    vector<vector<int>> adj(n+1);
    for(const auto& a_arc:A){
         if (a_arc.u < 0 || a_arc.u > n || a_arc.v < 0 || a_arc.v > n) {
             std::cerr << "[ERROR FATAL BUILDY - PRE-ADJ-POP] (T=" << current_step 
                       << ") Arco con indices fuera de rango para adj! Arc: ("
                       << a_arc.u << "," << a_arc.v << ") donde n=" << n << std::endl;
             continue;
         }
        adj[a_arc.u].push_back(a_arc.v);
        adj[a_arc.v].push_back(a_arc.u);
    }

    vector<int> depth(n+1,-1); queue<int> q_bfs;
    if (n >= 0) {
        depth[0]=0; q_bfs.push(0);
    }
    while(!q_bfs.empty()){
        int u=q_bfs.front(); q_bfs.pop();
        if (u < 0 || u > n) continue;
        for(int v:adj[u]) {
            if (v < 0 || v > n) continue;
            if(depth[v]==-1){ depth[v]=depth[u]+1; q_bfs.push(v); }
        }
    }
    if (detailed_debug_active) DBG_PRINT("buildY (T=" << current_step << "): BFS para depth completado.");

    vector<char> inY(N.size(),0); queue<int> yq;
    for(const auto& a_arc:A) {
        if(a_arc.plus && a_arc.strong){ 
            if (a_arc.v >= 0 && (size_t)a_arc.v < inY.size()) { 
                 if(inY[a_arc.v] == 0){ 
                    inY[a_arc.v]=1; 
                    yq.push(a_arc.v);
                 }
            } else {
                std::cerr << "[ERROR FATAL BUILDY ASSIGN INY] (T=" << current_step << ") Indice a.v=" << a_arc.v
                          << " invalido para arco plus/strong. N.size()=" << N.size() << std::endl;
            }
        }
    }

    while(!yq.empty()){
        int u=yq.front(); yq.pop();
        if (u < 0 || (size_t)u >= N.size()) continue; 
        for(int w:adj[u]) {
            if (w < 0 || (size_t)w >= inY.size() || (size_t)w >= depth.size() || (size_t)u >= depth.size() ) { 
                 std::cerr << "[ERROR FATAL BUILDY LOOP INY] (T=" << current_step << ") Indice w=" << w << " o u=" << u
                           << " invalido. N.size()=" << N.size() << " depth.size()=" << depth.size() << std::endl;
                 continue;
            }
            if (depth[w]>depth[u] && !inY[w]){
                inY[w]=1; yq.push(w);
            }
        }
    }
    if (detailed_debug_active) DBG_PRINT("buildY (T=" << current_step << "): Cola de expansion para inY procesada.");
    return inY;
}

/* ---------- *** MOD 3: cierre = Suma pesos de nodos en Y ---------- */
double closure(const vector<Arc>&A, const vector<Node>&N, int current_step, bool detailed_debug_active)
{
    if (detailed_debug_active) DBG_PRINT("closure (T=" << current_step << "): Calculando...");
    
    vector<char> inY = buildY(A,N, current_step, detailed_debug_active);
    double s=0;
    int nodes_in_y_count = 0;
    for(size_t i=1;i<N.size();++i) { 
        if(i < inY.size() && inY[i]) { 
             if (i < N.size()) {
                 s+=N[i].w;
                 nodes_in_y_count++;
             }
        }
    }
    // El DBG_SUMMARY en main se encargará de esto.
    // if (detailed_debug_active) DBG_PRINT("closure (T=" << current_step << "): Valor=" << s << ", Nodos en Y=" << nodes_in_y_count);
    return s;
}

/* ====================  main  ==================== */
int main(int argc,char*argv[]){
    std::ios_base::sync_with_stdio(false); 
    std::cin.tie(NULL);                    
    std::cout << std::fixed << std::setprecision(2); // Para imprimir doubles con 2 decimales
    std::cerr << std::fixed << std::setprecision(2);


    if(argc<2){ cerr<<"Uso: "<<argv[0]<<" archivo.csv\n"; return 1; }

    unordered_map<Key3,int,Key3Hash> idx;
    vector<Node> N   = readCSV(argv[1], idx);
    if (N.size() <= 1) { 
        cerr << "No se leyeron datos de bloques del CSV o el archivo esta vacio." << endl;
        return 1;
    }

    vector<pair<int,int>> feas = arcs1_5(N, idx);
    sort(feas.begin(),feas.end(),
         [&](const pair<int,int>&a, const pair<int,int>&b){ 
            bool idx_a_valid = (a.first >= 0 && (size_t)a.first < N.size());
            bool idx_b_valid = (b.first >= 0 && (size_t)b.first < N.size());
            if (idx_a_valid && idx_b_valid) return N[a.first].w > N[b.first].w;
            if (idx_a_valid) return true; 
            return false;
          });


    vector<Arc> A;
    for(size_t i=1;i<N.size();++i) A.push_back({0,(int)i,true,false});
    vector<double> acc;
    int normalize_loops_t0 = normalize(N,A,acc, 0, false); // T0, no depuracion detallada
    double closure_t0 = closure(A,N, 0, false);
    DBG_PRINT("main: Normalizacion inicial (T0) completada. A.size=" << A.size() << ", Loops Norm=" << normalize_loops_t0 << ", Cierre T0=" << closure_t0);


    std::unordered_set<int> nodes_with_dummy_arc_to_root; 
    for(const auto& arc_iter : A) if(arc_iter.u == 0) nodes_with_dummy_arc_to_root.insert(arc_iter.v);

    auto removeDummyArc = [&](int node_v, vector<Arc>& current_A, std::unordered_set<int>& current_dummies) {
        if (current_dummies.count(node_v)) { // El nodo esta en el set de dummies
            // Buscar y borrar el arco (0, node_v) de A
            for (size_t k = 0; k < current_A.size(); ++k) {
                if (current_A[k].u == 0 && current_A[k].v == node_v) {
                    current_A.erase(current_A.begin() + k);
                    current_dummies.erase(node_v); // Quitar del set
                    return true; 
                }
            }
        }
        return false; 
    };


    int step=1;
    DBG_PRINT("main: Iniciando bucle principal de iteraciones.");
    if (DEBUG_SUMMARY_ENABLED) {
        std::cerr << "[SUM] Iter  | Cierre      | A.size | N.Loops | Arco Anadido (u->v) | u_inY | v_inY | N[u].w | N[v].w" << std::endl;
        std::cerr << "[SUM] ------|-------------|--------|---------|---------------------|-------|-------|--------|--------" << std::endl;
    }

    int last_u_added = -1, last_v_added = -1;
    int repeat_add_count = 0;

    while(true){
        // ... (calculo de detailed_debug_this_step como antes) ...
        bool detailed_debug_this_step = false;
        #ifdef DEBUG
        if (DEBUG_DETAILED_IN_WINDOW && 
            DEBUG_TARGET_ITERATION > 0 &&
            step >= DEBUG_TARGET_ITERATION - DEBUG_WINDOW &&
            step <= DEBUG_TARGET_ITERATION + DEBUG_WINDOW) {
            detailed_debug_this_step = true;
        }
        #endif

        if (!detailed_debug_this_step && step % 50 == 0 && step > 0) { // Ajusta la frecuencia del "latido"
            DBG_PRINT("... Procesando Iteracion T=" << step << " / A.size=" << A.size() << " ...");
        }

        vector<char> inY = buildY(A,N, step, detailed_debug_this_step);
        // if(detailed_debug_this_step) DBG_PRINT("  (T=" << step << ") buildY completado. inY.size()=" << inY.size());
        
        bool added_arc_this_step=false; 
        int u_current_add = -1, v_current_add = -1;
        int normalize_loops_this_step = 0;

        int feas_processed_count = 0;
        for(const auto& current_feas_pair : feas){ 
            feas_processed_count++;
            int u_node = current_feas_pair.first;
            int v_node = current_feas_pair.second; 
                                                  
            if (u_node <= 0 || (size_t)u_node >= N.size() || v_node < 0 || (size_t)v_node >= N.size()) {
                 std::cerr << "[ERROR FATAL MAIN LOOP INVALID FEAS NODEID] (T=" << step << ") u o v en feas con ID invalido. u=" << u_node << ", v=" << v_node
                           << " N.size()=" << N.size() << std::endl;
                 continue;
            }
            if ((size_t)u_node >= inY.size() || (size_t)v_node >= inY.size()) { 
                 std::cerr << "[ERROR FATAL MAIN LOOP ACCESS INY] (T=" << step << ") u o v fuera de rango para inY. u=" << u_node << ", v=" << v_node
                           << " inY.size()=" << inY.size() << std::endl;
                 continue;
            }

            // El DBG detallado de "Probando arco feas" solo si detailed_debug_this_step es true
            if (detailed_debug_this_step && feas_processed_count % 1000 == 0) { 
                DBG_PRINT("  (T=" << step << ") Probando arco feas #" << feas_processed_count << ": " << u_node << "(abajo) ->" << v_node << "(arriba)"
                    << " inY[u]=" << int(inY[u_node])
                    << " inY[v]=" << int(inY[v_node]));
            }
            
            if(inY[u_node] && !inY[v_node]){ // Condicion original
                u_current_add = u_node; v_current_add = v_node; 

                if (detailed_debug_this_step) {
                     DBG_PRINT("  (T=" << step << ") >> Condicion cumplida para feas: " << u_node << "(en Y) ->" << v_node << "(no en Y).");
                }


                bool dummy_u_removed = removeDummyArc(u_node, A, nodes_with_dummy_arc_to_root);
                if (!dummy_u_removed) { 
                    removeDummyArc(v_node, A, nodes_with_dummy_arc_to_root);
                }
                
                A.push_back({u_node, v_node, true, false}); 
                
                normalize_loops_this_step = normalize(N,A,acc, step, detailed_debug_this_step); 

                nodes_with_dummy_arc_to_root.clear(); 
                for(const auto& arc_iter : A) if(arc_iter.u == 0) nodes_with_dummy_arc_to_root.insert(arc_iter.v);
                
                added_arc_this_step=true; break;
            }
        } // Fin del bucle feas

        double current_closure_val = closure(A,N, step, detailed_debug_this_step);
        
        // --- IMPRESION DE RESUMEN POR ITERACION ---
        std::stringstream summary_ss;
        summary_ss << std::setw(7) << ("T" + std::to_string(step)) << "| "
                   << std::setw(11) << current_closure_val << " | "
                   << std::setw(6) << A.size() << " | "
                   << std::setw(7) << (added_arc_this_step ? normalize_loops_this_step : 0) << " | ";
        if(added_arc_this_step){
            summary_ss << std::setw(5) << u_current_add << "->" << std::setw(5) << v_current_add << " | "
                       << std::setw(5) << int(inY[u_current_add]) << " | " // Valor de inY ANTES de añadir este arco y normalizar
                       << std::setw(5) << int(inY[v_current_add]) << " | "
                       << std::setw(6) << N[u_current_add].w << " | "
                       << std::setw(6) << N[v_current_add].w;

            if (u_current_add == last_u_added && v_current_add == last_v_added) {
                repeat_add_count++;
            } else {
                last_u_added = u_current_add;
                last_v_added = v_current_add;
                repeat_add_count = 1;
            }
            if (repeat_add_count > 5) { // Si el mismo arco se añade mas de 5 veces seguidas
                summary_ss << " [CICLO ARCO REPETIDO " << repeat_add_count << " veces!]";
            }

        } else {
            summary_ss << "No hay arco anadido";
             last_u_added = -1; last_v_added = -1; repeat_add_count = 0; // Resetear contador de ciclo
        }
        DBG_SUMMARY(summary_ss.str());

        // Impresion normal a cout solo si se anadio un arco
        if(added_arc_this_step){
             cout << "T" << step << "  cierre = " << current_closure_val << "\n";
        }

        if(!added_arc_this_step){
            DBG_PRINT("main: No se anadieron mas arcos en T=" << step << ". Terminando bucle while(true).");
            break;
        }
        
        if (step > std::max(40000, (int)N.size() * 2) && N.size() > 100 ) { 
             std::cerr << "[ADVERTENCIA MAIN] (T=" << step << ") Bucle principal (" << step << " iteraciones) parece excesivo, rompiendo." << std::endl;
             break;
        }
        if (repeat_add_count > 20 && current_closure_val < 13000) { // Condicion de parada por ciclo detectado y cierre bajo
             std::cerr << "[ADVERTENCIA MAIN] (T=" << step << ") Detectado posible ciclo con arco " << u_current_add << "->" << v_current_add
                       << " repitiendose " << repeat_add_count << " veces con cierre=" << current_closure_val << ". Terminando." << std::endl;
             break;
        }
        ++step;
    }

    cout << "\nMaximo cierre = " << closure(A,N, step, false) << "\n"; 
    return 0;
}