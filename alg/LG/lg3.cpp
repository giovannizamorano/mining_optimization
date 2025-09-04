// g++ -O3 -std=c++17 lg2.cpp -o lg2
// ---------------------------------------------------------------

#ifdef DEBUG
  #define DBG(x) do { std::cerr << "[DBG] " << x << std::endl; } while(0)
#else
  #define DBG(x) do {} while(0)
#endif

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
#include <string>
#include <queue>
#include <cmath>
#include <stack>
#include <unordered_set>
#include <iomanip>
#include <cmath>
#include <cstdint>
using namespace std;

/* =======================   ARISTAS PROCESADAS   ======================= */
inline uint64_t arcKey(int u,int v){ return ( (uint64_t)u << 32 ) | (uint32_t)v; }
/* aristas que ya se resolvieron (FORZAR o PODAR) */
static std::unordered_set<uint64_t> procesadas;


/* ---------- clave 3-D + hash ---------- */
struct Key3{ double x,y,z; bool operator==(const Key3&o)const noexcept{return x==o.x&&y==o.y&&z==o.z;}};
struct Key3Hash{
    size_t operator()(const Key3&k)const noexcept{
        size_t h=0x9e3779b97f4a7c15ULL; hash<double> hd;
        auto mix=[&](double v){h^=hd(v)+0x9e3779b9+(h<<6)+(h>>2);};
        mix(k.x); mix(k.y); mix(k.z); return h;}};

/* ---------- estructuras ---------- */
struct Node{double w,x,y,z; };
struct Arc { int u,v; bool plus=false,strong=false; };

/* ---------- auxiliar: rehacer lista de adyacencia sin realocar ---------- */
static void buildAdj(const vector<Arc>& A, vector<vector<int>>& adj)
{
    for (auto& row : adj) row.clear();
    for (const auto& a : A) {adj[a.u].push_back(a.v); adj[a.v].push_back(a.u);
    } // Asume que A es un grafo no dirigido
}

/* ---------- auxiliar: adyacencia dirigida ---------- */
static void buildAdjDir(const vector<Arc>& A,
                        vector<vector<int>>& child,   // sólo hijos (plus)
                        vector<int>& parent)          // padre de cada nodo
{
    for (auto& row : child) row.clear();
    fill(parent.begin(), parent.end(), -1);

    for (const auto& a : A) {
        if (a.plus) {                 // orientación padre → hijo
            child[a.u].push_back(a.v);
            parent[a.v] = a.u;
        } else {
            child[a.v].push_back(a.u);
            parent[a.u] = a.v;
        }
    }
    parent[0] = -1;                   // la raíz no tiene padre
}


/* ---------- lectura del archivo ---------- */
vector<Node>
readBlocks(const string& file,unordered_map<Key3,int,Key3Hash>& coord2idx, unordered_map<int,int>& id2idx)
{
    std::ifstream f(file);
    if (!f) { std::cerr << "No pude abrir " << file << '\n'; exit(1); }

    std::vector<Node> N;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '%') continue;      // comentario

        // Split usando espacios, tabs o ':'
        std::stringstream ss(line);
        std::string tok;
        std::vector<std::string> t;
        while (std::getline(ss, tok, ' '))                 // primer split
            if (!tok.empty()) t.push_back(tok);

        // Si algunos campos vienen pegados con ':', haz un segundo split
        std::vector<std::string> toks;
        for (auto& s : t) {
            size_t pos;
            while ((pos = s.find(':')) != std::string::npos) {
                toks.push_back(s.substr(0, pos));
                s.erase(0, pos + 1);
            }
            toks.push_back(s);
        }

        if (toks.size() < 5) continue; // línea malformada

        int    id = std::stoi(toks[0]);
        double X  = std::stod(toks[1]);
        double Y  = std::stod(toks[2]);
        double Z  = std::stod(toks[3]);
        double w  = std::stod(toks[9]);
        //double value_ex   = std::stod(toks[8]);
        //double value_proc = std::stod(toks[9]);
        //double w          = value_proc - value_ex;;

/*
to read marvin web example
        double tonn = std::stod(toks[4]);
        double pp   = std::stod(toks[7]);   // proc_profit (US$/t)
        double w    = (pp - 0.9) * tonn; 
*/
        int idx = (int)N.size();
        id2idx[id]          = idx;
        coord2idx[{X,Y,Z}]  = idx;
        N.push_back({w, X, Y, Z});
    }
    return N;
}


/* ---------- diferencia mínima positiva ---------- */
template<class T>
double minDiff(std::vector<T> v){
    std::sort(v.begin(), v.end());
    double d = std::numeric_limits<double>::infinity();
    for(size_t i = 1; i < v.size(); ++i){
        double t = v[i] - v[i-1];
        if (t > 1e-9 && t < d) 
            d = t;
    }
    return d;
}

/* ---------- aristas 1-5 ---------- */
vector<pair<int,int>> arcs1_5( const vector<Node>&N, const unordered_map<Key3,int,Key3Hash>&idx){

    vector<double>xs,ys,zs;
    for(size_t i=1;i<N.size();++i){xs.push_back(N[i].x); ys.push_back(N[i].y); zs.push_back(N[i].z);}
    double dx=minDiff(xs), dy=minDiff(ys), dz=minDiff(zs);
    if (!std::isfinite(dy) || dy > 1e50) {
        DBG("dy no válida (infinito o muy grande), ajustando dy = dx");
        dy = dx;
    }
 
    DBG("dx=" << dx << "  dy=" << dy << "  dz=" << dz);

    vector<pair<int,int>>E;
    for(size_t id=1; id<N.size(); ++id){
        const auto&b=N[id]; double up=b.z+dz;
//        Key3 p[9] = {
//            { b.x - dx, b.y - dy, up },  // esquina NW
//            { b.x      , b.y - dy, up },  // norte
//            { b.x + dx, b.y - dy, up },  // esquina NE
//            { b.x - dx, b.y     , up },  // oeste
//            { b.x     , b.y     , up },  // centro (directamente arriba)
//            { b.x + dx, b.y     , up },  // este
//            { b.x - dx, b.y + dy, up },  // esquina SW
//            { b.x     , b.y + dy, up },  // sur
//            { b.x + dx, b.y + dy, up }   // esquina SE
//        };

        Key3 p[5]={{b.x,b.y,up},{b.x-dx,b.y,up},{b.x+dx,b.y,up},
                   {b.x,b.y-dy,up},{b.x,b.y+dy,up}};
/*
        Key3 p[1] = { { b.x, b.y, up } };
*/
        for(auto q:p){ 
            auto it=idx.find(q);
                DBG("Buscando predecesor en ("
                << q.x << "," << q.y << "," << q.z << ") => "
                << (it != idx.end() ? "OK" : "NO"));
            if(it!=idx.end()) E.push_back({id,it->second}); }
    }
    DBG("Arcos 1-5 generados: " << E.size());
    return E;
}

/* ---------- etiquetar---------- */
void label(const std::vector<Node>& N,
           const std::vector<std::vector<int>>& adj, // Esto es adj_cache desde normalize
           std::vector<Arc>&        A,             // Usado para asignar plus/strong al final
           std::vector<double>&     acc)           // Vector de pesos acumulados a rellenar
{
    int n = N.size() - 1;
    acc.assign(n + 1, 0.0);

    DBG("-- label: N nodos=" << N.size() - 1
        << "  A arcos=" << A.size());

    /* --- BFS para profundidad --- */
    static std::vector<int> depth;
    static std::vector<int> bfs;
    static std::vector<int> parent;

    if (parent.size() != N.size()) parent.resize(N.size());

    if (depth.size() != N.size()) {depth.resize(N.size()); bfs.resize(N.size());}

    std::fill(depth.begin(), depth.end(), -1);
    int head = 0, tail = 0;
    depth[0] = 0;
    parent[0] = -1; // Asignar padre ficticio para la raíz
    bfs[tail++] = 0;

    while (head < tail) {
        int u = bfs[head++];
        for (int v : adj[u])
            if (depth[v] == -1) {
                depth[v] = depth[u] + 1;
                parent[v] = u; // Asignar padre
                bfs[tail++] = v;
            }
    }

    DBG("Profundidades de nodos:");
    for (int i = 0; i <= n; ++i)
        DBG(" depth[" << i << "] = " << depth[i]);

    /* --- DFS para peso acumulado de cada sub-árbol --- */
    // === INICIO DEL NUEVO BLOQUE DFS RECURSIVO ===

    static std::vector<bool> visited_dfs_acc; 
    if (visited_dfs_acc.size() != N.size()) { visited_dfs_acc.resize(N.size());}
    std::fill(visited_dfs_acc.begin(), visited_dfs_acc.end(), false);

    // 2. Definir la función lambda recursiva para el DFS.

    std::function<double(int, int)> calculate_subtree_weight_dfs =
        [&](int u, int p_parent) -> double {
        visited_dfs_acc[u] = true; // Marcar el nodo actual como visitado
        double current_sum_weights = N[u].w;

        for (int v_neighbor : adj[u]) {
            if (v_neighbor == p_parent) continue; // No regresar al padre inmediato

            if (!visited_dfs_acc[v_neighbor]) { // Si el vecino no ha sido visitado
                current_sum_weights += calculate_subtree_weight_dfs(v_neighbor, u); // Sumar el peso acumulado del subárbol del vecino
            }
        }
        return acc[u] = current_sum_weights;
    };

    // 3. Llamada inicial al DFS.
    //    Se asume que el nodo 0 es la raíz del árbol.
    //    Verifica que N no esté vacío y que el nodo 0 exista antes de llamar.
    if (!N.empty() && 0 < N.size() && !visited_dfs_acc[0]) { // Asegurar que el nodo 0 es un índice válido
         calculate_subtree_weight_dfs(0, -1); // Iniciar DFS desde el nodo 0, con -1 como padre ficticio
    }
    // Si tu grafo pudiera tener múltiples componentes desconectados y necesitas
    // calcular 'acc' para todos ellos, o si el nodo 0 no siempre es la raíz:
    // for (int i = 0; i <= n; ++i) { // Iterar sobre todos los nodos posibles
    //     // Podrías necesitar una condición para N[i] si algunos índices no son nodos válidos
    //     if (i < N.size() && !visited_dfs_acc[i]) {
    //         calculate_subtree_weight_dfs(i, -1);
    //     }
    // }
    // === FIN DEL NUEVO BLOQUE DFS RECURSIVO ===

    /* --- imprimir pesos acumulados --- */
    DBG("Acc (peso sub-árbol) de cada nodo:");
    for (int i = 0; i <= n; ++i)
        DBG(" acc[" << i << "] = " << acc[i]);

    /* --- asignar plus/minus y strong/weak--- */
    for (auto& arc_ref : A) {
        arc_ref.plus = (parent[arc_ref.v] == arc_ref.u);
        double branch_weight = acc[arc_ref.plus ? arc_ref.v : arc_ref.u];
        arc_ref.strong = ((arc_ref.plus && branch_weight >0.0) ||
                          (!arc_ref.plus && branch_weight <=0.0));
    }
}
/************  UTILIDADES DUMMY  +  BACKUP DE ARISTAS  ************/
int posDummy(const std::vector<Arc>& A, int node){
    for(size_t k=0;k<A.size();++k)
        if(A[k].u==0 && A[k].v==node) return static_cast<int>(k);
    return -1;
}

void eliminaDummy(std::vector<Arc>& A,int node){
    int k = posDummy(A,node);
    if(k!=-1) A.erase(A.begin()+k);
}

// Inserta dummy (0 -> node) **solo si** no existe ya
void addDummy(std::vector<Arc>& A, int node){
    if(posDummy(A,node)==-1)
        A.push_back({0,node,true,false});
}

/* ---------- normalizar (tree-cut) ---------- */

void normalize(const vector<Node>& N,
               vector<Arc>&        A,
               vector<double>&     acc,
               int main_iter_step, // Parámetro para la iteración actual del bucle main
               int debug_start_iter) // Parámetro para saber cuándo empezar el debug detallado
{
    int n = static_cast<int>(N.size()) - 1;
    static std::vector<std::vector<int>> adj_cache;
    if (adj_cache.empty()) adj_cache.resize(n + 1);

    bool again = true;

    if (main_iter_step >= debug_start_iter) {
        std::cout << "[DETAILED_DEBUG T" << main_iter_step << " normalize] Inicio. A.size()=" << A.size() << std::endl;
    }

    int normalize_internal_loop_count = 0;
    while (again) {
        normalize_internal_loop_count++;
        again = false;
        if (main_iter_step >= debug_start_iter) {
            std::cout << "[DETAILED_DEBUG T" << main_iter_step << " normalize] Ciclo interno NORM #" << normalize_internal_loop_count 
                      << ". A.size()=" << A.size() << std::endl;
        }

        buildAdj(A, adj_cache);
        label(N, adj_cache, A, acc);

        for (size_t i = 0; i < A.size(); ++i) {        
            const Arc& current_arc = A[i];    

            /* Action 1 (strong-minus) — **cortar TODA strong-minus**  */
            if (current_arc.strong && !current_arc.plus) {   // ← condición corregida (sin acc[u] < 0)
                if (main_iter_step >= debug_start_iter) {
                    std::cout << "[DETAILED_DEBUG T" << main_iter_step << " normalize] Action 1 para A[" << i
                            << "]:(" << current_arc.u << "," << current_arc.v << ") strong=" << current_arc.strong
                            << ", plus=" << current_arc.plus
                            << ". acc_branch=" << acc[current_arc.u]     // branch en el extremo !plus
                            << std::endl;
                    std::cout << "[DETAILED_DEBUG T" << main_iter_step << " normalize] -> Eliminando A[" << i
                            << "]:(" << current_arc.u << "," << current_arc.v << ")" << std::endl;
                }

                int q_node = current_arc.u;            // Nodo para el dummy (extremo que apunta a x0)
                A.erase(A.begin() + i);
                addDummy(A, q_node);

                again = true;
                break;                                 // El while(again) continuará
            }

            /* Action 2 (strong-plus fuera del root) */
            if (current_arc.strong && current_arc.plus && current_arc.u != 0) {
                if (main_iter_step >= debug_start_iter) {
                    std::cout << "[DETAILED_DEBUG T" << main_iter_step << " normalize] Action 2 para A[" << i
                            << "]:(" << current_arc.u << "," << current_arc.v << ") u=" << current_arc.u
                            << ", strong=" << current_arc.strong << ", plus=" << current_arc.plus
                            << ". acc_branch=" << acc[current_arc.v]     // branch en el extremo plus
                            << std::endl;
                    std::cout << "[DETAILED_DEBUG T" << main_iter_step << " normalize] -> Eliminando A[" << i
                            << "]:(" << current_arc.u << "," << current_arc.v << ")" << std::endl;
                }

                int v_node = current_arc.v;            // Nodo para el dummy
                A.erase(A.begin() + i);
                addDummy(A, v_node);

                again = true;
                break;                                 // El while(again) continuará
            }
        }

        // Si 'again' es true, el 'while' continuará. Si es false, saldrá.
    }
    if (main_iter_step >= debug_start_iter) {
        std::cout << "[DETAILED_DEBUG T" << main_iter_step << " normalize] Fin. A.size()=" << A.size() 
                  << ". Total ciclos internos NORM: " << normalize_internal_loop_count << std::endl;
    }
}

/* ---------- *** MOD 2: construir Y completo (sub-rama) ---------- */
vector<char> buildY(const vector<Arc>&A, const vector<Node>&N)
{
    int n=N.size()-1;
    vector<vector<int>> adj(n+1);
    for(const auto& a:A){
        adj[a.u].push_back(a.v);
        adj[a.v].push_back(a.u); 
    }

    /* profundidad desde root */
    vector<int> depth(n+1,-1); queue<int>q; depth[0]=0; q.push(0);
    while(!q.empty()){
        int u=q.front(); q.pop();
        for(int v:adj[u]) if(depth[v]==-1){ depth[v]=depth[u]+1; q.push(v); }
    }

    vector<char> inY(N.size(),0); queue<int> yq;
    for(const auto& a:A) if(a.plus && a.strong){
        inY[a.v]=1; yq.push(a.v);               // raíz de la sub-rama
    }
    while(!yq.empty()){
        int u=yq.front(); yq.pop();
        for(int w:adj[u])
            if(depth[w]>depth[u] && !inY[w]){   // sólo “hacia abajo”
                inY[w]=1; yq.push(w);
            }
    }
    return inY;
}

/* ---------- *** MOD 3: cierre = Σ pesos de nodos en Y ---------- */
double closure(const vector<Arc>&A, const vector<Node>&N){
    vector<char> inY = buildY(A,N);
    double s=0;
    for(size_t i=1;i<N.size();++i) if(inY[i]) s+=N[i].w;
    return s;
}
/* ---------- NUEVO: contar bloques del cierre ---------- */
size_t closureBlocks(const vector<Arc>& A, const vector<Node>& N){
    vector<char> inY = buildY(A, N);
    size_t n = 0;
    for(size_t i = 1; i < N.size(); ++i)
        if(inY[i]) ++n;
    return n;
}


/* respaldo del vector A */
static std::vector<Arc> A_bak;
inline void backupA(const std::vector<Arc>& A){ A_bak = A; }
inline void restoreA(std::vector<Arc>& A){ A = A_bak; }


void decideArc(int u, int v,
               vector<Arc>& A,
               const vector<Node>& N,
               vector<double>&acc,
               int step, int DBG_START,
               unordered_set<uint64_t>& procesadas)
{
    /* ———  evitar repetir el mismo arco ——— */
    uint64_t key = arcKey(u, v);
    if (procesadas.count(key)) return;


    /* ———  mantener la mejor alternativa ——— */
    backupA(A);
    double bestClosure = -numeric_limits<double>::infinity();
    vector<Arc> bestA;

    auto try_variant = [&](bool cutU)                
    {                                                
        restoreA(A);

        if (cutU)  eliminaDummy(A, u);
        else       eliminaDummy(A, v);

        /* añadir arco real (u ➜ v) */
        A.push_back({u, v, true, false});            // plus = provisional; se corrige en normalize
        normalize(N, A, acc, step, DBG_START);

        double cl = closure(A, N);
        if (cl > bestClosure) { bestClosure = cl; bestA = A; }
    };

    try_variant(true);   // FORZAR
    try_variant(false);  // PODAR

    /* ———  adoptar el árbol ganador ——— */
    A.swap(bestA);

    procesadas.insert(key);
}

// ---------- NUEVO: lectura de precedencias (u,v) ----------
vector<pair<int,int>> readPrec(const string& path,
                               const unordered_map<int,int>& id2idx)
{
    ifstream f(path);
    if(!f){ cerr << "No pude abrir " << path << '\n'; exit(1); }

    vector<pair<int,int>> E;
    string line;
    while (getline(f,line)){
        if(line.empty() || line[0]=='%') continue;

        vector<int> tok; int x;
        stringstream ss(line);
        while (ss >> x) tok.push_back(x);
        if (tok.size() < 2) continue;

        int b = tok[0];
        /* ----------- ESTA es la condición correcta ----------- */
        size_t first = (tok.size() - 2 == (size_t)tok[1]) ? 2 : 1;

        for (size_t i = first; i < tok.size(); ++i){
            auto itB = id2idx.find(b), itP = id2idx.find(tok[i]);
            if (itB != id2idx.end() && itP != id2idx.end())
                E.emplace_back(itB->second, itP->second);
        }
    }
    return E;                 // newman1 => 3 922 arcos exactamente
}


/* ====================  main  ==================== */
int main(int argc,char*argv[]){
    if (argc < 2) { cerr << "Uso: " << argv[0] << " bloques.csv [arcos.csv]\n"; return 1;}
    // En main()
    const int DETAILED_DEBUG_START_ITERATION = 99999; //numero de iteraciones
    int step = 1;

    unordered_map<int,int> id2idx;
    unordered_map<Key3,int,Key3Hash> coord2idx;
    vector<Node> N    = readBlocks(argv[1], coord2idx, id2idx);
    std::cout << "[INFO] Total nodos (N.size()): " << N.size() << std::endl;
    DBG("Nodos totales: " << N.size()-1);
    //version original:
  //vector<pair<int,int>> feas = arcs1_5(N, idx);

    //version que lee precedencias de archivo:
        vector<pair<int,int>> feas;
        if (argc >= 3) {
            feas = readPrec  (argv[2], id2idx);            // usa CSV explícito
            cout << "[INFO] Precedencias leidas: "
                << feas.size() << '\n';
        } else {
            feas = arcs1_5(N, coord2idx);                // patrón 1-5 si no hay archivo
        }

    DBG("Aristas factibles iniciales: " << feas.size());
    /* *** MOD 4: ordenar aristas factibles para reproducir ejemplo */
    sort(feas.begin(),feas.end(),
         [&](auto&a,auto&b){return N[a.first].w > N[b.first].w;});

    /* ---------- T0 ---------- */
    vector<Arc> A; // Usando 'A' como en tu código original
    for(size_t i=1;i<N.size();++i) A.push_back({0,(int)i,true,false});
    vector<double> acc; // Usando 'acc' como en tu código original

    // Definir un 'step' conceptual para la fase T0 si es necesario para los logs dentro de normalize.
    // O simplemente pasar valores que desactiven los logs detallados para esta llamada.
    int step_for_T0 = 0; // O cualquier valor que sepas que es < DETAILED_DEBUG_START_ITERATION

    // Llamada actualizada:
    normalize(N, A, acc, step_for_T0, DETAILED_DEBUG_START_ITERATION);

    DBG("ANTES DE T0: feas.size()=" << feas.size());
    DBG("ANTES DE T0: A.size()=" << A.size());
    for(int i = 0; i < min(5, (int)A.size()); ++i)
        DBG("  A["<<i<<"] = ("<<A[i].u<<","<<A[i].v<<")");

    /* util para encontrar (x0,·) 
    auto posDummy=[&](int node){
        for(size_t k=0;k<A.size();++k)
            if(A[k].u==0 && A[k].v==node) return (int)k;
        return -1;
    };
    */
    DBG("T0 creado con " << A.size() << " aristas y closure=" 
        << closure(A,N));
/* ===========  BUCLE PRINCIPAL L-G — cierre máximo  =========== */
    while (true)
    {
        if (step >= DETAILED_DEBUG_START_ITERATION)
            std::cout << "\n[DETAILED_DEBUG T" << step
                    << "] === Inicio iteración principal ===  A.size()=" << A.size() << '\n';

        /* 1) Etiqueta & normaliza el árbol Tk --------------------- */
        normalize(N, A, acc, step, DETAILED_DEBUG_START_ITERATION);

        /* 2) Construye Y  (nodos unidos a la raíz x0 sólo por strong-plus) */
        vector<char> inY = buildY(A, N);                     // ya definida

        /* 3) Busca el PRIMER arco del grafo (feas) que salga de Y hacia X\Y */
        int viol_u = -1, viol_v = -1;
        for (auto [u,v] : feas) {
            if (inY[u] && !inY[v] && !procesadas.count(arcKey(u,v))) {
                viol_u = u; viol_v = v;
                break;                                       // basta uno
            }
        }

        /* 4) ¿Termina?  — no queda ningún arco Y→X\Y en TODO el grafo */
        if (viol_u == -1) {
            if (step >= DETAILED_DEBUG_START_ITERATION)
                std::cout << "[DETAILED_DEBUG T" << step
                        << "]  No hay violaciones → cierre máximo alcanzado.\n";
            break;                                           // FIN
        }

        /* 5) Trata el arco violador con lógica FORZAR / PODAR ------- */
        decideArc(viol_u, viol_v,
                A, N, acc,
                step, DETAILED_DEBUG_START_ITERATION,
                procesadas);

        /* 6) Estadística y control ---------------------------------- */
        if (step < DETAILED_DEBUG_START_ITERATION)
            std::cout << "T" << step << "  cierre = " << closure(A, N) << '\n';

        ++step;
        if (step > 100000) { std::cerr << "Guard-step triggered\n"; break; }
    }
    double valor = closure(A, N);
    size_t bloques = closureBlocks(A, N);

    vector<char> inY_final = buildY(A, N);

    std::ofstream fout("pit_final.csv");
    fout << "id,X,Y,Z,Economic\n";                 // cabecera
    for (size_t i = 1; i < N.size(); ++i) {
        if (inY_final[i]) {                        // forma parte del pit
            fout << i              << ','
                 << N[i].x         << ','
                 << N[i].y         << ','
                 << N[i].z         << ','
                 << N[i].w         << '\n';
        }
    }
    fout.close();
    std::cout << "Archivo \"pit_final.csv\" escrito con "
              << bloques << " bloques del pit final.\n"; 

    DBG("Resultado final Maximo cierre = " << closure(A,N));
    std::cout << "\nMaximo cierre = " << closure(A,N) << "\n";
    std::cout << "Bloques en pit cierre maximo = " << bloques << '\n';
    std::cout << std::llround( closure(A, N) ) << '\n';
    return 0;
}