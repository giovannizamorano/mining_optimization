// lg2_approach2.cpp  –  Implementación íntegra del Approach 2 (cap. 5)
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

/* ---------- lectura CSV ---------- */
vector<Node> readCSV(const string&file,
                     unordered_map<Key3,int,Key3Hash>&idx,char sep = ','){
    ifstream f(file); if(!f){cerr<<"No file\n";exit(1);}
    string hdr; getline(f,hdr);
    vector<Node>N(1); 
    string ln, tmp;
    while(getline(f, ln)){
        if(ln.empty()) continue;
        stringstream ss(ln);
        vector<string> tok;
        // parte la línea en tantos campos como haya
        while(getline(ss, tmp, sep))
            tok.push_back(tmp);

        // convierte los que te interesan
        double X = stod(tok[0]);            // columna “X”
        double Y = stod(tok[1]);            // columna “Y”
        double Z = stod(tok[2]);            // columna “Z”
        double w = stod(tok[5]);            // columna “Economic”

        int id = N.size();
        N.push_back({ w, X, Y, Z });
        idx[{X,Y,Z}] = id;
    }
    DBG("CSV leído: " << N.size()-1 << " bloques");

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

void computeAccIterative(const vector<Node>& N,
                         const vector<vector<int>>& adj,
                         vector<double>& acc) {
    int n = (int)N.size() - 1;
    acc.assign(n+1, 0.0);
    vector<int> parent(n+1, -1);
    vector<int> order;
    order.reserve(n+1);

    // 1) Construye post-order con stack
    stack<int> st;
    st.push(0);
    while (!st.empty()) {
        int u = st.top(); st.pop();
        order.push_back(u);
        for (int v : adj[u]) {
            if (v == parent[u]) continue;
            parent[v] = u;
            st.push(v);
        }
    }

    // 2) Suma pesos de abajo hacia arriba
    for (auto it = order.rbegin(); it != order.rend(); ++it) {
        int u = *it;
        acc[u] += N[u].w;
        if (parent[u] >= 0)
            acc[parent[u]] += acc[u];
    }
}
/* ---------- aristas 1-5 ---------- */
vector<pair<int,int>> arcs1_5(
        const vector<Node>&N,
        const unordered_map<Key3,int,Key3Hash>&idx)
{
    vector<double>xs,ys,zs;
    for(size_t i=1;i<N.size();++i){
        xs.push_back(N[i].x); ys.push_back(N[i].y); zs.push_back(N[i].z);}
    double dx=minDiff(xs), dy=minDiff(ys), dz=minDiff(zs);
    if (!std::isfinite(dy) || dy > 1e50) {
        DBG("dy no válida (infinito o muy grande), ajustando dy = dx");
        dy = dx;
    }

    DBG("dx=" << dx << "  dy=" << dy << "  dz=" << dz);

    vector<pair<int,int>>E;
    for(size_t id=1; id<N.size(); ++id){
        const auto&b=N[id]; double up=b.z-dz;
        //const auto&b=N[id]; double up=b.z+dz;
        Key3 p[5]={{b.x,b.y,up},{b.x-dx,b.y,up},{b.x+dx,b.y,up},
                   {b.x,b.y-dy,up},{b.x,b.y+dy,up}};
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

/* ---------- etiquetar (Tabla 5.15) ---------- */
void label(const vector<Node>& N,
           vector<Arc>&        A,
           vector<double>&  acc)
{
    int n = N.size()-1;
    acc.assign(n+1,0);
    
    DBG("-- label: N nodos=" << N.size()-1
    << "  A arcos=" << A.size());

    /* --- construir adyacencias sin dirección --- */
    vector<vector<int>> adj(n+1);
    for(auto& a:A){ adj[a.u].push_back(a.v); adj[a.v].push_back(a.u); }

    /* --- BFS para profundidad --- */
    vector<int> depth(n+1,-1);
    queue<int> q; depth[0]=0; q.push(0);
    while(!q.empty()){
        int u=q.front(); q.pop();
        for(int v: adj[u]) if(depth[v]==-1){
            depth[v]=depth[u]+1; q.push(v);
        }
    }


    DBG("Profundidades de nodos:");
    for(int i = 0; i <= n; ++i)
        DBG(" depth[" << i << "] = " << depth[i]);

    /* --- DFS para peso acumulado de cada sub-árbol --- */
    // --- DFS ITERATIVO ---
    // if (in_debug_window) DBG("label (T=" << current_step << "): Iniciando DFS ITERATIVO para acc.");
    int max_dfs_simulated_depth = 0; 
    std::fill(acc.begin(), acc.end(), 0.0); 

    std::stack<pair<int, int>> dfs_iter_stack; 
    std::vector<int> post_order_nodes;   
    post_order_nodes.reserve(n + 1);     

    std::vector<int> dfs_parent_map(n + 1, -1); 

    if (n >= 0) {
        dfs_iter_stack.push({0, 0}); 
    }

    int current_simulated_depth = 0; 

    while(!dfs_iter_stack.empty()){
        int u = dfs_iter_stack.top().first;
        int &child_idx = dfs_iter_stack.top().second; 

        if (child_idx == 0 && u>=0 && (size_t)u < N.size()) { // Primera vez procesando 'u', simular entrada
            current_simulated_depth++;
            if (current_simulated_depth > max_dfs_simulated_depth) {
                max_dfs_simulated_depth = current_simulated_depth;
            }
        }

        bool processed_a_child_this_iteration = false;
        if(u>=0 && (size_t)u < adj.size()){ // Chequeo para adj[u]
            while(child_idx < adj[u].size()){
                int v = adj[u][child_idx];
                child_idx++; 

                if (v == dfs_parent_map[u]) continue;
                
                if(v>=0 && (size_t)v < depth.size() && u>=0 && (size_t)u < depth.size()){ // Chequeos para depth
                    if (depth[v] == depth[u] + 1) { 
                        dfs_parent_map[v] = u;      
                        dfs_iter_stack.push({v, 0});   
                        processed_a_child_this_iteration = true;
                        break; 
                    }
                } else {
                    std::cerr << "[ERROR FATAL DFS ITERATIVO] Indice v o u invalido: v=" << v << " u=" << u
                              << " depth.size()=" << depth.size() 
                              << " N.size()=" << N.size() << std::endl;
                    continue; // Salta este ciclo si hay error
                }
            }
        }


        if(!processed_a_child_this_iteration){
            post_order_nodes.push_back(u);
            if(child_idx == 0 && u>=0 && (size_t)u < N.size()){ // Solo decrementar si se incrementó
                // Esto es un poco impreciso para la profundidad real de la pila simulada si hay ramas vacías,
                // pero max_dfs_simulated_depth debería capturar el máximo anidamiento.
            }
            current_simulated_depth--; // Simular salida (mejorar esta logica de profundidad si es critica)
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
    //    if (in_debug_window) DBG("label (T=" << current_step << "): DFS (iterativo) para acc completado. Maxima profundidad simulada: " << max_dfs_simulated_depth);
    // if (in_debug_window) DBG("label (T=" << current_step << "): DFS (iterativo) para acc completado. Maxima profundidad simulada: " << max_dfs_simulated_depth);
    /* --- asignar plus/minus y strong/weak --- */
    for(auto& a:A){
        a.plus   = (depth[a.v] > depth[a.u]);
        double branch = acc[a.plus ? a.v : a.u];
        /* *** MOD 1: aceptar branch == 0 en arcos plus (variante libro) */
        a.strong = ( (a.plus && branch>=0) || (!a.plus && branch<=0) );
    }
}

/* ---------- normalizar (tree-cut) ---------- */
void normalize(const vector<Node>& N,
               vector<Arc>&        A,
               vector<double>&  acc)
{
    bool again = true;
    while (again) {
        again = false;
        DBG("normalize: inicio ciclo, A.size()=" << A.size());
        label(N, A, acc);                       // (re)etiquetar

        for (size_t i = 0; i < A.size(); ++i) {
            const Arc& a = A[i];

            /* Action 1 (strong-minus) */
            if (a.strong && !a.plus) {
                DBG(" normalize Action1: eliminando arco "
                    << a.u << "->" << a.v);
                int q = a.u;
                A.erase(A.begin() + i);
                bool dup = false;
                for(auto& d:A) if(d.u==0&&d.v==q){dup=true;break;}
                if(!dup) A.push_back({0,q,true,false});
                again = true; break;
            }

            /* Action 2 (strong-plus fuera del root) */
            if (a.strong && a.plus && a.u != 0) {
                DBG(" normalize Action2: eliminando arco "
                    << a.u << "->" << a.v
                    << " e insertando dummy 0->" << a.v);
                int v = a.v;
                A.erase(A.begin() + i);
                bool dup = false;
                for(auto& d:A) if(d.u==0&&d.v==v){dup=true;break;}
                if(!dup) A.push_back({0,v,true,false});
                again = true; break;
            }
        }
        DBG("normalize: fin ciclo, A.size()=" << A.size());
    }
}

/* ---------- *** MOD 2: construir Y completo (sub-rama) ---------- */
vector<char> buildY(const vector<Arc>&A, const vector<Node>&N)
{
    int n=N.size()-1;
    vector<vector<int>> adj(n+1);
    for(const auto& a:A){ adj[a.u].push_back(a.v); adj[a.v].push_back(a.u); }

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

/* ====================  main  ==================== */
int main(int argc,char*argv[]){
    if(argc<2){ cerr<<"Uso: "<<argv[0]<<" archivo.csv\n"; return 1; }

    unordered_map<Key3,int,Key3Hash> idx;
    vector<Node> N   = readCSV(argv[1], idx);
    DBG("Nodos totales: " << N.size()-1);
    vector<pair<int,int>> feas = arcs1_5(N, idx);
    DBG("Aristas factibles iniciales: " << feas.size());
    /* *** MOD 4: ordenar aristas factibles para reproducir ejemplo */
    sort(feas.begin(),feas.end(),
         [&](auto&a,auto&b){return N[a.first].w > N[b.first].w;});

    /* ---------- T0 ---------- */
    vector<Arc> A;
    for(size_t i=1;i<N.size();++i) A.push_back({0,(int)i,true,false});
    vector<double> acc;
    normalize(N,A,acc);
    DBG("ANTES DE T0: feas.size()=" << feas.size());
    DBG("ANTES DE T0: A.size()=" << A.size());
    for(int i = 0; i < min(5, (int)A.size()); ++i)
        DBG("  A["<<i<<"] = ("<<A[i].u<<","<<A[i].v<<")");

    /* util para encontrar (x0,·) */
    auto posDummy=[&](int node){
        for(size_t k=0;k<A.size();++k)
            if(A[k].u==0 && A[k].v==node) return (int)k;
        return -1;
    };
    DBG("T0 creado con " << A.size() << " aristas y closure=" 
        << closure(A,N));
    int step=1;
    while(true){
        DBG("=== Iteración T=" << step << " ===");
        vector<char> inY = buildY(A,N);      // --- Y completo ---
        DBG(" buildY: inY.size()=" << inY.size());
        bool added=false;

        for(auto [u,v] : feas){
            DBG(" Probando arco " << u << "->" << v
                << " inY[u]=" << int(inY[u])
                << " inY[v]=" << int(inY[v]));
            if(inY[u] && !inY[v]){
                DBG("  >> Añadiendo arco " << u << "->" << v);

                /* reglas de dummies */
                int d = posDummy(u);
                if(d!=-1) A.erase(A.begin()+d);
                else {
                        d = posDummy(v);
                        if(d!=-1) A.erase(A.begin()+d);
                    }

                A.push_back({u,v,true,false});   // nueva conexión
                normalize(N,A,acc);

                DBG("  after normalize: A.size()=" << A.size());
                cout << "T" << step << "  cierre = "
                     << closure(A,N) << "\n";
                added=true; break;
            }
        }
        if(!added){DBG(" No quedan arcos factibles, saliendo bucle");break;}
        ++step;                  // no hay más aristas factibles
    }

    DBG("Resultado final Maximo cierre = " << closure(A,N));
    cout << "\nMaximo cierre = " << closure(A,N) << "\n";
    return 0;
}
