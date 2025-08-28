/*  Lerchs & Grossmann clásico  –  versión clara y didáctica
    Entrada :
      nodes.csv  → id,x,y,z,valor,tonelaje,gradeCu,gradeAu
      arcs.csv   → from,to   (from encima de to)
    Salida  :
      Beneficio económico total, recuperación Cu/Au, strip ratio, tiempo,
      log de resultados y lista de IDs para 3-D
*/
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <queue>
#include <chrono>
#include <iomanip>
#include <functional>

using namespace std;

/* ---------- Estructuras de datos ---------- */
struct Block {
    int    id;
    double x, y, z;     // coordenadas 3D
    double value;       // beneficio neto
    double tonnage;     // toneladas
    double gradeCu;     // fracción de Cu
    double gradeAu;     // g/t de Au
};

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // 1) Leer nodes.csv
    ifstream fn("C:/Users/Giova/OneDrive/Escritorio/Mining/data/converted/pseudoflow_LG/nodes.csv");
    if(!fn){ cerr<<"No se encuentra nodes.csv\n"; return 1; }
    string line;
    getline(fn, line);  // cabecera

    vector<Block> blk;
    while(getline(fn, line)){
        if(line.empty()) continue;
        stringstream ss(line);
        Block b; string s;
        getline(ss, s, ','); b.id      = stoi(s);
        getline(ss, s, ','); b.x       = stod(s);
        getline(ss, s, ','); b.y       = stod(s);
        getline(ss, s, ','); b.z       = stod(s);
        getline(ss, s, ','); b.value   = stod(s);
        getline(ss, s, ','); b.tonnage = stod(s);
        getline(ss, s, ','); b.gradeCu = stod(s);
        getline(ss, s, ','); b.gradeAu = stod(s);
        blk.push_back(b);
    }
    fn.close();
    int n = (int)blk.size();

    // DEBUG
    cout<<"[DEBUG] Bloques cargados: "<<n<<"\n";
    int cntPos=0; for(int i=0;i<n;++i) if(blk[i].value>0) ++cntPos;
    cout<<"[DEBUG] Bloques con valor>0: "<<cntPos<<"\n";

    // 2) Leer arcs.csv y construir árbol de precedencias
    ifstream fa("C:/Users/Giova/OneDrive/Escritorio/Mining/data/converted/pseudoflow_LG/arcs.csv");
    if(!fa){ cerr<<"No se encuentra arcs.csv\n"; return 1; }
    getline(fa, line);  // cabecera

    vector<vector<int>> below(n);          // grafo de extracción: nodo→hijos
    unordered_map<int,int> idx; idx.reserve(n);
    for(int i=0;i<n;++i) idx[blk[i].id]=i;

    int arcLines=0;
    while(getline(fa,line)){
        if(line.empty()) continue;
        stringstream ss(line);
        string su, sv;
        getline(ss, su, ',');
        getline(ss, sv, ',');
        int uId = stoi(su), vId = stoi(sv);
        if(!idx.count(uId) || !idx.count(vId)) continue;
        int u = idx[uId], v = idx[vId];
        // en Lerchs & Grossmann, u→v significa que v cuelga de u (u padre)
        below[u].push_back(v);
        ++arcLines;
    }
    fa.close();

    long totalArcs=0;
    for(int i=0;i<n;++i) totalArcs += below[i].size();
    cout<<"[DEBUG] Líneas de arcs.csv: "<<arcLines<<"\n";
    cout<<"[DEBUG] Arcos en memoria: "<<totalArcs<<"\n";

    // 3) Construir árbol BFS desde raíz ficticia ρ
    //    Para simular ρ, incluimos todos los bloques con value>0 como hijos de ρ
    int root = n;  // índice ficticio
    vector<vector<int>> treeChildren(n+1);
    vector<char> inPit(n);
    // Inicializar inPit con valor>0
    for(int i=0;i<n;++i) inPit[i] = (blk[i].value > 0);
    // Añadir ρ→i para cada i con valor>0
    for(int i=0;i<n;++i) if(inPit[i]) treeChildren[root].push_back(i);

    // BFS para que todos queden conectados
    queue<int> q0;
    q0.push(root);
    while(!q0.empty()){
        int u = q0.front(); q0.pop();
        for(int v: treeChildren[u]){
            for(int w: below[v]){
                // si no estaba ya añadido bajo v, lo colgamos
                if(!inPit[w]){
                    inPit[w] = 1;
                    treeChildren[v].push_back(w);
                    q0.push(w);
                }
            }
        }
    }
    // Restaurar inPit inicial (solo blocks value>0)
    for(int i=0;i<n;++i) inPit[i] = (blk[i].value > 0);

    // 4) Fase MERGE (cierre) – expandir inPit
    bool changed=true;
    while(changed){
        changed=false;
        for(int v=0;v<n;++v) if(inPit[v]){
            // incluir todos los antecesores en el árbol BFS
            for(int u: treeChildren[root]) if(u==v){
                // ρ→v directo, nada que hacer
            }
            for(int u: below[v]){
                if(!inPit[u]){
                    inPit[u]=1;
                    changed=true;
                }
            }
        }
    }

    // 5) Fase PRUNING (poda ramas negativas) – post-orden DFS sobre treeChildren
    vector<double> mass(n+1, 0.0);
    vector<int> post;
    function<void(int)> dfs = [&](int u){
        for(int c: treeChildren[u]) dfs(c);
        post.push_back(u);
    };
    dfs(root);

    for(int u: post){
        if(u==root) continue;
        if(!inPit[u]) continue;
        double m = blk[u].value;
        for(int c: treeChildren[u]) if(inPit[c])
            m += mass[c];
        if(m<0){
            function<void(int)> drop = [&](int x){
                inPit[x]=0;
                for(int cc: treeChildren[x]) if(inPit[cc]) drop(cc);
            };
            drop(u);
        } else {
            mass[u] = m;
        }
    }

    // 6) Métricas y tiempo
    auto toc = chrono::steady_clock::now();
    double secs = chrono::duration<double>(toc - chrono::steady_clock::now()).count();
    double benefit=0, oreT=0, wasteT=0, cuT=0, auT=0;
    for(int i=0;i<n;++i) if(inPit[i]){
        benefit += blk[i].value;
        if(blk[i].value>=0){
            oreT += blk[i].tonnage;
            cuT  += blk[i].tonnage*blk[i].gradeCu;
            auT  += blk[i].tonnage*blk[i].gradeAu;
        } else {
            wasteT += blk[i].tonnage;
        }
    }
    double strip = oreT>0 ? wasteT/oreT : 0;
    double gCu = oreT>0 ? cuT/oreT : 0;
    double gAu = oreT>0 ? auT/oreT : 0;

    // 7) Imprimir
    cout<<fixed<<setprecision(3);
    cout<<"Beneficio economico: "<<benefit<<"\n";
    cout<<"Recuperacion Cu (t): "<<cuT<<" | Au (g): "<<auT<<"\n";
    cout<<"Strip ratio: "<<strip<<"\n";
    cout<<"Ley media: "<<gCu*100<<"% Cu, "<<gAu<<" g/t Au\n";
    cout<<"Tiempo: "<<secs<<" s\n";

    // 8) Guardar log
    ofstream fout("C:/Users/Giova/OneDrive/Escritorio/Mining/output/LG_output/lg_resultados.txt", ios::app);
    auto now = chrono::system_clock::to_time_t(chrono::system_clock::now());
    fout<<"----- "<<ctime(&now);
    fout<<fixed<<setprecision(3)
        <<"Beneficio: "<<benefit<<"\n"
        <<"Rec Cu: "<<cuT<<" | Au: "<<auT<<"\n"
        <<"Strip: "<<strip<<"\n"
        <<"Ley: "<<gCu*100<<"% Cu, "<<gAu<<" g/t Au\n"
        <<"Tiempo: "<<secs<<" s\n\n";
    fout.close();

    // 9) Guardar IDs para 3D
    ofstream plist("C:/Users/Giova/OneDrive/Escritorio/Mining/output/LG_output/pit_ids.csv");
    for(int i=0;i<n;++i)
        if(inPit[i]) plist<<blk[i].id<<"\n";
    plist.close();

    return 0;
}
