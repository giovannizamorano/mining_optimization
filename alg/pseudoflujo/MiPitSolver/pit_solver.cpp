// pit_solver.cpp
#include "mineflow.h"
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
using Cap = long long;
constexpr int SCALE = 1000;  // coincide con tu escala original

struct Block {
    int id;
    Cap valInt;
    double ton, gCu, gAu;
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // --- 1. Rutas fijas de tus archivos ---
    const string BASE    = "../../../data/converted/pseudoflow_LG/";
    const string NODES   = BASE + "nodes.csv";
    const string ARCS    = BASE + "arcs.csv";
    const string OUT_DIR = "../../../output/pseudoflujo/";

    // --- 2. Leer bloques y valores ---
    ifstream fn(NODES);
    if (!fn) {
        cerr << "[ERROR] No se puede abrir " << NODES << "\n";
        return 1;
    }
    string line;
    getline(fn, line);  // salta cabecera

    vector<Block> blocks;
    unordered_map<int,int> idx;
    while (getline(fn, line)) {
        stringstream ss(line);
        Block b;
        string s;
        getline(ss, s, ','); b.id     = stoi(s);
        for (int i = 0; i < 3; ++i) getline(ss, s, ','); // X,Y,Z
        getline(ss, s, ','); b.valInt = (Cap)llround(stod(s));
        getline(ss, s, ','); b.ton    = stod(s);
        getline(ss, s, ','); b.gCu    = stod(s);
        getline(ss, s, ','); b.gAu    = stod(s);

        idx[b.id] = blocks.size();
        blocks.push_back(b);
    }
    fn.close();

    int nBlocks = blocks.size();

    // --- 3. Construir precedencias ---
    auto precedence = make_shared<mvd::mineflow::ExplicitPrecedence>(nBlocks);
    ifstream fa(ARCS);
    if (!fa) {
        cerr << "[ERROR] No se puede abrir " << ARCS << "\n";
        return 1;
    }
    getline(fa, line);  // cabecera
    while (getline(fa, line)) {
        if (line.empty()) continue;
        auto p = line.find(',');
        int idU = stoi(line.substr(0, p));
        int idV = stoi(line.substr(p + 1));
        if (idx.count(idU) && idx.count(idV)) {
            precedence->AddPrecedenceConstraint(idx[idU], idx[idV]);
        }
    }
    fa.close();

    // --- 4. Cargar valores en MineFlow ---
    auto values = make_shared<mvd::mineflow::VecBlockValues>(nBlocks);
    for (int i = 0; i < nBlocks; ++i) {
        values->SetBlockValueSI(i, blocks[i].valInt);
    }

    // --- 5. Ejecutar el solver pseudoflow ---
    mvd::mineflow::PseudoSolver solver(precedence, values);
    mvd::mineflow::PseudoSolverSolveInfo info;
    solver.Solve(&info);

    // Extraer pit resultado
    vector<char> pit(nBlocks);
    for (int i = 0; i < nBlocks; ++i) {
        pit[i] = solver.InMinimumCut(i) ? 1 : 0;
    }

    // --- 6. Calcular métricas idénticas a tu código original ---
    Cap    benefitInt = 0;
    double oreT = 0, wasteT = 0, cuT = 0, auT = 0;
    for (int i = 0; i < nBlocks; ++i) if (pit[i]) {
        benefitInt += blocks[i].valInt;
        if (blocks[i].valInt >= 0) {
            oreT   += blocks[i].ton;
            cuT    += blocks[i].ton * blocks[i].gCu;
            auT    += blocks[i].ton * blocks[i].gAu;
        } else {
            wasteT += blocks[i].ton;
        }
    }
    double strip = oreT > 0 ? wasteT / oreT : 0;
    double gCu   = oreT > 0 ? cuT    / oreT : 0;
    double gAu   = oreT > 0 ? auT    / oreT : 0;

    // Contar precedencias violadas
    size_t bad = 0;
    ifstream fa2(ARCS);
    getline(fa2, line);
    while (getline(fa2, line)) {
        if (line.empty()) continue;
        auto p = line.find(',');
        int idU = stoi(line.substr(0, p));
        int idV = stoi(line.substr(p + 1));
        if (idx.count(idU) && idx.count(idV)) {
            bool inU = pit[idx[idU]];
            bool inV = pit[idx[idV]];
            if (inV && !inU) ++bad;
        }
    }
    fa2.close();

    // Tiempo total
    double secs = info.ElapsedSeconds;

    // --- 7. Mostrar sólo tus métricas en pantalla ---
    cout << "Precedencias violadas: " << bad << "\n"
         << fixed << setprecision(3)
         << "Beneficio: " << benefitInt / double(SCALE) << "\n"
         << "Rec Cu: "   << cuT << " | Au: " << auT << "\n"
         << "Strip: "    << strip << "\n"
         << "Ley: "      << gCu * 100 << "% Cu, " << gAu << " g/t Au\n"
         << "Tiempo: "   << secs << " s\n"
         << "Bloques en pit: " << count(pit.begin(), pit.end(), 1)
         << "\n";

    // --- 8. Guardar resultados en tus archivos originales ---
    // 8.1 IDs de pit
    ofstream ids(OUT_DIR + "pit_ids_pf.csv");
    for (int i = 0; i < nBlocks; ++i)
        if (pit[i])
            ids << blocks[i].id << "\n";
    ids.close();

    // 8.2 Resumen en pf_resultados.txt (append)
    ofstream fout(OUT_DIR + "pf_resultados.txt", ios::app);
    auto now = chrono::system_clock::to_time_t(chrono::system_clock::now());
    fout << "----- " << ctime(&now)
         << fixed << setprecision(3)
         << "Beneficio: " << benefitInt / double(SCALE) << "\n"
         << "Rec Cu: "   << cuT << " | Au: " << auT << "\n"
         << "Strip: "    << strip << "\n"
         << "Ley: "      << gCu * 100 << "% Cu, " << gAu << " g/t Au\n"
         << "Tiempo: "   << secs << " s\n\n";
    fout.close();

    return 0;
}
