#include <pagmo/algorithms/sga.hpp>
#include <pagmo/algorithms/pso.hpp>
#include <pagmo/population.hpp>
#include <iostream>
#include <chrono>
#include <fstream>
#include <vector>
#include <functional>
#include <numeric>
#include <iomanip>

#include "loader.hpp"
#include "fitness.hpp"

using namespace std;

struct PitProblem {
    BlockArray blocks;
    ArcAdj     pred;

    // Devuelve fitness = -NPV para minimización
    vector<double> fitness(const vector<double>& x) const {
        auto perm = perm_from_continuous(x);
        return { -pit_npv(perm, blocks, pred) };
    }
    pair<vector<double>, vector<double>> get_bounds() const {
        size_t n = blocks.size() - 1;
        return { vector<double>(n, 0.0), vector<double>(n, 1.0) };
    }
    bool has_gradient() const { return false; }
};

void write_results(const string &tag,
                   const pagmo::population &pop,
                   double runtime_s,
                   const BlockArray &B,
                   const ArcAdj &pred) {
    // Decodificar permutación
    auto bestx = pop.champion_x();
    vector<int> perm = perm_from_continuous(bestx);

    // 1) Clausura completa usando todas las precedencias
    size_t n = B.size();
    vector<char> pit(n, 0);
    function<void(int)> markFull = [&](int u) {
        if (pit[u]) return;
        pit[u] = 1;
        auto it = pred.find(u);
        if (it != pred.end()) {
            for (int p : it->second) markFull(p);
        }
    };
    // incluir sólo bloques con valor positivo
    for (int v : perm) {
        if (B[v].value > 0) markFull(v);
    }

    // 2) Calcular métricas
    int    bad    = 0;
    double benefit= 0;
    double oreT   = 0, wasteT = 0;
    double cuT    = 0, auT    = 0;
    for (size_t v = 1; v < n; ++v) {
        if (pit[v]) {
            // verificar precedencias
            auto it = pred.find(v);
            if (it != pred.end()) {
                for (int p : it->second) if (!pit[p]) { bad++; break; }
            }
            // acumular métricas
            benefit += B[v].value;
            oreT     += B[v].tonnage;
            cuT      += B[v].tonnage * B[v].cu_grade;
            auT      += B[v].tonnage * B[v].au_grade;
        } else {
            wasteT   += B[v].tonnage;
        }
    }
    double strip = oreT > 0 ? (wasteT / oreT) : 0;
    double gCu   = oreT > 0 ? (cuT / oreT)  : 0;
    double gAu   = oreT > 0 ? (auT / oreT)  : 0;
    int cntPit   = accumulate(pit.begin(), pit.end(), 0);

    // 3) Mostrar en pantalla
    cout << "--- " << tag << " results ---\n"
         << "Precedencias violadas: " << bad << "\n"
         << fixed << setprecision(3)
         << "NPV: "    << benefit << "\n"
         << "Rec Cu: " << cuT     << " | Au: " << auT << "\n"
         << "Strip: "  << strip   << "\n"
         << "Ley: "    << gCu*100 << "% Cu, " << gAu << " g/t Au\n"
         << "Tiempo: " << runtime_s << " s\n"
         << "Bloques en pit: " << cntPit << "\n\n";

    // 4) Guardar archivos de salida
    string out_dir = "C:/Users/Giova/OneDrive/Escritorio/Mining/output/GA_output/";  // ajusta si quieres otro destino
    {
        ofstream ids(out_dir + "pit_ids_" + tag + ".csv");
        for (size_t v = 1; v < n; ++v) if (pit[v]) ids << B[v].id << "\n";
    }
    {
        ofstream fout(out_dir + "pf_resultados_" + tag + ".txt", ios::app);
        auto now = chrono::system_clock::to_time_t(chrono::system_clock::now());
        fout << "----- " << ctime(&now)
             << fixed << setprecision(3)
             << "NPV: "    << benefit << "\n"
             << "Rec Cu: " << cuT     << " | Au: " << auT << "\n"
             << "Strip: "  << strip   << "\n"
             << "Ley: "    << gCu*100 << "% Cu, " << gAu << " g/t Au\n"
             << "Tiempo: " << runtime_s << " s\n"
             << "Bloques en pit: " << cntPit << "\n\n";
    }
}

int main(int argc, char **argv) {
    if (argc != 2) {
        cerr << "Uso: ./pit_opt C:/Users/Giova/OneDrive/Escritorio/Mining/data/converted/Pso_Ga/data.json\n";
        return 1;
    }
    auto blocks = load_blocks_json(argv[1]);
    auto pred   = load_arcs_json  (argv[1]);

    pagmo::problem prob{ PitProblem{blocks, pred} };

    // --- PSO con población independiente ---
    pagmo::population pop_pso(prob, 60);
    auto t2 = chrono::steady_clock::now();
    pagmo::algorithm pso( pagmo::pso(500, 0.72, 1.2, 1.2) );
    pop_pso = pso.evolve(pop_pso);
    auto t3 = chrono::steady_clock::now();
    write_results("pso", pop_pso, chrono::duration<double>(t3-t2).count(), blocks, pred);

    return 0;
}
