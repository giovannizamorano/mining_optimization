#include <iostream>
#include <chrono>
#include <fstream>
#include <vector>
#include <algorithm>
#include <random>
#include <functional>
#include <numeric>
#include <iomanip>

#include "loader.hpp"
#include "fitness.hpp"

using namespace std;

// Estructura para representar un individuo en GA
struct Individual {
    vector<int> permutation;
    double fitness;
    
    Individual(size_t n) : permutation(n), fitness(-1e9) {
        // Inicializar con IDs desde 1 hasta n
        iota(permutation.begin(), permutation.end(), 1);
        random_shuffle(permutation.begin(), permutation.end());
    }
};

// Algoritmo genético simple
class SimpleGA {
private:
    vector<Individual> population;
    size_t pop_size;
    double mutation_rate;
    double crossover_rate;
    mt19937 rng;
    
public:
    SimpleGA(size_t pop_size, double crossover_rate, double mutation_rate) 
        : pop_size(pop_size), mutation_rate(mutation_rate), crossover_rate(crossover_rate), rng(random_device{}()) {}
    
    // Evaluar fitness de un individuo
    void evaluate(Individual& ind, const BlockArray& blocks, const ArcAdj& pred) {
        ind.fitness = pit_npv(ind.permutation, blocks, pred);
    }
    
    // Crossover de orden (OX)
    Individual crossover(const Individual& parent1, const Individual& parent2) {
        size_t n = parent1.permutation.size();
        Individual child(n);
        
        // Seleccionar dos puntos de corte
        size_t start = uniform_int_distribution<size_t>(0, n-1)(rng);
        size_t end = uniform_int_distribution<size_t>(start, n-1)(rng);
        
        // Copiar segmento del padre 1
        vector<bool> used(n+1, false);  // +1 porque los IDs van de 1 a n
        for (size_t i = start; i <= end; i++) {
            child.permutation[i] = parent1.permutation[i];
            used[parent1.permutation[i]] = true;
        }
        
        // Llenar el resto con genes del padre 2 en orden
        size_t pos = 0;
        for (size_t i = 0; i < n; i++) {
            if (pos == start) pos = end + 1;
            if (pos >= n) break;
            
            int gene = parent2.permutation[i];
            if (!used[gene]) {
                child.permutation[pos] = gene;
                pos++;
            }
        }
        
        return child;
    }
    
    // Mutación por intercambio
    void mutate(Individual& ind) {
        if (uniform_real_distribution<double>(0, 1)(rng) < mutation_rate) {
            size_t n = ind.permutation.size();
            size_t i = uniform_int_distribution<size_t>(0, n-1)(rng);
            size_t j = uniform_int_distribution<size_t>(0, n-1)(rng);
            swap(ind.permutation[i], ind.permutation[j]);
        }
    }
    
    // Selección por torneo
    Individual tournament_selection(const vector<Individual>& pop, size_t tournament_size = 3) {
        Individual best = pop[uniform_int_distribution<size_t>(0, pop.size()-1)(rng)];
        
        for (size_t i = 1; i < tournament_size; i++) {
            Individual candidate = pop[uniform_int_distribution<size_t>(0, pop.size()-1)(rng)];
            if (candidate.fitness > best.fitness) {
                best = candidate;
            }
        }
        return best;
    }
    
    // Ejecutar evolución
    void evolve(size_t generations, const BlockArray& blocks, const ArcAdj& pred) {
        // Inicializar población
        population.clear();
        for (size_t i = 0; i < pop_size; i++) {
            population.emplace_back(blocks.size() - 1);  // -1 porque blocks[0] no se usa
            evaluate(population.back(), blocks, pred);
        }
        
        cout << "Generación 0 - Mejor fitness: " << get_best().fitness << endl;
        
        // Evolucionar
        for (size_t gen = 1; gen <= generations; gen++) {
            vector<Individual> new_population;
            
            // Elitismo: mantener al mejor
            new_population.push_back(get_best());
            
            // Generar resto de la nueva población
            while (new_population.size() < pop_size) {
                Individual parent1 = tournament_selection(population);
                Individual parent2 = tournament_selection(population);
                
                Individual child1 = crossover(parent1, parent2);
                Individual child2 = crossover(parent2, parent1);
                
                mutate(child1);
                mutate(child2);
                
                evaluate(child1, blocks, pred);
                evaluate(child2, blocks, pred);
                
                new_population.push_back(child1);
                if (new_population.size() < pop_size) {
                    new_population.push_back(child2);
                }
            }
            
            population = new_population;
            
            if (gen % 50 == 0) {
                cout << "Generación " << gen << " - Mejor fitness: " << get_best().fitness << endl;
            }
        }
    }
    
    Individual get_best() const {
        return *max_element(population.begin(), population.end(), 
                          [](const Individual& a, const Individual& b) {
                              return a.fitness < b.fitness;
                          });
    }
};

void write_results(const string &tag,
                   const Individual& best,
                   double runtime_s,
                   const BlockArray &B,
                   const ArcAdj &pred) {
    
    // Usar la permutación directamente
    vector<int>& perm = const_cast<vector<int>&>(best.permutation);

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
    string out_dir = "C:/Users/Giova/OneDrive/Escritorio/Mining/output/GA_output/";
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
        cerr << "Uso: ./pit_simple <ruta_al_archivo_data.json>\n";
        cerr << "Ejemplo: ./pit_simple C:/Users/Giova/OneDrive/Escritorio/Mining/data/converted/Pso_Ga/data.json\n";
        return 1;
    }
    
    cout << "Cargando datos desde: " << argv[1] << endl;
    auto blocks = load_blocks_json(argv[1]);
    auto pred   = load_arcs_json(argv[1]);
    
    cout << "Datos cargados: " << blocks.size()-1 << " bloques" << endl;
    
    // Parámetros del GA
    size_t population_size = 60;
    size_t generations = 200;
    double crossover_rate = 0.9;
    double mutation_rate = 0.02;
    
    cout << "Iniciando Algoritmo Genético..." << endl;
    cout << "Población: " << population_size << ", Generaciones: " << generations << endl;
    
    SimpleGA ga(population_size, crossover_rate, mutation_rate);
    
    auto t0 = chrono::steady_clock::now();
    ga.evolve(generations, blocks, pred);
    auto t1 = chrono::steady_clock::now();
    
    Individual best = ga.get_best();
    double runtime = chrono::duration<double>(t1-t0).count();
    
    write_results("ga", best, runtime, blocks, pred);

    return 0;
}
