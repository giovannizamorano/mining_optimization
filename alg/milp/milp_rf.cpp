#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <mlpack/core.hpp>
#include <mlpack/methods/random_forest/random_forest.hpp>

using namespace mlpack;
using namespace arma;

// Carga de datos desde MARVIN_BM.csv
// Columnas esperadas por fila:
// 0: X, 1: Y, 2: Z, 3: Rock Type, 4: Density, 5: Economic, 6: CU, 7: AU, 8: Perfin
void cargarDatos(const std::string& archivo,
                 arma::mat& caracteristicas,
                 arma::Row<size_t>& etiquetas,
                 arma::vec& valores,
                 arma::vec& costos,
                 arma::vec& volumenes)
{
    std::ifstream in(archivo);
    std::string linea;

    std::vector<size_t> etiquetasDatos;
    size_t lineCount = 0;

    while (std::getline(in, linea))
    {
        ++lineCount;
        std::stringstream ss(linea);
        std::vector<double> fila;
        std::string campo;

        // Separar por comas
        while (std::getline(ss, campo, ','))
            fila.push_back(std::stod(campo));

        if (fila.size() < 9)
        {
            std::cerr << "Advertencia: línea " << lineCount
                      << " ignorada (menos de 9 campos).\n";
            continue;
        }

        // 1) ECONOMIC (col 5) → función objetivo
        valores.insert_rows(valores.n_elem, 1);
        valores(valores.n_elem - 1) = fila[5];

        // 2) DENSITY (col 4) → volumen/tonelaje
        volumenes.insert_rows(volumenes.n_elem, 1);
        volumenes(volumenes.n_elem - 1) = fila[4];

        // 3) COSTOS: como no hay columna dedicada, lo ponemos igual a density
        costos.insert_rows(costos.n_elem, 1);
        costos(costos.n_elem - 1) = fila[4];

        // 4) LABELS: bloque valioso si Economic > 0
        etiquetasDatos.push_back(fila[5] > 0 ? 1 : 0);

        // 5) CARACTERÍSTICAS: [RockType, Density, CU, AU, Perfin]
        arma::vec bloque(5);
        bloque(0) = fila[3];
        bloque(1) = fila[4];
        bloque(2) = fila[6];
        bloque(3) = fila[7];
        bloque(4) = fila[8];
        caracteristicas.insert_cols(caracteristicas.n_cols, bloque);
    }

    // Convertir vector de etiquetas a arma::Row<size_t>
    etiquetas = arma::Row<size_t>(etiquetasDatos);
}

// Función principal
int main()
{
    arma::mat caracteristicas;
    arma::Row<size_t> etiquetas;
    arma::vec valores, costos, volumenes;

    std::string archivo = "C:/Users/Giova/OneDrive/Escritorio/Mining/data/MARVIN_BM.csv";
    cargarDatos(archivo, caracteristicas, etiquetas, valores, costos, volumenes);

    std::cout << "Datos cargados: " << caracteristicas.n_cols
              << " bloques con " << caracteristicas.n_rows
              << " variables cada uno.\n";

    // Entrenar Random Forest (2 clases: valioso / no valioso)
    RandomForest<> modelo;
    modelo.Train(caracteristicas, etiquetas, 2);

    // Guía: seleccionar solo bloques con predicción = 1
    std::vector<size_t> seleccion;
    arma::mat featSel;
    arma::vec valSel, costSel, volSel;

    for (size_t i = 0; i < caracteristicas.n_cols; ++i)
    {
        arma::vec col = caracteristicas.col(i);
        arma::Row<size_t> pred;
        modelo.Classify(col, pred);

        if (pred(0) == 1)
        {
            seleccion.push_back(i);
            featSel.insert_cols(featSel.n_cols, col);
            valSel.insert_rows(valSel.n_elem, 1);
            costSel.insert_rows(costSel.n_elem, 1);
            volSel.insert_rows(volSel.n_elem, 1);
            valSel(valSel.n_elem - 1) = valores(i);
            costSel(costSel.n_elem - 1) = costos(i);
            volSel(volSel.n_elem - 1) = volumenes(i);
        }
    }

    // Escribir pit_final.csv directamente (opcional)
    std::ofstream out("pit_final.csv");
    out << "idx,Economic,Cost,Density\n";
    for (size_t j = 0; j < seleccion.size(); ++j)
    {
        out << seleccion[j] << ","
            << valSel(j)       << ","
            << costSel(j)      << ","
            << volSel(j)       << "\n";
    }
    out.close();

    std::cout << "Selección final: " << seleccion.size()
              << " bloques. Archivo pit_final.csv generado.\n";

    return 0;
}
