#include <string>
#include <iostream>
#include <fstream>
#include "include/json.hpp"
#include <vector>
#include <cmath>
#include <limits>

using namespace std;

// double f_en_tramo(double x0, double y0, double x1, double y1, sdouble x){
    
//     double a = (((y1-y0)/(x1-x0))*(x-x0)) + y0;

//     return a;
// };
std::vector<double> f_en_tramo(double x0, double y0, double x1, double y1, const std::vector<double>& x) {
    std::vector<double> prediccion;
    
    double pendiente = (y1 - y0) / (x1 - x0);
    
    for (double xi : x) {
        prediccion.push_back(pendiente * (xi - x0) + y0);
    }
    
    return prediccion;
}

// double estimar_error_y(sol,x,y){

//     int i = 0;
//     double error = 0;

//     while (i < (sol.size()-1)){

//         vector<double> sub_x, sub_y;
//         tie(sub_x, sub_y) = subconjunto(x, y, sol[i].first, sol[i+1].first);
//         vector<double> prediccion = f_en_tramo(sol[i].first, sol[i].second, sol[i+1].first, sol[i+1].second, sub_x);
//         error = error + calcular_error(prediccion, sub_y);
//         i = i + 1;
//     }

//     return error;
    
// };

// Función para estimar el error total de la recta
double estimar_error_y(const std::vector<std::pair<double, double>>& sol,std::vector<double>& x, std::vector<double>& y) {
    double error = 0.0;
    
    for (size_t i = 0; i < sol.size() - 1; ++i) {
        std::vector<double> sub_x, sub_y;
        std::tie(sub_x, sub_y) = subconjunto(x, y, sol[i].first, sol[i + 1].first);
        std::vector<double> sub_x_np(sub_x.begin(), sub_x.end());
        std::vector<double> prediccion = f_en_tramo(sol[i].first, sol[i].second, sol[i + 1].first, sol[i + 1].second, sub_x_np);
        error += calcular_error(prediccion, sub_y);
    }
    
    return error;
}

double calcular_error(std::vector<double>& vector1, std::vector<double>& vector2) {
    if (vector1.size() != vector2.size()) { 
        std::cerr << "Error: los vectores tienen diferentes longitudes." << std::endl;
        return -1; 
    }

    double error = 0.0;
    for (size_t i = 0; i < vector1.size(); ++i) {
        error += std::abs(vector1[i] - vector2[i]);
    }
    return error;
}

std::pair<std::vector<double>, std::vector<double>> subconjunto(std::vector<double>& x, std::vector<double>& y, double x0, double x1) {
    std::vector<double> sub_X;
    std::vector<double> sub_Y;
    size_t indice_inferior = 0;
    
    // Generamos subconjunto de x entre x0 y x1
    for (size_t i = 0; i < x.size(); ++i) {
        if (x[i] >= x0 && x[i] <= x1) {
            sub_X.push_back(x[i]);
            if (sub_X.size() == 1) { // Guardamos el índice inferior
            indice_inferior = i;
        }
        }
    }

    // Calculamos el índice superior
    size_t indice_superior = indice_inferior + sub_X.size();
    
    // Generamos subconjunto de y respecto al subconjunto de x
    for (size_t i = indice_inferior; i < indice_superior; ++i) {
        sub_Y.push_back(y[i]);
    }
    
    return std::make_pair(sub_X, sub_Y); // Retornamos ambos subconjuntos
}

// Función para estimar el error y guardar en un mapa memo para evitar recálculos
// std::pair<double, std::unordered_map<std::vector<std::pair<double, double>>, double>> estimar_error_y_pd(
//     const std::vector<std::pair<double, double>>& sol,
//     std::vector<double>& x,
//     std::vector<double>& y,
//     std::unordered_map<std::vector<std::pair<double, double>>, double>& memo) {
    
//     std::vector<std::pair<double, double>> key = sol;
    
//     if (memo.find(key) != memo.end()) {
//         return {memo[key], memo};
//     }
    
//     double error = 0.0;
//     for (size_t i = 0; i < sol.size() - 1; ++i) {
//         std::vector<double> sub_x, sub_y;
//         std::tie(sub_x, sub_y) = subconjunto(x, y, sol[i].first, sol[i + 1].first);
//         std::vector<double> sub_x_np(sub_x.begin(), sub_x.end());
//         std::vector<double> prediccion = f_en_tramo(sol[i].first, sol[i].second, sol[i + 1].first, sol[i + 1].second, sub_x_np);
//         error += calcular_error(prediccion, sub_y);
//     }
    
//     memo[key] = error;
//     return {error, memo};
// }

// std::pair<double, std::vector<std::pair<double, double>>> fuerza_bruta(const std::vector<double>& grid_x, const std::vector<double>& grid_y, const std::vector<double>& x, const std::vector<double>& y, size_t N, std::vector<std::pair<double, double>>& sol_parcial) {
//     if (grid_x.size() < N - sol_parcial.size()) {
//         return {std::numeric_limits<double>::max(), {}};
//     } else if (sol_parcial.size() == N) {
//         double error_actual = estimar_error_y(sol_parcial, x, y);
//         return {error_actual, sol_parcial};
//     } else {
//         std::pair<double, std::vector<std::pair<double, double>>> sol_global = {std::numeric_limits<double>::max(), {}};
//         if (N - sol_parcial.size() == 1) {
//             grid_x = {grid_x.back()};
//         }
//         for (double i : grid_y) {
//             sol_parcial.push_back({grid_x[0], i});
//             auto parcial = fuerza_bruta({grid_x.begin() + 1, grid_x.end()}, grid_y, x, y, N, sol_parcial);
//             if (parcial.first < sol_global.first) {
//                 sol_global = parcial;
//             }
//             sol_parcial.pop_back();
//         }
//         if (sol_parcial.size() > 0) {
//             auto parcial = fuerza_bruta({grid_x.begin() + 1, grid_x.end()}, grid_y, x, y, N, sol_parcial);
//             if (parcial.first < sol_global.first) {
//                 sol_global = parcial;
//             }
//         }
//         return sol_global;
//     }
// }

// Para libreria de JSON.
using namespace nlohmann;

int main(int argc, char** argv) {
    std::string instance_name = "../../data/titanium.json";
    std::cout << "Reading file " << instance_name << std::endl;
    std::ifstream input(instance_name);

    json instance;
    input >> instance;
    input.close();

    int K = instance["n"];
    int m = 6;
    int n = 6;
    int N = 5;

    std::cout << K << std::endl;

    // Aca empieza la magia.

    // Ejemplo para guardar json.
    // Probamos guardando el mismo JSON de instance, pero en otro archivo.
    std::ofstream output("test_output.out");

    ////////////////////FUERZA BRUTA///////////
    // std::vector<std::pair<double, double>> sol_parcial;
    // auto resultado = fuerza_bruta(grid_x, grid_y, x, y, N, sol_parcial);
    // std::cout << "Error mínimo estimado: " << resultado.first << std::endl;
    // std::cout << "Puntos correspondientes: ";
    // for (const auto& punto : resultado.second) {
    //     std::cout << "(" << punto.first << ", " << punto.second << ") ";
    // }
    // std::cout << std::endl;

    output << instance;
    output.close();

    return 0;
}
