#include <string>
#include <iostream>
#include <fstream>
#include "include/json.hpp"
#include <vector>
#include <cmath>
#include <limits>
#include <ctime>

using namespace std;


vector<double> f_en_tramo(double x0, double y0, double x1, double y1, const vector<double>& x) {
    vector<double> prediccion;
    
    double pendiente = (y1 - y0) / (x1 - x0);
    
    for (double xi : x) {
        prediccion.push_back(pendiente * (xi - x0) + y0);
    }
    
    return prediccion;
}



// Función para estimar el error total de la recta

double calcular_error(vector<double>& vector1, vector<double>& vector2) {
    if (vector1.size() != vector2.size()) { 
        cerr << "Error: los vectores tienen diferentes longitudes." << endl;
        return -1; 
    }

    double error = 0.0;
    for (size_t i = 0; i < vector1.size(); ++i) {
        error += abs(vector1[i] - vector2[i]);
    }
    return error;
}

pair<vector<double>, vector<double>> subconjunto(vector<double>& x, vector<double>& y, double x0, double x1) {
    vector<double> sub_X;
    vector<double> sub_Y;
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
    
    return make_pair(sub_X, sub_Y); // Retornamos ambos subconjuntos
}

// Función para estimar el error y guardar en un mapa memo para evitar recálculos
// pair<double, unordered_map<vector<pair<double, double>>, double>> estimar_error_y_pd(
//     const vector<pair<double, double>>& sol,
//     vector<double>& x,
//     vector<double>& y,
//     unordered_map<vector<pair<double, double>>, double>& memo) {
    
//     vector<pair<double, double>> key = sol;
    
//     if (memo.find(key) != memo.end()) {
//         return {memo[key], memo};
//     }
    
//     double error = 0.0;
//     for (size_t i = 0; i < sol.size() - 1; ++i) {
//         vector<double> sub_x, sub_y;
//         tie(sub_x, sub_y) = subconjunto(x, y, sol[i].first, sol[i + 1].first);
//         vector<double> sub_x_np(sub_x.begin(), sub_x.end());
//         vector<double> prediccion = f_en_tramo(sol[i].first, sol[i].second, sol[i + 1].first, sol[i + 1].second, sub_x_np);
//         error += calcular_error(prediccion, sub_y);
//     }
    
//     memo[key] = error;
//     return {error, memo};
// }

double estimar_error_y(const vector<pair<double, double>>& sol,vector<double>& x, vector<double>& y) {
    double error = 0.0;
    
    for (size_t i = 0; i < sol.size() - 1; ++i) {
        vector<double> sub_x, sub_y;
        tie(sub_x, sub_y) = subconjunto(x, y, sol[i].first, sol[i + 1].first);
        vector<double> sub_x_np(sub_x.begin(), sub_x.end());
        vector<double> prediccion = f_en_tramo(sol[i].first, sol[i].second, sol[i + 1].first, sol[i + 1].second, sub_x_np);
        error += calcular_error(prediccion, sub_y);
    }
    
    return error;
}

pair<double, vector<pair<double, double>>> fuerza_bruta(vector<double>& grid_x, vector<double>& grid_y, vector<double>& x, vector<double>& y, size_t N, vector<pair<double, double>>& sol_parcial) {
    if (grid_x.size() < N - sol_parcial.size()) {
        return {numeric_limits<double>::max(), {}};
    } else if (sol_parcial.size() == N) {
        double error_actual = estimar_error_y(sol_parcial, x, y);
        return {error_actual, sol_parcial};
    } else {
        pair<double, vector<pair<double, double>>> sol_global = {numeric_limits<double>::max(), {}};
        if (N - sol_parcial.size() == 1) {
            grid_x = {(grid_x.back())};
        }
        vector<double> grid_x_sliced = {grid_x.begin() + 1, grid_x.end()};
        for (double i : grid_y) {
            sol_parcial.push_back({grid_x[0], i});
            auto parcial = fuerza_bruta(grid_x_sliced, grid_y, x, y, N, sol_parcial);
            if (parcial.first < sol_global.first) {
                sol_global = parcial;
            }
            sol_parcial.pop_back();
        }
        if (sol_parcial.size() > 0) {
            auto parcial = fuerza_bruta(grid_x_sliced, grid_y, x, y, N, sol_parcial);
            if (parcial.first < sol_global.first) {
                sol_global = parcial;
            }
        }
        return sol_global;
    }
}



// Para libreria de JSON.
using namespace nlohmann;

int main(int argc, char** argv) {
    string instance_name = "../../data/titanium.json";
    cout << "Reading file " << instance_name << endl;
    ifstream input(instance_name);

    json instance;
    input >> instance;
    input.close();

    int K = instance["n"];
    int m = 6;
    int n = 6;
    int N = 5;

    //cout << K << endl;

    // Aca empieza la magia.

    // Ejemplo para guardar json.
    // Probamos guardando el mismo JSON de instance, pero en otro archivo.
    ofstream output("test_output.out");

    ////////////////////FUERZA BRUTA///////////
    // vector<pair<double, double>> sol_parcial;
    // auto resultado = fuerza_bruta(grid_x, grid_y, x, y, N, sol_parcial);
    // cout << "Error mínimo estimado: " << resultado.first << endl;
    // cout << "Puntos correspondientes: ";
    // for (const auto& punto : resultado.second) {
    //     cout << "(" << punto.first << ", " << punto.second << ") ";
    // }
    // cout << endl;
    vector<pair<double, double>> sol = { {1.0, 1.0}, {3.0, 3.0}, {6.0, 6.0} };
    vector<double> x = {1.0, 2.0, 3.0, 4.5, 5.0};
    vector<double> y = {1.0, 0.0, 3.0, 3.0, 5.2};
    double error = estimar_error_y(sol, x, y);
    cout << "Error estimado: " << error << endl;


    output << instance;
    output.close();

    return 0;
}
