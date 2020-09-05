#include <iostream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <string>
#include <chrono> 

using namespace std;
typedef long long T;

struct cluster_vectors{	
	T cantidad_vec, dimensiones;
	vector<vector<T>> vector_cluster;
	vector<T> vec_filas;
	vector<T> vec_dimensiones;

	cluster_vectors(vector<T> vec_filas, vector<T> vec_dimensiones){
		this->vec_dimensiones = vec_dimensiones;
		this->vec_filas = vec_filas;
		for (size_t i = 0; i < vec_filas.size(); i++){
			for (size_t j = 0; j < vec_dimensiones.size(); j++){
				cantidad_vec = vec_filas[i];
				dimensiones = vec_dimensiones[j];
				vector_cluster.resize(cantidad_vec);
				rellenar(i, j,vector_cluster);
				imprimir();
			}
		}
	}

	void rellenar(T i,T j,vector<vector<T>>& vec_vec){
		for (size_t i = 0; i < cantidad_vec; i++){
			vec_vec[i].resize(dimensiones);
			for (size_t j = 0; j < dimensiones; j++){
				vec_vec[i][j] = rand() % 20;
			}
		}
	}

	void distancia_euclides(vector<T> vec_1, vector<T> vec_2){
		T total = 0;
		for (size_t i = 0; i < vec_1.size(); i++){
			total += pow(vec_2[i] - vec_1[i], 2);
		}
		sqrt(total);
	}
	
	void imprimir(){
		cout << "Cantidad de Datos: " << cantidad_vec << endl;
		cout << "Dimension de Vectores: " << dimensiones << endl;
		auto start = chrono::high_resolution_clock::now();
		for (size_t i = 0; i < cantidad_vec; i++){
			for (size_t j = i + 1; j < cantidad_vec; j++){
				distancia_euclides(vector_cluster[i], vector_cluster[j]);
			}
		}
		auto stop = chrono::high_resolution_clock::now();
		auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
		cout << "Tiempo tomado: " << duration.count() << " microsegundos" << endl<<endl;
	}
};


int main(){
	srand(time(NULL));
	vector<T> dimensiones{ 4,6,8,10,12,18,20};
	vector<T> cantidad_vec{10000,15000,20000,25000};
	cluster_vectors task(cantidad_vec, dimensiones);
	return 0;
}