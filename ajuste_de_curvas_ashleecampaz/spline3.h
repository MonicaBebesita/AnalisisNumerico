/**
* @file spline3.h
* @brief Clase para realizar interpolación cúbica por tramos utilizando splines cúbicos.
* 
* Esta clase proporciona una implementación de splines cúbicos para la interpolación
* de puntos en un conjunto de datos.
*/

#ifndef SPLINE3_H
#define SPLINE3_H

#include <cmath>
#include <vector>
#include <iostream>

using std::vector;
using std::cout;
using std::cin;
using std::endl;

/**
* @namespace interpolacion
* @brief Espacio de nombres que contiene las funciones y clases para la interpolación.
*/
namespace interpolacion {
	
	/**
	* @class spline3
	* @brief Clase para realizar interpolación cúbica por tramos.
	* 
	* La clase utiliza splines cúbicos para interpolar valores de datos discretos.
	*/
	class spline3 {
	public:
		/**
		* @brief Constructor que inicializa los puntos de interpolación.
		* @param x Vector de valores de la variable independiente.
		* @param y Vector de valores de la variable dependiente.
		*/
		spline3(vector<double> x, vector<double> y) : x(x), y(y) {
			construir();
		}
		
		/**
		* @brief Verifica si la construcción del spline fue exitosa.
		* @return true si el spline es válido, false en caso contrario.
		*/
		bool es_valido() {
			return valido;
		}
		
		/**
		* @brief Interpola un valor dado utilizando el spline cúbico.
		* @param x_int El valor de la variable independiente a interpolar.
		* @return El valor interpolado, o NAN si el spline no es válido o el valor está fuera de rango.
		*/
		double interpolar(double x_int) {
			int i;
			for (i = 0; i < (int)x.size() && x_int > x[i]; i++);
			
			// Imprimir los límites generales del intervalo
			cout << "Intervalo general: [" << x.front() << ", " << x.back() << "]" << endl;
			
			cout << "Intervalo en el que se encuentra el dato " << x_int << " : " << i << endl;
			
			// Usando los nodos intermedios
			if (i == 0 || i == (int)x.size()) {
				return NAN;
			}
			
			// Imprimir el polinomio para el intervalo específico
			imprimir_polinomios();
			
			return evaluar(i - 1, x_int);
		}
		
		/**
		* @brief Imprime el error relativo y porcentual entre el valor real y el valor interpolado.
		* @param valor_real El valor real del dato.
		* @param valor_interpolado El valor interpolado calculado.
		*/
		void imprimir_error(double valor_real, double valor_interpolado) {
			if (valor_real == 0) {
				cout << "El valor real es 0, el error relativo no se puede calcular." << endl;
				return;
			}
			
			// Calculamos el error absoluto y relativo correctamente
			double error_relativo = fabs((valor_real - valor_interpolado) / fabs(valor_real));
			double error_relativo_porcentual = error_relativo * 100;
			
			cout << "Error Relativo: " << error_relativo << endl;
			cout << "Error Relativo Porcentual: " << error_relativo_porcentual << "%" << endl;
		}
		
	private:
			vector<double> x; ///< Valores de la variable independiente.
			vector<double> y; ///< Valores de la variable dependiente.
			vector<double> f2; ///< Valores de las segundas derivadas.
			bool valido = false; ///< Indica si el spline es válido.
			
			/**
			* @brief Construye el spline cúbico calculando las segundas derivadas en los nodos.
			*/
			void construir() {
				int n = (int)x.size();
				if (x.size() != y.size() || n < 3) {
					valido = false;
					cout << "Interpolación inválida" << endl;
					return;
				}
				vector<double> h(n - 1), alpha(n - 1);
				for (int i = 0; i < n - 1; i++) {
					h[i] = x[i + 1] - x[i];
					if (i > 0) {
						alpha[i] = (3.0 / h[i]) * (y[i + 1] - y[i]) - (3.0 / h[i - 1]) * (y[i] - y[i - 1]);
					}
				}
				
				vector<double> l(n), mu(n), z(n);
				l[0] = 1.0;
				mu[0] = z[0] = 0.0;
				
				for (int i = 1; i < n - 1; i++) {
					l[i] = 2.0 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
					mu[i] = h[i] / l[i];
					z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
				}
				
				l[n - 1] = 1.0;
				z[n - 1] = 0.0;
				f2.resize(n);
				f2[n - 1] = 0.0;
				
				for (int j = n - 2; j >= 0; j--) {
					f2[j] = z[j] - mu[j] * f2[j + 1];
				}
				
				valido = true;
			}
			
			/**
			* @brief Evalúa el spline en un intervalo específico.
			* @param i Índice del intervalo.
			* @param x_int Valor de la variable independiente a interpolar.
			* @return El valor interpolado en el intervalo.
			*/
			double evaluar(int i, double x_int) {
				if (!valido) {
					return NAN;
				}
				double h = x[i + 1] - x[i];
				double A = (f2[i + 1] - f2[i]) / (6.0 * h);
				double B = f2[i] / 2.0;
				double C = (y[i + 1] - y[i]) / h - (f2[i + 1] + 2.0 * f2[i]) * h / 6.0;
				double D = y[i];
				double dx = x_int - x[i];
				return A * dx * dx * dx + B * dx * dx + C * dx + D;
			}
			
			/**
			* @brief Imprime los polinomios cúbicos correspondientes a cada intervalo.
			*/
			void imprimir_polinomios() {
				for (int i = 0; i < (int)x.size() - 1; i++) {
					double h = x[i + 1] - x[i];
					double A = (f2[i + 1] - f2[i]) / (6.0 * h);
					double B = f2[i] / 2.0;
					double C = (y[i + 1] - y[i]) / h - (f2[i + 1] + 2.0 * f2[i]) * h / 6.0;
					double D = y[i];
					
					cout << "Polinomio correspondiente en el intervalo [" << x[i] << ", " << x[i + 1] << "]:\n";
					cout << "S(x) = " << A << "(x - " << x[i] << ")^3 + " 
						<< B << "(x - " << x[i] << ")^2 + " 
						<< C << "(x - " << x[i] << ") + " 
						<< D << endl;
				}
			}
	};
}

#endif // SPLINE3_H
