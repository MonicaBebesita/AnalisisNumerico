#ifndef NEWTHON_H
#define NEWTHON_H

//metodo de interpolacion
#include <vector>
#include <cmath>
using std::vector; 

/**
*@file
*@brief contiene la clase newton y sus rutinas
*@author MONICA ALEJANDRA CASTELLANOS MENDEZ <monicacastellanos@unicauca.edu.co>
*@author ASHLEE VANESSA CAMPAZ VALENCIA <ashleecampaz@unicauca.edu.co>
*@copyright MIT License
*/
namespace interpolacion{
	
	class newton{
		
	public: 
		newton(vector<double>x, vector<double> y){
			if(x.size()==0 || x.size()!=y.size()){
				valido = false;
				return; 
			}
			this->x = x;
			this->y = y; 
			construir();
		}
		
		bool es_valido(){
			return valido; 
		}
		
		double interpolar(double x_int){
			if(!valido){
				return NAN; 
			}
			return evaluar(x_int); 
		}
		
		void imprimir(ostream & os){
			if(es_valido()==0) {
				return;
			}
			
			os<<"y = "; 
			os<< (b[0]<0? " - ": " + ")<<fabs(b[0]);
			for(size_t i = 1; i < b.size(); i++){
				os<< (b[i]<0? " - ": " + ")<<fabs(b[i]);
				for(size_t j = 0; j<i; j++){
					os<< "(x "
						<<(x[j]>0?" - ":" + ")
						<<fabs(x[j])
						<<")";
				}
			}
			os<< endl;
		} 
		double obtenerUltimoCoeficiente(){
			return b[x.size()-1];
		}	
	private:
				
		void construir(){
			double n = x.size();
			//creamos la matriz con n filas
			vector<vector<double>> F(n); 
			
			for(size_t i=0; i<n;i++){
				F[i].resize(n-i); 
			}
			//rellenar la primera columna de la matriz
			for(size_t i=0; i<n;i++){
				F[i][0] = y[i]; 
			}
			for(size_t j=1; j<n;j++){
				for (size_t i=0; i<n -j;i++){
					F[i][j]= (F[i+1][j-1]-F[i][j-1])/(x[i+j] -x[i]); 
					
				}
			}
			b= F[0];
			
			valido = true;
		}
		double evaluar(double int_x){
			double suma = b[0]; 
			double n = x.size(); 
			for(int i=1; i<n; i++){
				double prod = 1.0f;
				for(int j=0; j<i; j++){
					prod*= (int_x - x[j]);
				}
				suma += b[i]*prod;
			}
			
			return suma; 
		}
		
		
		vector<double> x;
		vector<double> y;
		vector<double>  b; 
		bool valido = false; 
	};
	
	
}
#endif
