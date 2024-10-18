#ifndef INTERPOLACION_INVERSA_H
#define INTERPOLACION_INVERSA_H

#include "mm.h"
#include <string>
using raices::mm;

namespace interpolacion{
	
	class intInversa{
		
	public:
		intInversa (vector<double> x, vector<double> y){
			if(x.size()==0 || x.size()!=y.size()){
				valido = false;
				return; 
			}
			this->x = x;
			this->y = y; 
			this->intervalo_y.resize(3);
			this->intervalo_x.resize(3);
			
			
		}
		
		bool es_valido(){
			return valido; 
		}
		double interpolar(double y_int){
			if(!valido){
				return NAN; 
			}
			
			return evaluar(y_int); 
		}
		void construirPolinomio(double int_y){
			using std::to_string;
			if(es_valido()==0) {
				return;
			}
			polinomio = (b[0]<0? " ~ ":"") + to_string(fabs(b[0]));
			for(size_t i = 1; i < b.size(); i++){
				polinomio += (b[i]<0? " - ": " + ") + to_string(fabs(b[i]));
				for(size_t j = 0; j<i; j++){
					polinomio +="*(x "
						+string((intervalo_x[j]>0? " - " : " + "))
						+to_string(fabs(intervalo_x[j]))
						+")";
				}
			}
			polinomio +=(int_y<0? " + ":" - ") + to_string(fabs(int_y)); 
		}
			vector<double> getIntervalo(){
				return intervalo_y; 
			}
		string getPolinomio(){
			return polinomio;
		}
	private:
	void construir(){
		double n = intervalo_x.size();
		//creamos la matriz con n filas
		vector<vector<double>> F(n); 
		
		for(size_t i=0; i<n;i++){
			F[i].resize(n-i); 
		}
		//rellenar la primera columna de la matriz
		for(size_t i=0; i<n;i++){
			F[i][0] = intervalo_y[i]; 
		}
		for(size_t j=1; j<n;j++){
			for (size_t i=0; i<n -j;i++){
				F[i][j]= (F[i+1][j-1]-F[i][j-1])/(intervalo_x[i+j] -intervalo_x[i]); 
				
			}
		}
		b= F[0];
		
		valido = true;
	}		
	void hallar_intevalo(double int_y){
		int num_datos = 3; 
		intervalo_x.resize(num_datos);
		intervalo_y.resize(num_datos);
		int contador = 0; 
		//este metodo asume que x esta ordenado en orden ascendente
		int indice_lim_inf = hallar_limiteInferior(int_y); 
		int indice_lim_sup = hallar_limiteSuperior(int_y); 
		int numeros_arriba_lim_inf = (x.size() -1)-indice_lim_inf;
		int numeros_debajo_lim_sup =  indice_lim_sup; 
		int i;
		int j; 

		int num_dato1=(num_datos-1)/2;
		int num_dato2 = num_datos - num_dato1; 
		
		//calcula el numero de datos que se deben tomar por debajo y por arriba de x_int
		if(numeros_arriba_lim_inf>=num_dato1-1){
			
			if(numeros_debajo_lim_sup>=num_dato2-1){
				i= indice_lim_inf - num_dato1 +1;
				//la mayor parte de los datos se toman por arriba
				j= indice_lim_sup + num_dato2 -1; 
				
				int lim_inf1= i; 
				int lim_sup2 = j; 
				if(numeros_arriba_lim_inf>=num_dato2-1){
					
					i= indice_lim_inf - num_dato2 +1;
					//la mayor parte de los datos se toman por arriba
					j= indice_lim_sup + num_dato1 -1; 
				}
				
				if(y[j]-y[i]>y[lim_sup2]-y[lim_inf1]){
					i=lim_inf1;
					j = lim_sup2;
				}
			}else{
				i= indice_lim_inf - num_dato2 +1;
				//la mayor parte de los datos se toman por arriba
				j= indice_lim_sup + num_dato1 -1;
			}
		}
		else{
			i=0;
			//la parte superior compensa los datos que no se pueden
			//tomar por debajo de x_int
			j = indice_lim_sup + num_dato2 + num_dato1 -1;
		}
		
		while (i<= indice_lim_inf){
			intervalo_x[contador]=x[i];
			intervalo_y[contador]=y[i]; 
			i++;
			contador++; 
		}
		i=indice_lim_sup;
		while (i<=j && i<x.size()){
			intervalo_x[contador]=x[i];
			intervalo_y[contador]=y[i]; 
			i++;
			contador++; 
		}
		construir();
	}	
	int hallar_limiteInferior(double int_y){
		//devuelve el indice de y en donde se encuentra el limite inferior
		int i=1; 
		double limite_sup = y[0];
		while(y[i]>int_y && y[i]<limite_sup){
			limite_sup = y[i]; 
			i++; 
		}
		return i-1; 
	}
		int hallar_limiteSuperior(double int_y){
			
			//devuelve el indice de y en donde se encuentra el limite superior
			int i=y.size()-2; 
			double limite_inf = y[y.size()-1];
			while(y[i]<int_y && y[i]>limite_inf){
				limite_inf = y[i]; 
				i--; 
			}
			return i+1; 
		}
		double evaluar(double int_y){
			hallar_intevalo( int_y); 
			construirPolinomio(int_y);
			mm muller(polinomio);
			double erp =1e-5;
			int n = 100;
			return muller.encontrar(intervalo_x[0],intervalo_x[1],intervalo_x[2],erp,n);
		}	
		vector<double> x;
		vector<double> y;
		vector<double> intervalo_x;
		vector<double> intervalo_y;
		vector<double>  b; 
		string polinomio; 
		bool valido = true;
	};
}
#endif
