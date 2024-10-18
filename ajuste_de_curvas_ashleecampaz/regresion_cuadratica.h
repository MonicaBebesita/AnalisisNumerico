#ifndef REGRESION_CUADRATICA_H
#define REGRESION_CUADRATICA_H

#include <vector>
#include "gauss.h"
#include <cmath>
using std::vector;

namespace regresion{
	struct modelo_cuadratico{
		vector <double> x;/*!<Variable independiente*/
		vector <double> y;/*!<Variable dependiente*/
		vector<double> a ; /*!<arreglo de coeficientes*/
		double st;
		double sr;
		double sy;
		double r2;
		double syx;
		double y_prom;
		bool valido = false;
		
		modelo_cuadratico(vector <double> x,vector <double> y):x(x), y(y){
			construir();
		};
		
		void construir(){
			double s_x = 0.0f;
			double s_y = 0.0f;
			double s_xy = 0.0f;
			double s_x2 = 0.0f;
			double s_x3 = 0.0f;
			double s_x4 = 0.0f;
			double s_x2y = 0.0f;
			
			size_t n = x.size();
			if(n!= y.size()){
				return;
			}
			
			if(n<3){
				return; 
			}
			for(size_t i =0; i<n; i++){
				s_x += x[i];
				s_y +=y[i];
				s_xy += x[i]*y[i];
				s_x2 += x[i]*x[i];
				s_x3 += x[i] * x[i] * x[i]; // x[i]^3
				s_x4 += x[i] * x[i] * x[i] * x[i]; // x[i]^4
				s_x2y += x[i] * x[i] * y[i]; // x[i]^2 * y[i]
			}
			
			vector<vector<double>> m ={
				{(double)n,s_x,s_x2,s_y},
				{s_x,s_x2,s_x3,s_xy},
				{s_x2,s_x3,s_x4,s_x2y}
			};
			 a = gauss(m);
			
			valido = (a[0]!=NAN && a[1]!=NAN && a[2]!=NAN);
			sr = 0.0f;
			st= 0.0f;
			y_prom = s_y/(double)n;
			for(size_t i =0; i<n; i++){
				st +=pow(y[i] -y_prom,2);
				sr += pow(y[i] -estimar(x[i]),2);
				
			}
			sy = sqrt(st/(double)(n-1));
			syx = sqrt(sr/(double)(n-3));
			r2= (st-sr)/st;
			valido = true; 
		}
		
		double estimar(double x_est){
			if(!valido){return NAN;}
			return a[0] + a[1]*x_est + a[2]*x_est*x_est;
		}
			
		friend ostream& operator <<(ostream & os, const modelo_cuadratico & m){
			if (!m.valido){
				os << "El modelo no es valido"<<endl;
			}
			os<<"Ecuacion cuadratica: "
				<< "y="
				<<m.a[2]<<" *x^2"
				<<(m.a[1]<0.0f?"- ":" +")
				<<fabs(m.a[1])<<"* x"
				<<(m.a[0]<0.0f?"- ":" +")
				<<fabs(m.a[0])
				<<endl;
			
			os << "\n-Datos estadisticos-\n"<<endl;
			os	<<"r2 =" <<m.r2 <<endl;
			os << "Desv. estandar:"
				<<m.sy
				<<endl;
			
			os<< "Error estandar de aproximacion:"
				<<m.syx
				<<endl;
			return os; 
		}	
	};
	
	class cuadratica{
		
	public:
		cuadratica(vector <double> x,vector <double> y):modelo(x,y){}
			
		double estimar(double x_est){
			if (!modelo.valido){
				return NAN;
			}
			return modelo.estimar(x_est);
		}
			
		modelo_cuadratico obtener_modelo(){
			return modelo;
		}
	private:
			modelo_cuadratico modelo;
	};
	
}

#endif
