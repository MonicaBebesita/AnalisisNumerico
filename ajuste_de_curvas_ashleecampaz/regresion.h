#ifndef REGRESION_H
#define REGRESION_H
#include <iostream>
#include <cmath>
#include <vector>
using std::ostream;
using std::endl;
using std::cout;
using std::cin;
using std::vector;

namespace regresion{
	
	struct modelo_lineal{
		vector <double> x;/*!<Variable independiente*/
		vector <double> y;/*!<Variable dependiente*/
		double b0;/*!<Termino independiente de la ecuacion de la recta de regresion*/
		double b1;/*!<Termino independiente de la ecuacion de la recta de regresion*/
		double st;
		double sr;
		double sy;
		double r2;
		double syx;
		double y_prom;
		bool valido = false;
		modelo_lineal(vector <double> x,vector <double> y):x(x), y(y){
			contruir();
		};
		
		void contruir(){
			double s_x = 0.0f;
			double s_y = 0.0f;
			double s_xy = 0.0f;
			double s_x2 = 0.0f;
			double x_prom;
			double y_prom;
			
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
			}
			
			x_prom = s_x/(double)n;
			y_prom = s_y/(double)n;
			sr = 0.0f;
			st= 0.0f;
			r2= 0.0f;
			sy =0.0f;
			syx = 0.0f;
			b1= (s_xy -(y_prom*s_x))/(s_x2 -x_prom*s_x);
			b0= y_prom -b1*x_prom;
			for(size_t i =0; i<n; i++){
				st +=pow(y[i] - y_prom,2);
				sr += pow(y[i] -b0 -(b1*x[i]),2);
				
			}
			sy = sqrt(st/(double)(n-1));
			syx = sqrt(sr/(double)(n-2));
			r2= (st-sr)/st;
			valido = true; 
		}
			
		double estimar(double x){
			return b1*x + b0;
		}
		friend ostream& operator <<(ostream & os, const modelo_lineal & m){
			if (!m.valido){
				os << "El modelo no es malido"<<endl;
			}
			os<<"Recta de regresion: "
			<< "y="
			<<m.b1<<" *x"
			<<(m.b1<0.0f?"- ":" +")
			<<fabs(m.b0)
			<<endl;
			
			os << "r2 =" <<m.r2 <<endl;
			os << "Desv. estandar:"
			<<m.sy
			<<endl;
			
			os<< "Error estandar de aproximacion:"
			<<m.syx
			<<endl;
			if(m.syx<m.sy){
			os<<"syx<sy"<<endl;
			
				
			}
			return os; 
		}
	};
	
	class lineal_simple{
		public:
			lineal_simple(vector <double> x,vector <double> y):modelo(x,y){
				
			}
			
			double estimar(double x){
				//To do
				if(!modelo.valido){
					return NAN;
				}
				return modelo.estimar(x);
			}
			modelo_lineal obtenerModelo(){
				return modelo;
			}
		private: 
			modelo_lineal modelo;
		
	};
	
	
	
	vector<double>  ln(vector<double> v){
		for (auto & x:v){
			x= log(x);
		}
		return v;
	}
		
	vector<double> logD(vector<double> v){
		for (auto & x:v){
			x= log10(x);
		}
		return v;
	}
	/**
	@brief Funcion potencia y=c*x^a
	*/
	struct modelo_potencia{
		double c;
		double a;
		double valido;
		struct modelo_lineal lineal;
		modelo_potencia(vector<double> x, vector <double>y):lineal(logD(x),logD(y)){
			valido = lineal.valido;
			c = pow(10.0f,lineal.b0);
			a = lineal.b1;
			
		}
		
		friend ostream& operator <<(ostream & os, const modelo_potencia &m){
			if (!m.lineal.valido){
				os << "El modelo no es valido"<<endl;
			}
			os<<"Funcion potencia: "
				<< "y="<<m.c
				<<"*x^"<<m.a
				<<endl;
			
			os << "r2 =" <<m.lineal.r2 <<endl;
			os << "Desv. estandar:"
				<<m.lineal.sy
				<<endl;
			
			os<< "Error estandar de aproximacion:"
				<<m.lineal.syx
				<<endl;
			return os;
		}
		
		double estimar (double x_est){
			return c*pow(x_est,a);
		}
	};
	
	class lineal_potencia{
		public:
			lineal_potencia(vector <double> x,vector <double> y):modelo(x,y){
				
			}
			double estimar(double x){
				if(!modelo.valido){
					return NAN;
				}
				return modelo.estimar(x);
			}
				modelo_potencia obtenerModelo(){
					return modelo;
				}
		private: 
			modelo_potencia modelo;
	};
}
#endif
