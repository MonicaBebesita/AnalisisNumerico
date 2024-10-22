#include <iostream>
#include <vector>
#include "caso.h"
#include "regresion.h"
#include "regresion_cuadratica.h"
#include "newthon.h"
#include "lagrange.h"
#include "interpolacion_inversa.h"
#include "spline3.h"

using std::cout;
using std::endl;
using std::vector;
using interpolacion::newton; 
using interpolacion::lagrange; 
using interpolacion::intInversa;
using regresion::lineal_simple;
using regresion::lineal_potencia;
using regresion::lineal_exponencial;
using regresion::modelo_lineal;
using regresion::modelo_potencia;
using regresion::logD;
using interpolacion::spline3;
using regresion::modelo_cuadratico;
using regresion::modelo_exponencial;
using regresion::cuadratica;

/**
* @file 
* @brief contiene la implementacion de las rutinas de los casos
*/
void caso_1_regresion(){
	cout << "Caso de 1 regresion lineal"<< endl;
	
	vector<double> x {0.0f,		200.0f,	350.0f,	550.0f,	800.0f,	1000.0f,	1250.0f,	1500.0f,	1700.0f	};
	vector<double> y {12.90f,	13.0f,	13.10f,	13.30f,	13.90f,	14.20f,		14.30f,		14.60f,		14.80f	};
	
	double xi = 500.0f;
	double y_estimado;
	lineal_simple r(x,y); 
	

	y_estimado = r.estimar(xi);
	modelo_lineal modelo = r.obtenerModelo();
	
	cout<<modelo;
	cout<<"El valor estimado para x="<<xi
		<<" es y="<<y_estimado<<endl;
}

	void caso_1_funcion_potencia(){
		/**Remplazamos el 0 por 10^-10 para que no cause problemas con los logaritmos*/
		cout<<"Caso 1 regresion modelo potencial."<<endl;
		vector<double> x {0.0000000001f,		200.0f,	350.0f,	550.0f,	800.0f,	1000.0f,	1250.0f,	1500.0f,	1700.0f	};
		vector<double> y {12.90f,	13.0f,	13.10f,	13.30f,	13.90f,	14.20f,		14.30f,		14.60f,		14.80f	};
		
		vector<double> xlog10 = logD(x);
		double xi = 500.0f;
		double y_estimado;
		lineal_potencia r(x,y); 
		
		
		y_estimado = r.estimar(xi);
		modelo_potencia modelo = r.obtenerModelo();
		
		cout<<modelo;
		cout<<"El valor estimado para x="<<xi
			<<" es y="<<y_estimado<<endl;
	}
	
		
		
		void caso_1_modelo_exponencial(){
			cout<<"Caso 1 regresion modelo exponencial."<<endl;
			vector<double> x {0.0f,		200.0f,	350.0f,	550.0f,	800.0f,	1000.0f,	1250.0f,	1500.0f,	1700.0f	};
			vector<double> y {12.90f,	13.0f,	13.10f,	13.30f,	13.90f,	14.20f,		14.30f,		14.60f,		14.80f	};
			
			
			double xi = 500.0f;
			double y_estimado;
			lineal_exponencial r(x,y); 
			
			
			y_estimado = r.estimar(xi);
			modelo_exponencial modelo = r.obtenerModelo();
			
			cout<<modelo;
			cout<<"El valor estimado para x="<<xi
				<<" es y="<<y_estimado<<endl;
		}	
		
		
		
	void caso_1_regresion_cuadratica(){
		cout<<"Caso 1 regresion_cuadratica ."<<endl;
		vector<double> x {0.0f,		200.0f,	350.0f,	550.0f,	800.0f,	1000.0f,	1250.0f,	1500.0f,	1700.0f	};
		vector<double> y {12.90f,	13.0f,	13.10f,	13.30f,	13.90f,	14.20f,		14.30f,		14.60f,		14.80f	};
		
		vector<double> xlog10 = logD(x);
		double xi = 500.0f;
		double y_estimado;
		cuadratica c(x,y); 

		y_estimado = c.estimar(xi);
		modelo_cuadratico modelo = c.obtener_modelo();
		cout<<modelo;
		cout<<"\nEl valor estimado para x="<<xi
			<<" es y="<<y_estimado<<endl;
	}
	
	void caso_1_interpolacion_newton(){
		cout<<"Caso 1 interpolacion de Newton"<<endl;
		vector<double>  x = {100.0f,200.0f,300.0f,400.0f,500.0f};
		vector<double> y = {-160.0f,-35.0f,-4.2f,9.0f,16.9f};
		
		newton n(x,y);
		
		double int_x = 140.0f;
		double y_est;
		
		if(!n.es_valido()){
			cout<<"El polinomio no es valido."<<endl;
		}
		else{
			y_est = n.interpolar(int_x);
			n.imprimir(cout);
			cout<<"\nPara x = "<<int_x <<" la interpolacion es  y="<<y_est<<endl; 
		}
	}
	
	void caso_1_interpolacion_lagrange(){
		cout<<"Caso 1 interpolacion de lagrange"<<endl;
		vector<double>  x = {100.0f,200.0f,300.0f,400.0f,500.0f};
		vector<double> y = {-160.0f,-35.0f,-4.2f,9.0f,16.9f};
		
		lagrange l(x,y);
		
		double int_x = 140.0f;
		double y_est;
		
		if(!l.es_valido()){
			cout<<"El polinomio no es valido."<<endl;
		}
		else{
			y_est = l.interpolar(int_x);
			cout<<"\nPara x = "<<int_x <<" la interpolacion es  y="<<y_est<<endl; 
		}
		
	}	
	void caso_2_interpolacion_lagrange(){
		cout<<"Caso 2 interpolacion de lagrange"<<endl;
		vector<double>  x = {2.0f,2.2f,2.4f,2.6f,2.8f};
		vector<double> y = {0.5103757f,0.5207843f,0.5104147f,0.4813306f,0.4359160f};
		
		lagrange l(x,y);
		
		double int_x = 2.5f;
		double y_est;
		
		if(!l.es_valido()){
			cout<<"El polinomio no es valido."<<endl;
		}
		else{
			y_est = l.interpolar(int_x);
			cout<<"\nPara x = "<<int_x <<" la interpolacion es  y="<<y_est<<endl; 
		}
		
	}
	void caso_3_interpolacion_lagrange(){
		cout<<"\nCaso 3 interpolacion de lagrange con grado\n"<<endl;
		vector<double>  x = {2.0f,2.2f,2.4f,2.6f,2.8f};
		vector<double> y = {0.5103757f,0.5207843f,0.5104147f,0.4813306f,0.4359160f};
		
		lagrange l(x,y);
		
		double int_x = 2.5f;
		int grado = 2; 
		double y_est;
		
		if(!l.es_valido()){
			cout<<"El polinomio no es valido."<<endl;
		}
		else{
			y_est = l.interpolar(int_x,grado);
			cout<<"Empleado un polinomio de grado "<<grado
			<<"\npara x = "<<int_x <<" la interpolacion es  y="<<y_est<<endl;
			cout<<"Utilizando el intervalo:[";
			for(auto x:l.getIntervalo()){
				cout<<x<<" ";
			}
			cout<<"]"<<endl;
			cout<<"r2 = "<<l.getError()<<endl;
		}
		
	}
	void caso_1_interpolacion_inversa(){
		cout<<"\nCaso  interpolacion inversa\n"<<endl;
		vector<double>  x = {1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f};
		vector<double> y = {1.0f,0.5f,0.3333f,0.25f,0.2f,0.1667f,0.1429f};
		
		intInversa intv(x,y);
		
		double int_y = 0.3f;
		int grado = 2; 
		double x_est;
		
		if(!intv.es_valido()){
			cout<<"El polinomio no es valido."<<endl;
		}
		else{
			x_est = intv.interpolar(int_y);
			cout<<"Empleado un polinomio de grado "<<grado
				<<"\npara y = "<<int_y <<" la interpolacion es  x="<<x_est<<endl;
			cout<<"Utilizando el intervalo:[";
			for(auto x:intv.getIntervalo()){
				cout<<x<<" ";
			}
			cout<<"]"<<endl;
			cout<< intv.getPolinomio()<<endl;
		}
	}
		
	void caso_4_interpolacion_lagrange_g2(){
		cout<<"\nInterpolacion de lagrange con grado 2\n"<<endl;
		vector<double>  x = {5.0f,5.25f,5.50f,5.75f,6.0f,6.25f,6.50f,6.75f};
		vector<double> y = {143.9264189f,165.3421693f,174.3447736f,161.4457388f,114.5160168f,19.0198411f,-141.2134580f,-382.4545422f};
		
		lagrange l(x,y);
		
		double int_x = 5.17f;
		double int_x2 = 6.32f;
		int grado = 2; 
		double y_est;
		double y_est2;
		
		if(!l.es_valido()){
			cout<<"El polinomio no es valido."<<endl;
		}
		else{
			y_est = l.interpolar(int_x,grado);
			cout<<"Empleado un polinomio de grado "<<grado
				<<"\npara x = "<<int_x <<" la interpolacion es  y="<<y_est<<endl;
			cout<<"Utilizando el intervalo:[";
			for(auto x:l.getIntervalo()){
				cout<<x<<" ";
			}
			cout<<"]"<<endl;
			cout<<"r2 = "<<l.getError()<<endl;
			
			y_est2 = l.interpolar(int_x2,grado);
			cout<<"\nEmpleado un polinomio de grado "<<grado
				<<"\npara x = "<<int_x2 <<" la interpolacion es  y="<<y_est2<<endl;
			cout<<"Utilizando el intervalo:[";
			for(auto x:l.getIntervalo()){
				cout<<x<<" ";
			}
			cout<<"]"<<endl;
			cout<<"r2 = "<<l.getError()<<endl;
		}
		
	}	
	void caso_4_interpolacion_lagrange_g3(){
		cout<<"\nInterpolacion de lagrange con grado 3\n"<<endl;
		vector<double>  x = {5.0f,5.25f,5.50f,5.75f,6.0f,6.25f,6.50f,6.75f};
		vector<double> y = {143.9264189f,165.3421693f,174.3447736f,161.4457388f,114.5160168f,19.0198411f,-141.2134580f,-382.4545422f};
		
		lagrange l(x,y);
		
		double int_x = 5.17f;
		double int_x2 = 6.32f;
		int grado = 3; 
		double y_est;
		double y_est2;
		
		if(!l.es_valido()){
			cout<<"El polinomio no es valido."<<endl;
		}
		else{
			y_est = l.interpolar(int_x,grado);
			cout<<"Empleado un polinomio de grado "<<grado
				<<"\npara x = "<<int_x <<" la interpolacion es  y="<<y_est<<endl;
			cout<<"Utilizando el intervalo:[";
			for(auto x:l.getIntervalo()){
				cout<<x<<" ";
			}
			cout<<"]"<<endl;
			cout<<"r2 = "<<l.getError()<<endl;
			
			y_est2 = l.interpolar(int_x2,grado);
			cout<<"\nEmpleado un polinomio de grado "<<grado
				<<"\npara x = "<<int_x2 <<" la interpolacion es  y="<<y_est2<<endl;
			cout<<"Utilizando el intervalo:[";
			for(auto x:l.getIntervalo()){
				cout<<x<<" ";
			}
			cout<<"]"<<endl;
			cout<<"r2 = "<<l.getError()<<endl;
		}
		
	}	
	void trazador_cubico(){
		cout << "Trazador Cubico" << endl;
		vector<double> x = {6.00, 6.25, 6.50, 6.75};
		vector<double> y = {114.5160168, 19.0198411, -141.2134580, -382.4545422
		};
		spline3 sp(x,y);
		double x_int =  6.3125;
		double y_int;
		
		if(!sp.es_valido()){
			cout << "No es valido" << endl;
			return;
		}else{
			y_int = sp.interpolar(x_int);
			cout<<"\nPara x = "<<x_int<<" la interpolacion es y="<<y_int<<endl;
			sp.imprimir_error(0.775847881,y_int);
		}
	}
