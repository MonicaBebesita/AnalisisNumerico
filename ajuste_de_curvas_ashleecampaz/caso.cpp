#include <iostream>
#include <vector>
#include "caso.h"
#include "regresion.h"
#include "regresion_cuadratica.h"
#include "newthon.h"
#include "lagrange.h"
#include "interpolacion_inversa.h"
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

using regresion::modelo_cuadratico;
using regresion::modelo_exponencial;
using regresion::cuadratica;

void caso_1_regresion(){
	cout << "Caso de 1 regresion lineal"<< endl;
	
	vector<double> x {10.0f,20.0f,30.0f,40.0f,50.0f,60.0f,70.0f,80.0f};
	vector<double> y {1.06f,1.33f,1.52f,1.68f,1.81f,1.91f,2.01f,2.11f};
	
	double xi = 35.0f;
	double y_estimado;
	lineal_simple r(x,y); 
	

	y_estimado = r.estimar(xi);
	modelo_lineal modelo = r.obtenerModelo();
	
	cout<<modelo;
	cout<<"El valor estimado para x"<<xi
		<<"es y="<<y_estimado<<endl;
}

	void caso_1_funcion_potencia(){
		cout<<"Caso 1 regresion modelo potencial."<<endl;
		vector<double> x {10.0f,20.0f,30.0f,40.0f,50.0f,60.0f,70.0f,80.0f};
		vector<double> y {1.06f,1.33f,1.52f,1.68f,1.81f,1.91f,2.01f,2.11f};
		vector<double> xlog10 = logD(x);
		double xi = 35.0f;
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
			vector<double> x {10.0f,20.0f,30.0f,40.0f,50.0f,60.0f,70.0f,80.0f};
			vector<double> y {1.06f,1.33f,1.52f,1.68f,1.81f,1.91f,2.01f,2.11f};
			
			double xi = 35.0f;
			double y_estimado;
			lineal_exponencial r(x,y); 
			
			
			y_estimado = r.estimar(xi);
			modelo_exponencial modelo = r.obtenerModelo();
			
			cout<<modelo;
			cout<<"El valor estimado para x="<<xi
				<<" es y="<<y_estimado<<endl;
		}	
		
		
		
	void caso_1_regresion_cuadratica(){
		cout<<"Caso 1 linealizacion funcion ."<<endl;
		vector<double> x {10.0f,20.0f,30.0f,40.0f,50.0f,60.0f,70.0f,80.0f};
		vector<double> y {1.06f,1.33f,1.52f,1.68f,1.81f,1.91f,2.01f,2.11f};
		vector<double> xlog10 = logD(x);
		double xi = 35.0f;
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
