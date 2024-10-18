#include<iostream>
#include "caso.h"
#include <vector>
#include "gauss.h"
using std::cout;
using std::endl;
using std::cin;
using std::vector;
void prueba_gauss();
int main (int argc, char *argv[]) {
	
	caso_1_regresion();
	//caso_1_funcion_potencia();
	//prueba_gauss();
	//caso_1_regresion_cuadratica();
	//caso_1_interpolacion_newton();
	//caso_2_interpolacion_lagrange(); 
	//caso_3_interpolacion_lagrange();
	//caso_1_interpolacion_inversa();
	return 0;
}

void prueba_gauss(){
	
	vector<vector<double>> m_gauss = {
	{6.0f, 15.0f,55.0f,152.6f},	
	{15.0f, 55.0f, 225.0f,585.6f},	
	{55.0f,225.0f,979.0f,2488.8f}	
	};
	
	vector<double> a = gauss(m_gauss);
		
	for (size_t i = 0; i <a.size(); i++){
		
		cout<<"a["<< i <<"]="<<a[i]<<endl;
			
	}
}
