#ifndef MM_H
#define MM_H

#include "expresion.h"
namespace raices{
	
	bool es_cero(double x){
		return (fabs(x)<=DBL_EPSILON);
	}
	double error_relativo(double nuevo, double ant){
		return fabs((nuevo-ant)/nuevo);
	}
	/**
	*@brief Metodo de muller
	*@param f_str Funcion como cadena de caracteres
	*/
	class mm{
	public:
		mm(string f_str):f(f_str){}
		/**
		*@brief Encuentra una solucion a la funcion usando el metodo de muller
		*@param x0 valor inicial1
		*@param x1 valor inicial2
		*@param x2 valor inicial3
		*@param tol tolerancia
		*@param n numero de iteraciones
		*/
		double encontrar(double x0,double x1,double x2,double tol, int n){
			using raices::es_cero;
			using raices::error_relativo;
			
			//paso 0. validar si la f en los valores inciales es 0
			if(raices::es_cero(f(x0))){
				return x0;
			}
			if(raices::es_cero(f(x1))){
				return x1;
			}
			if(raices::es_cero(f(x2))){
				return x2;
			}
			//paso 1. se inicia el contador de iteraciones
			int i=2;
			while(i<=n){
				
				
				//paso3. se calcula la siguiente aproximacion
				double h1 = x1-x0;
				double h2 = x2-x1;
				double delta1 =(f(x1)-f(x0))/h1;
				double delta2 = (f(x2)-f(x1))/h2;
				//paso 3.1 se hallan los coeficientes
				double a = (delta2 - delta1)/(h2+h1);
				double b = delta2+(h2*a);
				double dis = sqrt(pow(b,2) -(4*f(x2)*a));
				
				//paso 3.2 
				double e;
				if(fabs(b-dis)<fabs(b+dis)){
					 e = b + dis;
				}
				else{
					 e = b - dis;
				}
				//paso 3.3 hallar la aproximacion
				double h=(-2*f(x2))/e; 
				double p = x2 + h;
				
				//paso 4. se calcula el error relativo en p y x2
				double er = raices::error_relativo(p,x2);
				double erp = er*100.0f;
				//paso 5. se verifica si f(p) es 0 o el erp es inferior a la tolerancia
				if(es_cero(f(p)) || erp<tol){
					return p;
				}
				
				//paso 6.paso aumenta el contador de iteraciones
				i++;
				//paso 7. se asigna el nuevo p0
				x0 = x1;
				x1= x2;
				x2=p;
				
			}
			//No se encontro la raiz
			return NAN;
		}
	private:
		expression f;
	};
	
		
}	
#endif
