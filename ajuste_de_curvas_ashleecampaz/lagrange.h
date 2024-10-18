#ifndef LAGRANGE_H
#define LAGRANGE_H

//metodo de interpolacion
namespace interpolacion{
	
	class lagrange{
	public:
		lagrange(vector<double>x, vector<double> y){
			if(x.size()==0 || x.size()!=y.size()){
				valido = false;
				return; 
			}
			valido = true;
			this->x = x;
			this->y = y; 
			
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
		
		double interpolar(double x_int, int grado){
			if(!valido){
				return NAN; 
			}
			return evaluar(x_int,grado); 
		}	
		vector<double> getIntervalo(){
			return intervalo_x; 
		}
		double getError(){
			return r2; 
		}	
	private: 	
		void hallar_intevalo(double int_x, int grado){
			int num_datos = grado +1; 
			intervalo_x.resize(num_datos);
			intervalo_y.resize(num_datos);
			int contador = 0; 
			//este metodo asume que x esta ordenado en orden ascendente
			int indice_lim_inf = hallar_limiteInferior(int_x); 
			int indice_lim_sup = hallar_limiteSuperior(int_x); 
			int numeros_debajo_lim_inf = indice_lim_inf;
			int numeros_arriba_lim_sup = (x.size() -1)- indice_lim_sup; 
			int i;
			int j; 
			if(grado%2==0){
				//numero de datos es impar
				int num_dato1=(num_datos-1)/2;
				int num_dato2 = num_datos - num_dato1; 
				
				//calcula el numero de datos que se deben tomar por debajo y por arriba de x_int
				if(numeros_debajo_lim_inf>=num_dato1-1){
					
					if(numeros_arriba_lim_sup>=num_dato2-1){
						i= indice_lim_inf - num_dato1 +1;
						//la mayor parte de los datos se toman por arriba
						j= indice_lim_sup + num_dato2 -1; 
						
						int lim_inf1= i; 
						int lim_sup2 = j; 
						if(numeros_debajo_lim_inf>=num_dato2-1){
							
							i= indice_lim_inf - num_dato2 +1;
							//la mayor parte de los datos se toman por arriba
							j= indice_lim_sup + num_dato1 -1; 
						}

						if(x[j]-x[i]>x[lim_sup2]-x[lim_inf1]){
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
			}
			else{
				//numero de datos es par
				//calcula el numero de datos que se deben tomar por debajo de x_int
				if(numeros_debajo_lim_inf>=(num_datos/2 -1)){
					if(numeros_arriba_lim_sup>=(num_datos/2 -1)){
						i= indice_lim_inf - (num_datos/2) +1;
					}
					else{
						i = indice_lim_inf - (num_datos/2) +1  - ((num_datos/2) -numeros_arriba_lim_sup-1); 
					}
				}else{
					i =0;
				}
				while (i<= indice_lim_inf){
					intervalo_x[contador]=x[i];
					intervalo_y[contador]=y[i]; 
					i++;
					contador++; 
				}
				//calcula el numero de datos que se deben tomar por encima de x_int
				if(numeros_debajo_lim_inf>=(num_datos/2 -1)){
					j = indice_lim_sup + (num_datos/2) -1; 
				}
				else{
					j = indice_lim_sup + num_datos  - numeros_debajo_lim_inf -1; 
				}
				i=indice_lim_sup;
				while (i<=j && i<x.size()){
					intervalo_x[contador]=x[i];
					intervalo_y[contador]=y[i]; 
					i++;
					contador++; 
				}
			}
			//obtiene los intervalos para calcular el error 
			vector<double> intervalo_x_error = intervalo_x;
			vector<double> intervalo_y_error = intervalo_y; 
			intervalo_x_error.resize(intervalo_x.size()+1);
			intervalo_y_error.resize(intervalo_x.size()+1); 
			if(numeros_arriba_lim_sup>0){
				intervalo_x_error[intervalo_x.size()-1] = x[indice_lim_sup+1];
				intervalo_y_error[intervalo_x.size()-1] = y[indice_lim_sup+1];
			}
			else{
				for(int i=1; i<=intervalo_x_error.size();i++){
					intervalo_x_error[i]= intervalo_x_error[i-1];
					intervalo_y_error[i]= intervalo_y_error[i-1];
				}
				intervalo_x_error[0]=x[indice_lim_inf-1];
				intervalo_y_error[0]=y[indice_lim_inf-1];
			}
			newton n(intervalo_x_error,intervalo_y_error); 
			calcular_error(int_x,n);
		}
		
		void calcular_error(int int_x,newton n){
			
			r2 = n.obtenerUltimoCoeficiente();
			for(int i=0; i<intervalo_x.size();i++){
				r2*= (int_x - intervalo_x[i]);
			}
			
		}	
			
		int hallar_limiteInferior(double int_x){
			//devuelve el indice de x en donde se encuentra el limite inferior
			int i=1; 
			double limite_inf = x[0];
			while(x[i]<int_x && x[i]>limite_inf){
				limite_inf = x[i]; 
				i++; 
			}
			return i-1; 
		}
		int hallar_limiteSuperior(double int_x){
			//devuelve el indice de x en donde se encuentra el limite superior
			int i=x.size()-2; 
			double limite_sup = x[x.size()-1];
			while(x[i]>int_x && x[i]<limite_sup){
				limite_sup = x[i]; 
				i--; 
			}
			return i+1; 
		}
			
		double evaluar(double int_x){
			double suma = 0; 
			double n = x.size(); 
			for(int i=0; i<n; i++){
				double denominador= 1.0f;
				double determinador= 1.0f;
				for(int j=0; j<n; j++){
					if(i!=j){
						denominador*= (int_x - x[j]);
						determinador*=(x[i]-x[j]); 
					}
					
				}
				suma += y[i]*(denominador/determinador);
			}
			
			return suma; 
		}
		
		double evaluar(double int_x, int grado){
			double suma = 0; 
			if(grado>x.size() -1){
				return NAN;
			}
			hallar_intevalo( int_x,grado); 
			
			double n = intervalo_x.size(); 
			
			for(int i=0; i<n; i++){
				double denominador= 1.0f;
				double determinador= 1.0f;
				for(int j=0; j<n; j++){
					if(i!=j){
						denominador*= (int_x - intervalo_x[j]);
						determinador*=(intervalo_x[i]-intervalo_x[j]); 
					}
					
				}
				suma += intervalo_y[i]*(denominador/determinador);
			}
			
			return suma; 
		}
				
		vector<double> x;
		vector<double> y;
		vector<double> intervalo_x;
		vector<double> intervalo_y;
		double r2; 
		bool valido = false;
	};
	
	
}


#endif
