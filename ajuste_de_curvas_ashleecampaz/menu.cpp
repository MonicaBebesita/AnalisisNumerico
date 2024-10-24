#include<iostream>
#include "caso.h"
#include "menu.h"
#include <conio.h>
using std::cout;
using std::endl;
using std::cin;

/**
*@file
*@brief implementacion de los menus
*@author MONICA ALEJANDRA CASTELLANOS MENDEZ <monicacastellanos@unicauca.edu.co>
*@author ASHLEE VANESSA CAMPAZ VALENCIA <ashleecampaz@unicauca.edu.co>
*@copyright MIT License
*/

void clear(){
	cout<<"Presione una tecla para ir al menu anterior."<<endl;
	getche(); 
	system("cls");
}	

void menuPrincipal(){
	
	int opcion = -1; 
	do{	
		cout<<"--Menu--\n1.Regresion\n2.Interpolacion\n3.Trazador\n4.salir"<<endl; 
		cout<<"ingrese un opcion: "<<endl; 
		scanf("%i",&opcion); 
		switch(opcion){
		case 1: 
			system("cls");
			menuRegresion();
			//clear();
			break;
		case 2: 
			system("cls");
			menuInterpolacion();
			//clear();
			break;
		case 3:
			system("cls");
			trazador_cubico();
			clear();
			break;
		case 4:
			system("cls");
			break;
		default:
			cout<<"Opcion incorrecta, intente nuevamente."<<endl;
		}
		
	}while(opcion!=4);	
}

void menuRegresion(){
	int opcion = -1; 
	do{
		cout<<"--Menu regresion--\n1.Modelo lineal\n2.Modelo de potencia\n3.Modelo exponencial\n4.Modelo cuadratico\n5.volver atras "<<endl; 
		cout<<"ingrese un opcion: "<<endl; 
		scanf("%i",&opcion); 
		switch(opcion){
		case 1: 
			system("cls");
			caso_1_regresion();
			clear();
			break;
		case 2: 
			system("cls");
			caso_1_funcion_potencia();
			clear();
			break; 
		case 3:
			system("cls");
			caso_1_modelo_exponencial();
			clear();
			break;
		case 4:
			system("cls");
			caso_1_regresion_cuadratica();
			clear();
			break;
		case 5:
			system("cls");
			break;
		default:
			cout<<"Opcion incorrecta, intente nuevamente."<<endl;
		}
		
	}while(opcion!=5);
}

void menuInterpolacion(){
	int opcion = -1; 
	do{
		cout<<"--Menu interpolacion--\n1.Polinomio de grado 2\n2.Polinomio de grado 3\n3.Volver atrás"<<endl;
		cout<<"ingrese un opcion"<<endl; 
		scanf("%i",&opcion); 
		switch(opcion){
		case 1:
			system("cls");
			caso_4_interpolacion_lagrange_g2();
			clear();
			break;
		case 2:	
			system("cls");
			caso_4_interpolacion_lagrange_g3();
			clear();
			break;
		case 3: 
			system("cls");
			break;
		default:
			cout<<"Opcion incorrecta, intente nuevamente."<<endl;
		}
		
	}while(opcion!=3);
}


	
