#include <iostream>
#include <cmath>

using namespace std;

class Carros {
  public : float t; 
  public: float v; 

  public: float posicion(){
    float x;
    for (int i = 0; i < 100; i=i+0.1) { 
      t=t+i;
      x=v*t;
    }
    
  }

    Carros(){
    cout<<"El valor inicial de valocidad del carro:";
    cin>>v;
    cout<< "Intervalo de tiempo en el que incia el carro:";
    cin>> t;
}
   
};


int main(){

Carros Carro1 = Carros();
Carros Carro2 = Carros(); 
cout<< Carro1.posicion()<< endl; 


   
  return 0; 
  
}
