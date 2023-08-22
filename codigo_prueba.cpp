#include <iostream>
#include <cmath>

using namespace std;

class Persona{

  public: string Nombre;
  public: int anoNac;
  public: float altura;
  public: float peso;

  public: void edad(){
      int edad;
      edad= 2023-anoNac;
      cout<<"'Tu edad es: "<<edad<<endl;

    }

  public: void indice(){
      float indice;
      indice=peso/(altura*altura);
      cout<<"Tu indice corporal es: "<<indice<<endl;
    }

  Persona(){
    cout<<"Nombre:";
    cin>>Nombre;
    cout<<"Año de nacimiento: ";
    cin>>anoNac;
    cout<<"Peso(kg) ";
    cin>>peso;
    cout<<"Estatura(m) ";
    cin>>altura;
}

};

int main(){
  Persona Jose=Persona();
  Jose.indice();
  Jose.edad();
  return 0;
}


