#include<iostream>
#include<fstream>
#include <cmath>
#include<Eigen/Dense>

class Transformacion_Lorentz {
    public : 
        const double c = 299792458;
        double v, gamma , beta, ct,x,y,z; 
        
    void Lorentz_rectangular(double v , double ct , double x , double y , double z ){
        beta = v/c;
        gamma = 1/ (sqrt(1 - pow(beta , 2)));
         Eigen::VectorXd p(4);
         p(0)= ct;
         p(1)=x;
         p(2)=y;
         p(3)=z;
        //mLlenar la matriz de las transfomraciones con las respectivas expresiones de transformación 
        Eigen::MatrixXd t(4,4);
         t(0,2)=t(0,3)=t(1,2)=t(1,3)=t(2,0)=t(2,1)=t(2,3)=t(3,0)=t(3,1)=t(3,1)=0;
         t(0,0)=t(1,1)= gamma;
         t(1,0)=t(0,1) = - (gamma * beta);
         t(2,2)=t(3,3)=1;

         if (v >= c){
           std :: cout <<"Wow, rompes las leyes de la física hermano ." << std :: endl ; 
         };
         if (gamma==1){
          std :: cout <<"Creo que tu problema no es tan especial para la relatividad especial xD" << std :: endl ; 

         };

         Eigen::VectorXd pt(4);
         pt = t*p ;

         std :: cout <<"La transformacion   \n"<< t<< std :: endl ;
         std :: cout <<"--------------------" << std :: endl ; 
         std :: cout <<"El vector que se transformo   \n"<<"ct= "<< p(0)<<"\nx= "<<p(1)<<"\ny = "<<p(2)<< "\n z= "<< p(3)<<std :: endl ;
         std :: cout <<"--------------------" << std :: endl ; 
         std :: cout <<"El vector transformado  \n"<<"ct= "<< pt(0)<<"\nx= "<<pt(1)<<"\ny = "<<pt(2)<< "\n z= "<< pt(3)<<std :: endl ;
         std :: cout <<"--------------------" << std :: endl ; 
         std :: cout <<"Con un factor gamma  de \n"<< gamma<< std :: endl ;
         std :: cout <<"--------------------" << std :: endl ;  
         std :: cout <<"Con un factor beta  de \n"<< beta << std :: endl ;
         std :: cout <<"--------------------" << std :: endl ; 
    }; 

Transformacion_Lorentz (){
  std:: cout << "Ingrese el valor de la coordenada ct ";
  std :: cin >> ct;
  std :: cout << "Ingrese el valor de la coordenada x ";
  std :: cin >> x;
  std :: cout << "Ingrese el valor de la coordenada y ";
  std :: cin >> y;
  std :: cout << "Ingrese el valor de la coordenada z ";
  std :: cin >> z;
  std :: cout << "Ingrese la velocidad con la que esta viajando su partícula ";
  std :: cin >> v; 
};
 
};

 
int main()
{
 Transformacion_Lorentz a = Transformacion_Lorentz();
 a.Lorentz_rectangular(a.v,a.ct ,a.x,a.y,a.z);


 return 0 ; 
}