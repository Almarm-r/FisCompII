#include<iostream>
#include<vector>
#include<cstdlib>
#include<fstream>
#include<Eigen/Dense>


void ingreso_datos(std:: vector<double>& x, std:: vector<double>&  y){
    int npoints;
    std:: cout <<"Introduce the amount of data " << std :: endl;
    std:: cin >> npoints ;
    x.resize(npoints);
    y.resize(npoints);

    std:: cout << "Intoroduce the data in  (x,y) format: " <<std:: endl;

    for (int i = 0 ; i < npoints ; ++i ){
        std:: cout << "Point"  <<  i+1 << " :\n" ;
        std :: cin >> x[i] >> y[i];

    }; 
};

class Interpolation {
   public:
    int target;
    double lineal_interpolation(const std::vector<double>& x, const std::vector<double>& y, double target){
    int n = x.size();
    if (target < x[0]|| target > x[n-1]){
    std:: cerr << "The value is out of the range" << std:: endl;
    return 0.0;
    
    }
    //Found the target's interval 
    int i =0;
    while (i< n && target > x[i]){
        i++;
    };
    // Calculate the lineal interpolation 
    if (i==0){
        return y[0];
    }else {
        double slope = (y[i]- y[i-1])/(x[i]-x[i-1]);
        return y[i-1]+ slope *(target - x[i-1]);
    };
};
double lagrange_iterpolation(const std::vector<double>& x, const std::vector<double>& y, double z ){
    double  z, valor = 0; 
    for(int i=0; i<n ;i++){
        l=y[i];
        for(int j=0; j<n; j++){
          if(i!=j){
             l=(l*(z-x[j]))/(x[i]-x[j]);
            }
        }
        std ::cout<<std::endl<<std::endl<<"El valor al polinomio de interpolacion en Z= "<<z <<" es : "<<valor;
       return valor=valor+l;
 };
};
double lagrange_polinomials(const std::vector<double>& x, const std::vector<double>& y, double z ){
    double  z, valorp = 0, pol;   
          if(i!=j){
             pol=(z-x[j])/(x[i]-x[j]);
              return valorp=*pol;
            };
 };



double interpolacion_newton(std::vector<double>& x, std::vector<double>& y, double xi){
    // x: apuntador vector de valores de x, y: puntador vector de valores de y, xi: apuntador valor a interpolar
    // n: numero de puntos, b: matriz de diferencias divididas, yi: valor interpolado
    int n = y.size();
    std::vector<std::vector<double>> b(n, std::vector<double>(n, 0.0));
    double yi;
    double xt = 1;  //variable temporal para almacenar el producto de los factores (xi - xj)
    //inicializacion de la matriz de diferencias divididas
    b[0] = y;
    //calculo de las diferencias divididas
    for( int j = 1; j < n; j++){
        for( int i = 0; i < n - j; i++){
            b[i][j] = (b[i + 1][j - 1] - b[i][j - 1]) / (x[i + j] - x[i]);
        }
    }
    //interpolacion con el polinomio de Newton
    yi = b[0][0];   //primer fila de la matriz de diferencias divididas es primer coeficiente

    for (int j = 0; j < n - 1; j++){
        xt *= (xi - x[j]);  //calculo de los factores (xi - xj)
        yi += b[0][j + 1] * xt; //calculo del polinomio de Newton
    };
        
    return yi;

};
//Constructor 
Interpolation(){
    
    std:: cout << "Introduce the value for do the interpolation: " ;
    std :: cin >> target;
};
};

class Regresiones{
    private:
    Eigen::VectorXd X;  // Datos de características
    Eigen::VectorXd y;  // Datos de etiquetas
    int grado;  // Grado del polinomio para regresión polinomial

public:
   Regresion(const Eigen::VectorXd& features, const Eigen::VectorXd& labels, int polyDegree = 1)
        : X(features), y(labels), grado(polyDegree) {}

    void regresionLineal() {
        // Calcular las sumatorias necesarias
        double sumX = 0.0, sumY = 0.0, sumXY = 0.0, sumX2 = 0.0;
        int n = X.size();

        sumX = X.sum();
        sumY = y.sum();
        sumXY = (X.array() * y.array()).sum();
        sumX2 = (X.array() * X.array()).sum();

        // Calcular los coeficientes de la regresión lineal
        double a = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
        double b = (sumY * sumX2 - sumX * sumXY) / (n * sumX2 - sumX * sumX);

        std::cout << "Coeficientes de la regresion lineal: a = " << a << ", b = " << b << "\n";
    }

    void regresionPolinomial() {
        int n = X.size();

        // Crear la matriz de características para la regresión polinomial
        Eigen::MatrixXd X_poly = Eigen::MatrixXd::Zero(grado + 1, grado + 1);
        Eigen::VectorXd Y = Eigen::VectorXd::Zero(grado + 1);

        for (int i = 0; i < n; i++){
            Eigen::VectorXd powers = Eigen::VectorXd::Zero(grado + 1);  //vector de potencias de x
            for (int j = 0; j <= grado; j++){
                powers(j) = std::pow(X(i), j);
            }
            X_poly += powers * powers.transpose();
            Y += y(i) * powers;
        }

        // Resolver el sistema de ecuaciones usando eliminación de Gauss

        for (int i = 0; i < grado + 1; i++){
            double divisor = X_poly(i, i);
            
            X_poly.row(i) /= divisor;
            Y(i) /= divisor;

            for(int k = i + 1; k < grado + 1; k++){
                double factor = X_poly(k, i);
                X_poly.row(k) -= factor * X_poly.row(i);
                Y(k) -= factor * Y(i);
            }
        }

        // Resolver hacia atrás
        Eigen::VectorXd coeficientes = Eigen::VectorXd::Zero(grado + 1);

        for (int i = grado; i >= 0; --i) {
            coeficientes[i] = Y(i);
            for (int j = i + 1; j <= grado; ++j) {
                coeficientes[i] -= X_poly(i,j) * coeficientes[j];
            }
        }

        std::cout << "Coeficientes de la regresion polinomial:\n";
        for (int i = 0; i <= grado; ++i) {
            std::cout << "a" << i << " = " << coeficientes[i] << "\n";
        }
    }

};




int main (){
cout << "Opcion 1 INTERPOLACION LINEAL \n ";
int opcion = 0;
cin >> opcion;

switch(opcion)
{
    case 1: cout << "----INTERPOLACION LINEAL----";
    std :: vector<double> x, y;
    Interpolation lineal = Interpolation();
    ingreso_datos(x,y);
    lineal.lineal_interpolation(x,y,lineal.target);

    std:: ofstream plot("lineal.gnu");

        plot << "set title 'Lineal Interpolation'\n" ;
        plot << "set xlabel'x'\n" ;
        plot << "set ylabel 'y' \n" ;
        plot << "set xrange [0 : 5] \n" ;

        plot << "plot '-' with lp title 'Points', '-' with points title 'Interpolation'\n";

        for (size_t i = 0 ; i< x.size(); ++i){
            plot << x[i] << "" << y[i] << "\n";
        }
    plot << "e\n" ;
    plot << lineal.target << " " << lineal.lineal_interpolation (x,y , lineal.target) << "\n";
    plot << "e\n";
    plot << "pause -1";
    plot.close();

    std:: system("gnuplot lineal.gnu");
    case 2: cout << "----INTERPOLACION DE LAGRANGE----";


    
    case 3: cout << "Usted ha seleccionado la opción 3";
    break;
    default: cout << "Usted ha ingresado una opción incorrecta";
}

return 0 ; 
};
