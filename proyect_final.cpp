#include <iostream>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <Eigen/Dense>


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
  }
};

class Interpolation {
public:
      double target;
    std::string id ;
    std::string funtion_name;
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
    }
    // Calculate the lineal interpolation 
    if (i==0){
        return y[0];
    }else {
        double slope = (y[i]- y[i-1])/(x[i]-x[i-1]);
        return y[i-1]+ slope *(target - x[i-1]);
    }
  };
 double lagrange_Interpolation(const std::vector<double>& x, const std::vector<double>& y, double target) {
        int n = x.size();
        double result = 0.0;

        for (int i = 0; i < n; i++) {
            double term = y[i];
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    term *= (target - x[j]) / (x[i] - x[j]);
                }
            }
            result += term;
        }

        return result;
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
    }   
    return yi;
  };
  void graficar(const std::vector<double>& x, const std::vector<double>& y, double target,  std::string id ) {
    std::ofstream grafica("grafica.gnu");
    std:: string hlp= id;

    grafica << "set title '"<<id<<"'" <<"\n";
    grafica << "set xlabel 'x'\n";
    grafica << "set ylabel 'y'\n";
    grafica << "set xrange [0:5] \n";
    grafica << "set yrange [0:5] \n";
    grafica << "plot '-' with lp title 'Puntos', '-' with points title 'Interpolación' \n";
    for (size_t i = 0; i < x.size(); ++i) {
        grafica << x[i] << " " << y[i] << "\n";
    }

    grafica <<"e\n";
    grafica << target << " " << funtion_name << "\n"; 
    grafica <<"e\n";
    grafica <<"pause -1";
    grafica.close();

    std::system("gnuplot grafica.gnu");
};
  //Constructor 
  Interpolation(){  
    std::cout << "Introduce the value for do the interpolation: " << std::endl;
    std::cin >> target;
    std::cout << "Introduce the name label for plot  " << std::endl;
    std::cin >> id;
    funtion_name = 'lineal_interpolation(x , y, lineal.target)';

  };
};

class Regresiones{
private:
    Eigen::VectorXd features;  // Datos de características
    Eigen::VectorXd labels;  // Datos de etiquetas
    int polyDegree;  // Grado del polinomio para regresión polinomial

public:
    Regresiones(Eigen::VectorXd& X, Eigen::VectorXd& y, int grado){
        features = X;
        labels = y;
        polyDegree = grado;
    };

    Eigen::VectorXd regresion_lineal(Eigen::VectorXd& X, Eigen::VectorXd& y) {
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
        Eigen::VectorXd coef(2);
        coef << a, b;
        return coef;
    }

    Eigen::VectorXd regresion_polinomial(Eigen::VectorXd& X, Eigen::VectorXd& y, int grado) {
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
        return coeficientes;
    }

    void plot_data_and_regression(const Eigen::VectorXd& x, const Eigen::VectorXd& y, const Eigen::VectorXd& coef, const std::string& regressionType) {
        std::ofstream plot("plot.gnu");
            plot << "set title 'Regresión " << regressionType << "'\n";
            plot << "set xlabel 'x'\n";
            plot << "set ylabel 'y'\n";
            plot << "set terminal wxt size 800,600\n";

            // Graficar los datos y la regresión en un solo comando "plot"
            plot << "plot '-' with points title 'Datos' lc rgb 'red', " << coef(0);
            for (int i = 1; i < coef.size(); ++i) {
                plot << " + " << coef(i) << " * x**" << i;
            }
            plot << " title 'Regresión " << regressionType << "' lc rgb 'blue'\n";

            for (int i = 0; i < x.size(); ++i) {
                plot << x[i] << " " << y[i] << "\n";
            }
            plot << "e\n";
    
            plot << "pause -1";
            plot.close();

        std::system("gnuplot plot.gnu");
    }

};




int main (){
  std::cout << "Opcion 1: INTERPOLACION LINEAL" << std::endl;
  std::cout << "Opcion 2: INTERPOLACION DE LAGRANGE" << std::endl;
  std::cout << "Opcion 3: INTERPOLACION DE NEWTON" << std::endl;
  std::cout << "Opcion 4: REGRESION LINEAL" << std::endl;
  std::cout << "Opcion 5: REGRESION POLINOMIAL" << std::endl;
  int opcion = 0;
  std::cin >> opcion;

switch(opcion)
{
  case 1: { std::cout << "----INTERPOLACION LINEAL----\n";
    std :: vector<double> x, y;
    Interpolation lineal  = Interpolation();
    ingreso_datos(x,y);
    lineal.funtion_name = 'ineal_interpolation(x,y,lineal.target)l';
    lineal.lineal_interpolation(x,y,lineal.target);
    lineal.graficar(x,y,lineal.target, lineal.id);
    break;
    };


    case 2: {std::cout << "----INTERPOLACION DE LAGRANGE----\n";
   
    std :: vector<double> x, y;
    Interpolation lagrange  = Interpolation();
    ingreso_datos(x,y);
    lagrange.funtion_name = 'lagrange_Interpolation(x,y,lagrange.target)';
    lagrange.lagrange_Interpolation(x,y,lagrange.target);
    lagrange.graficar(x,y,lagrange.target, lagrange.id);

    break;
    };


    
    case 3: {std:: cout << "----INTERPOLACION NEWTON----\n";
    std :: vector<double> x, y;
    ingreso_datos(x,y);
    Interpolation newton  = Interpolation();
    newton.funtion_name = 'interpolacion_newton(x , y, newton.target)';
    newton.interpolacion_newton(x , y, newton.target);
    newton.graficar(x,y,newton.target, newton.id);
    break;
    
    };

      case 4:{ 
    
            std::cout << "----REGRESION LINEAL----\n";

            std::cout << "Ingrese el número de datos:" << std::endl;
            int N;
            std::cin >> N;
            Eigen::VectorXd features(N);
            Eigen::VectorXd labels(N);

            std::cout << "Ingrese los datos en el formato (x,y):" << std::endl;
            for (int i = 0; i < N; ++i) {
                std::cout << "Dato " << i + 1 << ": ";
                std::cin >> features[i] >> labels[i];
            }

            int polyDegree = 1;

            Regresiones Regresion(features, labels, polyDegree);

            Eigen::VectorXd coeficientes = Regresion.regresion_lineal(features, labels);

            std::cout << "El coeficiente x¹ es: " << coeficientes(0) << "El coeficiente x⁰ es: " << coeficientes(1) << std::endl;

            Regresion.plot_data_and_regression(features, labels, coeficientes, "Lineal");

            break;
        }

        case 5:{
            std::cout << "----REGRESION POLINOMIAL----\n";
            std::cout << "Ingrese el número de datos:" << std::endl;
            int N;
            std::cin >> N;

            Eigen::VectorXd features(N);
            Eigen::VectorXd labels(N);

            std::cout << "Ingrese los datos en el formato (x,y):" << std::endl;
            for (int i = 0; i < N; ++i) {
                std::cout << "Dato " << i + 1 << ": ";
                std::cin >> features[i] >> labels[i];
            }

            int polyDegree;
            std::cout << "Ingrese el grado del polinomio: " << std::endl;
            std::cin >> polyDegree;

            Regresiones Regresion(features, labels, polyDegree);

            Eigen::VectorXd coeficientes = Regresion.regresion_polinomial(features, labels, polyDegree);

            for (int i = 0; i < polyDegree; i++){
                std::cout << "Para el termino x^ "<< i << ". El coeficiente es: " << coeficientes(i) << std::endl;
            }

            Regresion.plot_data_and_regression(features, labels, coeficientes, "Polinómica");

            break;

        } 

  default: std::cout << "Usted ha ingresado una opción incorrecta";
}

return 0 ; 
};