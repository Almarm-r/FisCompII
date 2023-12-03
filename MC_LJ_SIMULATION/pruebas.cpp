#include<iostream>
#include<math.h>
#include<cmath>
#include<cstdlib>
#include<ctime>
#include<fstream>
#include <stdio.h>
#include <iomanip>
#include<sstream>
#include<string>

#define _USE_MATH_DEFINES //Para traer las constantes 
double epsilon , rho , sigma , T, rc , kb;
int binNo ;
void read_coord(char *file_name, int n, double *x, double *y, double *z);
void initial_conditions(){
    std:: string line;
    char str[67];
    char *file_name = new char;


    std:: cout << "Epsilon parameter\n"<< std:: endl;
    std:: cin >> epsilon ;

    std:: cout << "Rho parameter\n"<< std:: endl;
    std:: cin >> rho ;

    std:: cout << "Sigma parameter\n"<< std:: endl;
    std:: cin >> sigma ;

    std:: cout << "Delta parameter\n"<< std:: endl;
    std:: cin >> rho ;

    std:: cout << "T parameter\n"<< std:: endl;
    std:: cin >> T ;

    //std:: cout << "rc parameter\n"<< std:: endl;
    //std:: cin >> rc ;

    std:: cout << "binNo parameter\n"<< std:: endl;
    std:: cin >> binNo ;

    std:: cout << "kb parameter\n"<< std:: endl;
    std:: cin >> kb ;


    //infile >> file_name;   // store the last row 2nd column element in a string called file_name
    //read_coord(file_name, n, x, y, z);
};

void read_coord(char *file_name, int n, double *x, double *y, double *z){

  double C[n+3][4];
  char str[67];
  file_name = "initial_equilibrated_config.xyz";
  std ::string line;
  std:: ifstream infile;
   infile.open("initial_equilibrated_config.xyz");   // open the file
  if (!infile) {
        std::cout << "Unable to open file" <<std:: endl;
        exit(1); // terminate with error
  }
  int counter = 2;

  while(std::getline(infile, line)){   // reading the rest of lines from xyz file
      for (int j = 0; j < 4; j++){
        if (j == 0){
          infile >> str;   // reads the string of each row 1st column
          //std:: cout << str << "   \t   ";
        }
        else {
          infile >> C[counter][j];   // reads the 2nd, 3rd,  and 4th column of each row
        }
      }

      int j = 3;
      // incorporating the values inside x, y, z
      x[counter - 2] = C[counter][j-2];   // stores the elements of 2nd column of each row in the variable "x"
      y[counter - 2] = C[counter][j-1];   // stores the elements of 3rd column of each row in the variable "y"
      z[counter - 2] = C[counter][j];   // stores the elements of 4th column of each row in the variable "z"
      //std:: cout << x[counter - 2] << "   \t   " << y[counter - 2] << "   \t   " << z[counter - 2] << std:: endl;
      counter++;
  }
  infile.close();
};

double periodic_boundary_condition(double pair_distance, double box_length){   // Periodic boundary condition implemented
 if (pair_distance < - box_length/2){
  pair_distance = pair_distance + box_length;
 }
 else if (pair_distance > box_length/2){
  pair_distance = pair_distance - box_length;
 }
 return pair_distance;
};


int main(){
    int n= 999 ;
double *x = new double[n];
double *y = new double[n];
double *z = new double[n];


char str[67];
char *file_name = new char;  
//initial_conditions();
read_coord(file_name, n, x, y, z);
int binNo;
double K, box_length, binSize, acc = 0, epsilon, rho, kb, T , delta, sigma, rc;
rc = 2.5;
epsilon = 1.0;
sigma = 1.0 ;
box_length =  cbrt(n/rho);   // determination of Box length, L
binSize = (box_length*0.5)/double(binNo);   // binsize calculation
//std :: cout << rc << " \n   "<< std:: endl; 
//---------------
double energy_total = 0.0, E_stepwise;
  for (int i = 0; i < (n-1); i++){
    double W_stepwise = 0.0;
    for (int j = i+1; j < n; j++){
        double x_diff = periodic_boundary_condition((x[i] - x[j]), box_length);
        //std :: cout << x_diff << " \n   "<< std:: endl; 
        double y_diff = periodic_boundary_condition((y[i] - y[j]), box_length);
        double z_diff = periodic_boundary_condition((z[i] - z[j]), box_length);
        double R2 = (x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);   // calculating the distance between a pair
        //std :: cout << sigma << " \n   "<< epsilon<<  std:: endl; 
        
        //std :: cout << R2 << " \t   "<< pow(rc,2) <<  std:: endl; 

        if (R2 <= pow(rc, 2))   // condition for checking the pair distance within the cut-off distance
        {
          E_stepwise = 4*epsilon*(pow(sigma, 12)/pow(R2, 6) - pow(sigma, 6)/pow(R2, 3));  // calculation of Lennard-Jones Potential
         
        }

        else {
          E_stepwise = 0;
        }
    W_stepwise +=  E_stepwise;   // calculates the total energy for one selected particle
    //std :: cout << W_stepwise << " \n   "<< std:: endl;
    }
    //std :: cout << E_stepwise << " \t   "<<W_stepwise<< std:: endl;
    energy_total +=  W_stepwise;   // calculates the total energy of the system
    //std :: cout << energy_total << " \n   "<< std:: endl;  
  };
  

  return 0;

};







