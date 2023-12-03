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

    std:: cout << "rc parameter\n"<< std:: endl;
    std:: cin >> rc ;

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

double system_total_energy(int n, double box_length, double epsilon, double sigma, double rc, double *x, double *y, double *z){   // function for calculating the total system energy
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
  return energy_total;
};
void rdf(int n, int binNo, double binSize, double box_length, double *count, double *x, double *y, double *z){
  double dR;
  for (int i = 0; i < (n-1); i++){
    for (int j = i+1; j < n; j++){
      double dx = periodic_boundary_condition((x[i] - x[j]), box_length);
      double dy = periodic_boundary_condition((y[i] - y[j]), box_length);
      double dz = periodic_boundary_condition((z[i] - z[j]), box_length);
      dR = pow((dx*dx + dy*dy + dz*dz), 0.5);   // calculating the magnitude of each pair-wise distances
      double initial = 0.0;
      for (int k = 0; k < binNo; k++){
        if (dR <= box_length/2 && dR > initial && dR <= (initial + binSize)){
          count[k] = count[k] + 1.0;
        }
        initial = initial + binSize;
      }
    }
  }
};



int main(){
//---------------------------------------
int n= 999 ;
double *x = new double[n];
double *y = new double[n];
double *z = new double[n];
char str[67];
char *file_name = new char;  
//initial_conditions();
//------------------------------------------
double K, box_length, binSize, acc = 0, epsilon, rho, kb, T , delta, sigma, rc;
epsilon   =1.0;
rho       =0.85;
binNo     =500;
kb        =1.0;
T         =0.723;
delta     =0.2;
sigma     =1.0;
rc        =2.5;

read_coord(file_name, n, x, y, z);


//----------------------------------------
int mcMAX = 1000, m = 100 , binNo;
srand(time(0));
double *count = new double[binNo]; 
 for (int i = 0; i < binNo; i++){
     count[i] = 0;
  };
box_length =  cbrt(n/rho);   // determination of Box length, L
binSize = (box_length*0.5)/double(binNo);   // binsize calculation
 
  //Opening the files to write
  std::ofstream outffile;
  outffile.open("coordinates.dat");
  std::ofstream outfile;
  outfile.open("energy_vs_mc_time_step_0.2_temp_0.723.dat");
  std::ofstream outFile;
  outFile.open("acceptance_ratio_0.2_temp_0.723.dat");
  std:: ofstream outFfile;
  outFfile.open("rdf_0.2_temp_0.723.dat");

  
  // calculating the total system energy before monte carlo move
  double energy_before_disp = system_total_energy(n, box_length, epsilon, sigma, rc, x, y, z);   // calculating the total system energy before monte carlo move
  outfile << "0" << " \t   " << energy_before_disp << std:: endl; 
    //std:: cout << "0" << " \t   " << energy_before_disp << std:: endl;   // print the initial total energy in the provided output file
    //std:: cout << "0" << " \t   " <<  energy_before_disp << std:: endl;   // print the initial total energy in the provided output file

    // print the initial total energy in the provided output file
  // Initializing Monte Carlo Loops
  for (int mc = 1; mc < mcMAX; mc++){
    int k = rand() % 1000;   // choosing a random particle which I will take as basis particle.
    // It will change in each MC loop
    double ranf = double(rand())/double(RAND_MAX);
    // Random number of dx, dy, dz implemented to give x, y, z random and different moves
    double dx = delta*(ranf - 0.5);
    double dy = delta*(ranf - 0.5);
    double dz = delta*(ranf - 0.5);

    // update the x, y, z co-ordinates after the move of the k-th particle
    x[k] = x[k] + dx;
    y[k] = y[k] + dy;
    z[k] = z[k] + dz;

    double energy_after_disp = system_total_energy(n, box_length, epsilon, sigma, rc, x, y, z);   // calculating the total system energy after monte carlo move

    // Applying 'Metropolis Method'
    if (energy_after_disp <= energy_before_disp){   // accept the move
      outfile << mc << "   \t    " << energy_after_disp << "   \t   " << k << "   \t   " << "1" << std::endl;
//std:: cout << mc << "   \t    " << energy_after_disp << "   \t   " << k << "   \t   " << "1" << std::endl;

      // 'k' means which particle is being sampled in the current MC step
      // '1' means accepted and '0' means rejected
      acc = acc + 1.0;   // increasing 'acc' by 1
      // update the co-ordinates i.e. x, y, z values after adding with dx, dy, dz
      // Taking energy_after_disp into energy_before_disp to calculate next loop.
      // Because if the move is accepted then the energy_after_disp will be the energy_before_disp for the next loop
      energy_before_disp = energy_after_disp;
    }
    else {
      double random = double(rand())/double(RAND_MAX);   //generating a random no between 0 and 1
      double energy_diff = energy_after_disp - energy_before_disp;  // calculating energy differnce
      double P = exp(-(energy_diff)/(kb*T));   // Bolzman factor
      if (random < P){   // accept the move
        outfile << mc << "   \t    " << energy_after_disp << "   \t   " << k << "   \t   " << "1" << std::endl;   // writing in a provided output file
     //std:: cout  << mc << "   \t    " << energy_after_disp << "   \t   " << k << "   \t   " << "1" << std::endl;   // writing in a provided output file

        acc = acc + random;  // increasing 'acc' by random number
        // Taking energy_after_disp into energy_before_disp to calculate next loop.
        energy_before_disp = energy_after_disp;
      }

      else {   // reject the move
        outfile << mc << "   \t    " << energy_before_disp << "   \t   " << k << "   \t   " << "0" << std::endl;
        //std::cout << mc << "   \t    " << energy_before_disp << "   \t   " << k << "   \t   " << "0" << std::endl;

        // As the move is rejected then x, y,z values are being downgraded to the previous values before starting the current MC loop
        x[k] = x[k] - dx;
        y[k] = y[k] - dy;
        z[k] = z[k] - dz;

      };
    };

    if(mc % int(m) == 0){
      for (int i = 0; i < n; i++){
        outffile << x[i] << "   \t   " << y[i] << "   \t    " << z[i] << std:: endl;   // printing the coordinates of all particles
        //std:: cout << x[i] << "   \t   " << y[i] << "   \t    " << z[i] << std:: endl;   // printing the coordinates of all particles

      }
      outffile << "\n";
    };
        if(mc % int(m) == 0){
      rdf(n, binNo, binSize, box_length, count, x, y, z);
    };
 };
   // print the average accepatnce ratio
  outFile << acc/double(mcMAX) << std::endl;
  //std:: cout  << acc/double(mcMAX) << std::endl;

  double r[binNo], g[binNo];

  r[0] = 0.0;
  // Now I have count[1], .............. count[N] total values for all the MC steps. That's why I have to divide final g[i] with total number of MC steps
  for (int i = 0; i < binNo; i++){
    g[i] = count[i];
    g[i] = g[i]/(((4.0*M_PI)/3.0)*(pow((r[i]+binSize),3) - pow(r[i],3))*rho);  // Normalization of g(r)
    outFfile << r[i] << "   \t   " << (g[i]*2*m)/(n*mcMAX) << std::endl;
    r[i+1] = r[i] + binSize;
  };

  std:: cout << "System total energy :\t" << energy_before_disp<< std::endl; 


  // closing all the files
  outffile.close();
  outfile.close();
  outFile.close();
  outFfile.close();

  //Plotting results 

  std::ofstream gnu;
  std:: ofstream gnu_2;

  gnu.open("rdf.gnu");
  gnu << "set terminal jpeg enhanced font 'arial,10' size 800,600" << std:: endl;
gnu << "set output 'grafica_rdf.jpg'"<< std:: endl;
  gnu << "set title 'Radial Distribution function g(r)'"<< std:: endl;
  gnu<< "set xlabel 'r[A]'" << std::endl;
gnu << "set ylabel 'g(r)'"<< std::endl;
gnu<< "plot 'rdf_0.2_temp_0.723.dat' using 1:2 " << std:: endl;
gnu << "pause -1"<<std::endl;
gnu_2.open("coord.gnu");
gnu_2 << "set view 60,30,1.0,1.5" << std::endl;
gnu << "set output 'grafica_coord.jpg'"<< std:: endl;
gnu_2 << "set title 'coord'"<< std:: endl;
gnu_2 << "set xlabel 'x'" << std::endl; 
gnu_2 << "set ylabel 'y'" << std::endl; 
gnu_2 << "set zlabel 'z'" << std::endl;
gnu_2<< "splot 'coordinates.dat' using 1:2:3"<< std:: endl; 
gnu_2 << "pause -1"<<std::endl;

gnu.close();
gnu_2.close();
std::system("rdf.gnu");
std::system("coord.gnu");


      return 0;

};
