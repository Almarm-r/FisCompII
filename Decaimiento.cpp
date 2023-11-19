#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>  // Necesario para std::count
#include <ctime>      // Necesario para std::time
#include <cstdlib>    // Necesario para std::rand y std::srand
#include <fstream>
#include <random>
//Pasar todas las condiciones del problema original 
const double t_half_rad   = 14.8;
 const double t_half_act   = 10.0;
int N0           = 1000;
double t1           = 100.0;
int n_timepoints = 1000; // PARA EL MÉTODO NUMÉRICO
//número de átomos iniciales del radón
int N0_rad = 1000;
int N0_act = 0;

const double l1 = std::log(2) / t_half_rad;
const double l2 = std::log(2) / t_half_act;
    //Semilla para numeros aleatorios utilizando la del porblema original 1
std::mt19937 generator(1); 
    //Struc para organizar juntos todos los datos que incluiremos en el metodo MC
struct Resultados_ProbablesMonteCarlo {
    std::vector<double> count_radium;
    std::vector<double> count_actinium;
    std::vector<double>  count_beta;
    std::vector<double>  rand_numb;
    std::vector<double> atoms;
    double p_decay_rad;
    double p_decay_act;
};
//Para las familias de soluciones análiticas
double analytic_rad(int N0 ,double times){
    
    return N0*std::exp(-times /t_half_rad * std::log(2));
};

double analytic_ac(int N0_act, double times){

    double l1 = std::log(2)/t_half_rad;
    double l2 = std::log(2)/t_half_act;

    return ((l1/(l2-l1))*N0_rad *(std::exp(-times*l1)  -std::exp(-times*l2))+ N0_act* std::exp(-times*l2));

};


Resultados_ProbablesMonteCarlo montecarlo_method (int N0, double t1, int n_timepoints){

    double dt = t1 / n_timepoints;
//Arreglos a rellenar 

    std::vector<double> count_radium(n_timepoints, 0.0);
    std::vector<double> count_actinium(n_timepoints, 0.0);
    std::vector<double>  count_beta(n_timepoints, 0.0);
    std::vector<double>  rand_numb(n_timepoints, 0.0);
    std::vector<double> atoms(N0, 1.0);

     
    double p_decay_rad = 1.0 - std::exp(-dt / t_half_rad  * std::log(2.0));
    double p_decay_act = 1.0 - std::exp(-dt / t_half_act * std::log(2.0));
//Aquí las condiciones orginales del problema para el metodo MC  
    for (int idxTime = 0; idxTime < n_timepoints; ++idxTime) {
        // Contar los átomos de cada tipo
        count_radium[idxTime] = std::count(atoms.begin(), atoms.end(), 1.0);
        count_actinium[idxTime] = std::count(atoms.begin(), atoms.end(), 2.0);


        for (int idxAtom = 0; idxAtom < N0; ++idxAtom) {
            //Condiciones de la tabla
            if (atoms[idxAtom] == 1.0) {
                double a = static_cast<double>(rand()) / RAND_MAX;
                if (a <= p_decay_rad) {
                    atoms[idxAtom] = 2.0;
                } else {
                    atoms[idxAtom] = 1.0;
                }
            } else if (atoms[idxAtom] == 2.0) {
                double b = static_cast<double>(rand()) / RAND_MAX;
                if (b <= p_decay_act) {
                    atoms[idxAtom] = 3.0;
                } else {
                    atoms[idxAtom] = 2.0;
                }
            }
        }
    


}
return Resultados_ProbablesMonteCarlo{count_radium,count_actinium, atoms,  count_beta};
};
//Las soluciones de las EDO 

double EDO_rad(double N_rad, double N_act) {
    return -l1* N_rad;
}
double EDO_act(double N_rad, double N_act) {
    return l1* N_rad - l2 * N_act;
}
//Runge Kutta orden 4 

void Runge_Kutta_Method(std::vector<double>& N_rad, std::vector<double>& N_act, double dt){
for (int i = 0; i < n_timepoints - 1; ++i) {
        double K1_rad = dt * EDO_rad(N_rad[i], N_act[i]);
        double K1_act = dt * EDO_act(N_rad[i], N_act[i]);

        double K2_rad = dt * EDO_rad(N_rad[i] + 0.5 * K1_rad, N_act[i] + 0.5 * K1_act);
        double K2_act = dt * EDO_act(N_rad[i] + 0.5 * K1_rad, N_act[i] + 0.5 * K1_act);

        double K3_rad = dt * EDO_rad(N_rad[i] + 0.5 * K2_rad, N_act[i] + 0.5 * K2_act);
        double K3_act = dt * EDO_act(N_rad[i] + 0.5 * K2_rad, N_act[i] + 0.5 * K2_act);

        double K4_rad = dt * EDO_rad(N_rad[i] + K3_rad, N_act[i] + K3_act);
        double K4_act = dt * EDO_act(N_rad[i] + K3_rad, N_act[i] + K3_act);

        N_rad[i + 1] = N_rad[i] + (K1_rad + 2 * K2_rad + 2 * K3_rad + K4_rad) / 6.0;
        N_act[i + 1] = N_act[i] + (K1_act + 2 * K2_act + 2 * K3_act + K4_act) / 6.0;
    }
};

int main (){
    // Inicializando condiciones del problema orginal 
    const double t1 = 100.0;
    const int n_timepoints = 1000;
    const double dt = t1 / n_timepoints;
    const double t_half_rad = 14.8;
    const double t_half_act = 10.0;
    const int N0 = 1000;
//Datos del  vector tiempo 
    std::vector<double> times;
    for (int i = 0; i < n_timepoints; ++i) {
        times.push_back(i * dt);
    }

    //Llamando las Soluciones analiticas
    std::vector<double> n_analytic_ra(n_timepoints);
    std::vector<double> n_analytic_ac(n_timepoints);
    for(int i = 0; i < times.size(); ++i){
        
        n_analytic_ra[i]= analytic_rad(N0, times[i]);                                            
        n_analytic_ac[i] = analytic_ac(0, times[i]);    
    };
    //Solucionrd numericas RUNGE_KUTTA 
    std::vector<double> N_rad_numeric(n_timepoints);
    std::vector<double> N_act_numeric(n_timepoints);
    //Aplicar las condciones iniciales. 
    N_rad_numeric[0] = N0_rad;
    N_act_numeric[0] = N0_act;
    Runge_Kutta_Method(N_rad_numeric, N_act_numeric, dt);



    //Solucion Montercarlo 

    Resultados_ProbablesMonteCarlo Simulation_Data = montecarlo_method(N0, t1, n_timepoints);
        //Definir la data en vectores apropiados
    std::vector<double>& n_rad = Simulation_Data.count_radium;
    std::vector<double>& n_act = Simulation_Data.count_actinium;



//Script de la grafica 
    std::ofstream plot("de_plot.gnu");
    plot << "set terminal wxt size 1000,850\n";
    plot << "set title 'Decaimiento de Radium y Actinium'\n";
    plot << "set xlabel 't'\n";
    plot << "set ylabel 'N'\n";
    plot << "set key left top\n";

    // Definir el estilo de línea punteada
    plot << "set style line 1 lt 1 lc rgb 'red' dashtype 2\n";
    plot << "set style line 2 lt 1 lc rgb 'orange' dashtype 2\n";
    plot << "set style line 3 lt 1 lc rgb 'brown'\n";
    plot << "set style line 4 lt 1 lc rgb 'black'\n";
    plot << "set style line 5 lt 1 lc rgb 'green'\n";
    plot << "set style line 6 lt 1 lc rgb 'blue'\n";

    // Graficar los resultados de Monte Carlo en una sola gráfica
    plot << "plot '-' with lines linestyle 3 title 'Montecarlo Radium', "
        << "'-' with lines linestyle 4 title 'Montecarlo Actinium', "
        << "'-' with lines linestyle 6 title 'Numeric Solution Actinium', "
        << "'-' with lines linestyle 5 title 'Numeric Solution Radium', "
        << "'-' with lines linestyle 2 title 'Analytical Solution Radium', "
        << "'-' with lines linestyle 1 title 'Analytical Solution Actinium'\n";
 //Obtener los datos para cada metodo    
    // Datos de Monte Carlo para Radium
    for (int i = 0; i < n_timepoints; ++i) {
        plot << i * dt << " " << n_rad[i] << "\n";
    }
    plot << "e\n";

    // Datos de Monte Carlo para Actinium
    for (int i = 0; i < n_timepoints; ++i) {
        plot << i * dt << " " << n_act[i] << "\n";
    }
    plot << "e\n";

    // Graficar los resultados numéricos
    for (int i = 0; i < n_timepoints; ++i) {
        plot << times[i] << " " << N_rad_numeric[i] << "\n";
    }
    plot << "e\n";
    for (int i = 0; i < n_timepoints; ++i) {
        plot << times[i] << " " << N_act_numeric[i] << "\n";
    }
    plot << "e\n";

    // Datos analíticos del actinio
    for (double t = 0; t <= t1; t += dt) {
        plot << t << " " << analytic_ac(N0_act, t) << "\n";
    }
    plot << "e\n";
    // Datos analíticos del Radium
    for (double t = 0; t <= t1; t += dt) {
        plot << t << " " << analytic_rad(N0_rad, t) << "\n";
    }
    plot << "e\n";


    plot << "pause -1";
    plot.close(); 

    std::system("gnuplot de_plot.gnu");

    return 0;
};