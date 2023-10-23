#include <iostream>
#include <vector>
#include <random>
#include<cmath>
#include<fstream>

// Congruential Multiplicative PRNG
class CongruentialMethod {
public:
    CongruentialMethod(unsigned long long seed, unsigned long long a, unsigned long long m) : seed(seed), a(a), m(m) {}

    double generateS() {
        seed = (a * seed) % m;
        return static_cast<double>(seed) / static_cast<double>(m);
    }

private:
    unsigned long long seed;
    unsigned long long a;
    unsigned long long m;
};

int main() {
    // Initialize PRNG with appropriate parameters (Builder)
    unsigned long long seed = 281; // Initial seed
    unsigned long long a = pow(7,5); // Sumand 
    unsigned long long m = pow(2,31)-1; // Modulus (2^32-1)
    unsigned long long m_2= pow(2,10); // Modulus (2^32-1)

    CongruentialMethod rng(seed, a, m);
    CongruentialMethod rng_2(seed, a, m_2);

    // Generate pseudorandom numbers and store them in a vector
    std:: ofstream file("Table.txt");
    std::ofstream gnu("Grafica_R.txt");
 std::ofstream  gnu_2("Grafica_C.txt");
 int k =0; 
file << "pi" << " \t"<< "err" << " \t"<< "pi_2" << " \t"<< "err_2 "<< std:: endl; 


for (int k=0; k<=6; k++ ){
    double num_samples;
    std:: cout<< "Introduce the size for the random serie " << std:: endl;
    std:: cin >> num_samples;



    std::vector<double> random_numbersx , random_numbersy,random_numbersx_2 , random_numbersy_2;
    for (int i = 0; i < num_samples; i++) {
        random_numbersx.push_back(rng.generateS());
        random_numbersy.push_back(rng.generateS());
        std:: cout << random_numbersx[i] <<"  \t"<< random_numbersy[i]<<std:: endl;
        random_numbersx_2.push_back(rng_2.generateS());
        random_numbersy_2.push_back(rng_2.generateS());
       
        gnu << random_numbersx_2[i] <<"  \t"<< random_numbersy_2[i]<<std:: endl;
        
    };
// Contar los puntos dentro del circulo 
double r , r_2;
double contp =0 , contp_2= 0;
for (int j =0 ; j<= random_numbersx.size(); j++) {
    r = sqrt((random_numbersx[j]*random_numbersx[j])+(random_numbersy[j]*random_numbersy[j]));
    r_2 = sqrt((random_numbersx_2[j]*random_numbersx_2[j])+(random_numbersy_2[j]*random_numbersy_2[j]));

    if (r<=1){
        
        contp+=1;

    };
    if (r_2<=1){
        contp_2+=1;
        
        gnu_2 << random_numbersx_2[j] <<"  \t"<< random_numbersy_2[j]<<std:: endl;
    };
//std:: cout << "------------------------------------------------------------" << std:: endl;
//std:: cout << r << std:: endl;
};
std:: cout << "----------------------NUMERO DE PUNTOS EN EL CIRCULO ---------------------------------" << std:: endl;
std:: cout << contp << std:: endl;
std:: cout << contp_2 << std:: endl;

double pi, pi_2, err, err_2;
pi = 4* contp/num_samples;
pi_2 = 4* contp_2/num_samples;

err = pi - M_PI;
err_2 = pi_2 - M_PI;
std:: cout << "---------------------PI----------------------------" << std:: endl;
std:: cout << pi << std:: endl;
std:: cout << pi_2 << std:: endl;
file << pi << "\t"<< err << "\t"<< pi_2 << "\t"<< err_2 << "\n"; 

};
file.close(), gnu.close(), gnu_2.close();
    return 0;
}