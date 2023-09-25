//Laura Sofía Cortes 20202107049
// María Alejandra Moreno Ramírez 20202107023 
//--Fisica Computacional II---

//Include the libraries 
#include <iostream>
#include <iomanip>
#include <ios>
#include <cmath>
#include <vector>
#include <fstream>

using namespace std;

class Carro{
//Def the atributes 
    public:
        double  v,t_0,t=0.0;
        double x_0, x= 0.0;
        string  a ;
        vector <double> p, tmp;
       
    
//Method for get the position and storage in a vector 
    void  position ( double x_0 ,double t_0, double v){
        cout<< "Tiempo\tPosicion\n";
        
        for (double i=0.0;i<=25.0; i=i+0.1){
            x = x_0+ (v*(i+t_0));
            
            cout<<t<<" "<<x<<endl;
            t= i+t_0;
            //setpresicion(4);
            p.push_back(x);
            tmp.push_back(t);
            
            }
        }
    //Write the method for storage the p and tmp array 
    void write_files (string a) {
    string c= a+ ".txt";
    ofstream file(c);
        
    	 for (size_t j = 0; j < tmp.size(); ++j) {
        file << tmp.at(j) << "\t"<<p.at(j)<< "\n";        	
    };
    file.close();
    }; 




// The builder 

    Carro (){
            
        cout << "Object ID (Use char r or u for plot )"<<endl;
        cin >> a;
        cout << "Set the initial conditions "<<endl;
        cout << "Initial position x_0"<<endl;
        cin >> x_0;
        cout << "Initial time t_0"<<endl;
        cin >> t_0;
        cout << "Constant velocity"<<endl;
        cin >> v;            
        }   
 };

int main (){
Carro w= Carro();
cout<<" --------------------------"<<endl;
w.position(w.x_0,w.t_0,w.v);
cout<<" --------------------------"<<endl;
w.write_files(w.a);


Carro h= Carro();
cout<<" --------------------------"<<endl;
h.position(h.x_0,h.t_0,h.v);
cout<<" --------------------------"<<endl;
h.write_files(h.a);

//Make the loops for  get the position where the objects meet each other 
if(w.v == h.v){
    cout<< "Objects won't ever meet each other"<<endl;
}else{
    for (int k=0 ;k<= w.tmp.size(); k++){
        if (w.p[k]== h.p[k]){
            cout<< "The objects will meet each other at " <<w.p[k] << " meters at "<< h.tmp[k]<<" seconds"<<endl;
        }else{
            cout<< "Objects won't ever meet each other in this time inverval [250 seconds :'(]"<<endl;
        break;
        };
    };

};

//Script for plot the simulation 
    ofstream gnuarchivo("grafica.gp");
     gnuarchivo << "set title 'Posición vs. Tiempo'\n";
     gnuarchivo << "set xlabel 'Tiempo'\n";
     gnuarchivo << "set ylabel 'Posicion'\n";
     gnuarchivo << "plot 'r.txt' with lines title 'Posicion_r'\n\n";
     gnuarchivo << "replot 'u.txt' with lines title 'Posicion_u'\n\n";
     gnuarchivo << "set terminal png\n";
     gnuarchivo << "set output 'posicion.png'\n";
     gnuarchivo << "pause -1";
     gnuarchivo.close();
     system("gnuplot grafica.gp");





return 0;
};
