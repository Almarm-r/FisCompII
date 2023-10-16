#include <iostream>
#include <vector>
#include <random>

// Congruential Multiplicative PRNG
class CongruentialMultiplicative {
public:
    CongruentialMultiplicative(unsigned long long seed, unsigned long long a, unsigned long long m) : seed(seed), a(a), m(m) {}

    double generate() {
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
    unsigned long long seed = 12345; // Initial seed
    unsigned long long a = 1664525; // Multiplier
    unsigned long long m = 4294967296; // Modulus (2^32)

    CongruentialMultiplicative rng(seed, a, m);

    // Generate pseudorandom numbers and store them in a vector

    int num_samples;
    std:: cout<< "Introduce the size for the random serie" << std:: endl;
    std:: cin >> num_samples;
    std::vector<double> random_numbers;
    for (int i = 0; i < num_samples; i++) {
        random_numbers.push_back(rng.generate());
        std:: cout << random_numbers[i]<< std:: endl;
    }

    return 0;
}