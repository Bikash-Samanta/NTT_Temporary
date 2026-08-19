#include <iostream>
#include <vector>

const int Q = 3329;
const int N = 256;
const int ZETA = 17;

int modQ(int a) {
    int r = a % Q;
    if (r < 0) {
        r += Q;
    }
    return r;
}

// (base^exp) % Q
int powerMod(int base, int exp) {
    int res = 1;
    base = modQ(base);
    while (exp > 0) {
        if (exp % 2 == 1) res = modQ(res * base);
        base = modQ(base * base);
        exp /= 2;
    }
    return res;
}

// a^(Q-2) = a^-1 mod Q
int modInverse(int a) {
    return powerMod(a, Q - 2);
}

// Reverse the bits of a 7-bit integer
int bitRev7(int v) {
    int res = 0;
    for (int i = 0; i < 7; ++i) {
        res |= ((v >> i) & 1) << (6 - i);
    }
    return res;
}


// In-place Kyber NTT (Cooley-Tukey)
void kyberNTT(std::vector<int>& f) {
    int i = 1; 
    
    for (int len = 128; len >= 2; len /= 2) {
        for (int start = 0; start < 256; start += 2 * len) {
            
            int zeta = powerMod(ZETA, bitRev7(i));
            i++;
            
            for (int j = start; j < start + len; ++j) {
                int t = modQ(zeta * f[j + len]);
                
                f[j + len] = modQ(f[j] - t);
                f[j] = modQ(f[j] + t);
            }
        }
    }
}

// In-place Kyber Inverse NTT (Gentleman-Sande)
void kyberINTT(std::vector<int>& f) {
    for (int len = 2; len <= 128; len *= 2) {
        for (int start = 0; start < 256; start += 2 * len) {
            
            int i = (128 / len) + (start / (2 * len));
            
            int zeta = powerMod(ZETA, bitRev7(i));
            int zeta_inv = modInverse(zeta);
            
            for (int j = start; j < start + len; ++j) {
                int t = f[j];
                
                f[j] = modQ(t + f[j + len]);
                f[j + len] = modQ(t - f[j + len]);
                f[j + len] = modQ(f[j + len] * zeta_inv);
            }
        }
    }
    
    // Multiply all coefficients by 128^-1 mod 3329
    // 128 * 3303 = 422784 = 1 mod 3329
    const int SCALE_FACTOR = 3303; 
    for (int j = 0; j < 256; ++j) {
        f[j] = modQ(f[j] * SCALE_FACTOR);
    }
}


void printSample(const std::vector<int>& f, int count = 8) {
    std::cout << "[ ";
    for (int i = 0; i < count; ++i) {
        std::cout << f[i] << (i == count - 1 ? "" : ", ");
    }
    std::cout << " ... ]\n";
}

int main() {
    // Initialize a polynomial 'f' with arbitrary values
    std::vector<int> original_f(N);
    for (int i = 0; i < N; ++i)
        original_f[i] = (i * i + 33) % Q;
    
    std::vector<int> f = original_f;

    std::cout << "=== KYBER NTT VERIFICATION ===\n\n";
    
    std::cout << "1. Original Polynomial (first 8 coeffs):\n";
    printSample(original_f);

    kyberNTT(f);
    std::cout << "\n2. After Forward NTT (first 8 coeffs):\n";
    printSample(f);

    kyberINTT(f);
    std::cout << "\n3. After Inverse NTT (INTT) (first 8 coeffs):\n";
    printSample(f);

    // Verify identical arrays
    bool success = true;
    for (int i = 0; i < N; ++i) {
        if (f[i] != original_f[i]) {
            success = false;
            std::cout << "\n[!] Mismatch at index " << i 
                      << ": Expected " << original_f[i] << ", got " << f[i] << "\n";
            break;
        }
    }

    std::cout << "\n====================================\n";
    if (success) {
        std::cout << "SUCCESS: INTT(NTT(f)) == f\n";
        std::cout << "The implementation perfectly reverses the polynomial.\n";
    } else {
        std::cout << "FAILURE: The arrays do not match.\n";
    }
    std::cout << "====================================\n";

    return 0;
}

/*
[akaza@akaza NTT]$ g++ -O3 NTT_Kyber.cpp -o NTT_Kyber
[akaza@akaza NTT]$ ./NTT_Kyber
=== KYBER NTT VERIFICATION ===

1. Original Polynomial (first 8 coeffs):
[ 33, 34, 37, 42, 49, 58, 69, 82 ... ]

2. After Forward NTT (first 8 coeffs):
[ 1702, 318, 1381, 2601, 1587, 1479, 83, 738 ... ]

3. After Inverse NTT (INTT) (first 8 coeffs):
[ 33, 34, 37, 42, 49, 58, 69, 82 ... ]

====================================
SUCCESS: INTT(NTT(f)) == f
The implementation perfectly reverses the polynomial.
====================================
*/