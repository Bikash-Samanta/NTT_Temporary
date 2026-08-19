#include <iostream>
#include <vector>
#include <cstdint>

const int32_t Q = 8380417;
const int N = 256;


uint8_t bit_rev8(uint8_t m) {
    uint8_t res = 0;
    for (int i = 0; i < 8; ++i) {
        if ((m >> i) & 1) {
            res |= (1 << (7 - i));
        }
    }
    return res;
}

int32_t power_mod(int64_t base, int64_t exp, int32_t mod) {
    int64_t res = 1;
    base = base % mod;
    
    while (exp > 0) {
        if (exp % 2 == 1) {
            res = (res * base) % mod;
        }
        base = (base * base) % mod;
        exp /= 2;
    }
    return static_cast<int32_t>(res);
}

std::vector<int32_t> generate_zetas() {
    const int32_t Q = 8380417;
    const int32_t ZETA = 1753; 
    
    std::vector<int32_t> zetas(256);
    
    // Calculate zetas[m] = (ZETA ^ BitRev8(m)) mod Q
    for (int m = 0; m < 256; ++m) {
        uint8_t reversed_m = bit_rev8(static_cast<uint8_t>(m));
        zetas[m] = power_mod(ZETA, reversed_m, Q);
    }
    
    return zetas;
}

static const std::vector<int32_t> ZETAS_TABLE = generate_zetas();

const std::vector<int32_t> zetas = ZETAS_TABLE; 

void dilithium_ntt(std::vector<int32_t>& w_hat, const std::vector<int32_t>& zetas) {
    if (w_hat.size() != N) return; 

    int m = 0;
    int len = 128; 

    while (len >= 1) {
        int start = 0;
        while (start < N) {
            m = m + 1;
            int32_t z = zetas[m]; 
            
            for (int j = start; j < start + len; ++j) {
                int64_t t = (static_cast<int64_t>(z) * w_hat[j + len]) % Q;
                w_hat[j + len] = (w_hat[j] - t + Q) % Q;
                w_hat[j] = (w_hat[j] + t) % Q;
            }
            start = start + 2 * len;
        }
        len = len / 2;
    }
}


void dilithium_intt(std::vector<int32_t>& w_hat, const std::vector<int32_t>& zetas) {
    if (w_hat.size() != 256) return;

    int k = 256;
    int len = 1;

    while (len < 256) {
        int start = 0;
        
        while (start < 256) {
            k = k - 1;
            
            int32_t z = (Q - zetas[k]) % Q; 
            for (int j = start; j < start + len; ++j) {
                int32_t t = w_hat[j];
                
                w_hat[j] = (t + w_hat[j + len]) % Q;
                
                int64_t diff = (t - w_hat[j + len] + Q) % Q;
                w_hat[j + len] = (diff * static_cast<int64_t>(z)) % Q;
            }
            start = start + 2 * len;
        }
        len = len * 2;
    }

    // Multiply every coefficient by 256^-1 mod Q
    // For Dilithium (q = 8380417), the inverse of 256 is 8347681
    const int64_t F = 8347681;
    for (int i = 0; i < 256; ++i) {
        w_hat[i] = (w_hat[i] * F) % Q;
    }
}

void print_poly(const std::string& label, const std::vector<int32_t>& poly, int count = 10) {
    std::cout << label << ": [";
    for (int i = 0; i < count; ++i) {
        std::cout << poly[i] << (i == count - 1 ? "" : ", ");
    }
    std::cout << " ...]" << std::endl;
}


int main() {
    std::cout << "Precomputing zetas table..." << std::endl;
    const std::vector<int32_t> zetas = generate_zetas();

    std::vector<int32_t> poly(256);
    for (int i = 0; i < 256; ++i) {
        poly[i] = i%12 + 12;
    }

    std::cout << "\n--- Original Polynomial ---" << std::endl;
    print_poly("Input", poly, 10);

    // Run Forward NTT
    dilithium_ntt(poly, zetas);
    std::cout << "\n--- After Forward NTT ---" << std::endl;
    print_poly("NTT  ", poly, 10);

    // Run Inverse NTT
    dilithium_intt(poly, zetas);
    std::cout << "\n--- After Inverse NTT ---" << std::endl;
    print_poly("INTT ", poly, 10);

    return 0;
}

/*
[akaza@akaza NTT]$ g++ -O3 NTT_Dilithium.cpp -o NTT_Dilithium
[akaza@akaza NTT]$ ./NTT_Dilithium
Precomputing zetas table...

--- Original Polynomial ---
Input: [12, 13, 14, 15, 16, 17, 18, 19, 20, 21 ...]

--- After Forward NTT ---
NTT  : [7892715, 1191699, 5320725, 6790807, 554247, 5587440, 12356, 371242, 8175844, 1053012 ...]

--- After Inverse NTT ---
INTT : [12, 13, 14, 15, 16, 17, 18, 19, 20, 21 ...]
*/