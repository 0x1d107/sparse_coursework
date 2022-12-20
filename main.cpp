#include <fstream>
#include <iostream>
#include "compressed_matrix.hpp"
#include <ctime>



int main(){
    compressed_matrix m(5);
    compressed_matrix L(5);
    compressed_matrix U(5);
    int N = 1000;
    std::vector<double> x(N,0),y(N),b(N);
    double q = 0.20;

    m = compressed_matrix(N);
    b = std::vector<double>(N,1);
    m.generate(1338,q);
    L = compressed_matrix(N);
    U = compressed_matrix(N);
    auto lu_start = clock();
    m.LU_decomposition(L,U);
    L.solve_L(b,y);
    U.solve_U(y,x);
    auto lu_end = clock();
    std::cout<<"LU time:"<<(lu_end - lu_start*1.0)/CLOCKS_PER_SEC<<std::endl;
    if(m.check(b,x))
        std::cout<<"[OK]"<<std::endl;
    else
        std::cout<<"[FAIL]"<<std::endl;
    x = std::vector<double>(N,0);
    auto bcg_start = clock();
    double r = m.BiCGStab_solve(b,x,10000);
    auto bcg_end = clock();
    std::cout<<"BiCGStab time:"<<(bcg_end - bcg_start*1.0)/CLOCKS_PER_SEC<<std::endl;
    std::cout<<"Residue norm:"<<r<<std::endl;
    if(m.check(b,x))
        std::cout<<"[OK]"<<std::endl;
    else
        std::cout<<"[FAIL]"<<std::endl;

}
