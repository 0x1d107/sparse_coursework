#include <fstream>
#include <iostream>
#include "compressed_matrix.hpp"
#include <ctime>


int main(){
    compressed_matrix m(5);
    compressed_matrix L(5);
    compressed_matrix U(5);
    const int N = 500;
    std::vector<double> x(N,0),y(N),b(N);
	/*
    m.read("mat.txt");
    std::cout<<"M = "<<std::endl;
    m.print_matrix();
    std::cout<<std::endl;
    for(int i=0;i<5;i++){
        for(int j=0;j<5;j++){
            std::cout<<m.get(i,j)<<' ';
        }
        std::cout<<std::endl;
    }
    m.LU_decomposition(L, U);
    std::cout<<"L:\n";
    L.print_matrix();
    std::cout<<std::endl;
    std::cout<<"U:\n";
    U.print_matrix();
    std::cout<<std::endl;
    for(int i=0;i<m.row_num();i++)
        std::cin>>b[i];

    L.solve_L(b,y);
    U.solve_U(y,x);

    for(int i=0;i<m.row_num();i++){
        std::cout<< x[i]<<std::endl;
    }
    std::cout << "BiCStab"<<std::endl;
    x = std::vector<double>{1,1,1,1,1};
    m.BiCGStab_solve(b,x);
    for(int i=0;i<m.row_num();i++){
        std::cout<< x[i]<<std::endl;
    }*/
    std::cout<<"Generate matrix"<<std::endl;
    m = compressed_matrix(N);
    b = std::vector<double>(N,1);
    m.generate(1338,0.20);
    L = compressed_matrix(N);
    U = compressed_matrix(N);
    auto lu_start = clock();
    m.LU_decomposition(L,U);
    L.solve_L(b,y);
    U.solve_U(y,x);
    auto lu_end = clock();
    std::cout<<"LU: "<<(lu_end - lu_start*1.0)/CLOCKS_PER_SEC<<std::endl;
    for(int i=0;i<m.row_num();i++){
        //std::cout<< x[i]<<std::endl;
    }
    x = std::vector<double>(N,0);
    auto bcg_start = clock();
    double r = m.BiCGStab_solve(b,x,100000);
    auto bcg_end = clock();
    std::cout<<"ERR:"<<r<<std::endl;
    std::cout<<"BiCGStab: "<<(bcg_end - bcg_start*1.0)/CLOCKS_PER_SEC<<std::endl;
    for(int i=0;i<m.row_num();i++){
        //std::cout<< x[i]<<std::endl;
    }

}
