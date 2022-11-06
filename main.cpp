#include <fstream>
#include <iostream>
#include "compressed_matrix.hpp"

int main(){
    compressed_matrix m(5);
	m.read("mat.txt");
    std::cout<<"M = "<<std::endl;
    m.print_matrix();
    std::cout<<std::endl;
    compressed_matrix L(5);
    compressed_matrix U(5);
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

}
