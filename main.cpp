#include <fstream>
#include <iostream>
#include "compressed_matrix.hpp"

int main(){
    compressed_matrix m(5);
	m.read("mat.txt");
    std::cout<<"M = "<<std::endl;
    m.print_matrix();
    
    compressed_matrix L(5),U(5);
    // FIXME Incorrect result!!
    m.LU_decomposition(L,U);
    std::cout<<"L = "<<std::endl;
    L.print_matrix();

    std::cout<<"U = "<<std::endl;
    U.print_matrix();


}
