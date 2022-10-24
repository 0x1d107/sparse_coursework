#include <fstream>
#include <iostream>
#include "compressed_matrix.hpp"

int main(){
    compressed_matrix m(5);
	m.read("mat.txt");
	m.print_matrix();

}
