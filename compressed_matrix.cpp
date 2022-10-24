#include "compressed_matrix.hpp"
compressed_matrix::compressed_matrix(int n){
    rows.resize(n+1,0);
    rows[0] = 0;
    
}
int compressed_matrix::elem_num(){
    return value.size();
}
int compressed_matrix::row_num(){
    return rows.size();
}

void compressed_matrix::print_matrix(){
}
void compressed_matrix::set(int i, int j, double v){
    
}

compressed_matrix::~compressed_matrix(){
}
