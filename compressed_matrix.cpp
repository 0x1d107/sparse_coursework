#include "compressed_matrix.hpp"
#include <limits.h>
#include <utility>
#include <fstream>
#include <iostream>
#include <algorithm>
compressed_matrix::compressed_matrix(int n){
    set_size(n);
}
void compressed_matrix::set_size(int n){
	rows.clear();
	rows.resize(n+1,n*n);
    rows[0] = 0;
	cols.clear();
	value.clear();
}
int compressed_matrix::elem_num(){
    return value.size();
}
int compressed_matrix::row_num(){
    return rows.size()-1;
}

void compressed_matrix::print_matrix(){
	for(int i = 0;i<row_num();i++){
		auto begin = rows[i];
		auto end = rows[i+1];
		for(int j = 0; j< row_num();j++){	
			int a = 0;
			if(begin<end&&j == cols[begin]){
				a = value[begin];
				begin++;
			}


			std::cout<<a<<' ';
		}
		std::cout<<std::endl;
	}
}
void compressed_matrix::set_init(int i, int j, double v){
    rows[i] = std::min(rows[i],(int)value.size());
	value.push_back(v);
	cols.push_back(j);
}
void compressed_matrix::read(const std::string &filename){
	using namespace std;
	ifstream mat_file(filename);
	int n,k;
	mat_file>> k>> n;
	set_size(k);
	for(int i=0;i<n;i++){
		int ii,jj,v;
		mat_file >> ii>>jj>>v;
		set_init(ii,jj,v);
	}
}
double compressed_matrix::get(int i,int j){
    if(rows[i]>=cols.size())
        return 0.0;
    auto it = std::lower_bound(cols.begin()+rows[i],cols.begin()+std::min(rows[i+1],(int)cols.size()),j);
    if(*it == j){
        return value[it - cols.begin()];
    }
    return 0.0;
}
double compressed_matrix::mul_sub(const compressed_matrix& other,int i,int j,int k){
    double sum = 0.0;
    for(int ri=rows[i];ri<rows[i+1];ri++){
        if(ri>=cols.size())
            break;
        int rk = cols[ri];
        if(rk>=k)
            break;
        int rv = value[rk];
        auto it = std::lower_bound(other.cols.begin()+other.rows[rk],other.cols.begin()+other.rows[rk+1],j);
        if(*it == j){
            sum+=rv*other.value[it-other.cols.begin()];
        }

    }
    return sum;
}
void compressed_matrix::LU_decomposition(compressed_matrix &L, compressed_matrix &U){
    L.set_size(row_num());
    U.set_size(row_num());
    for(int i=0;i<row_num();i++){
        for(int j=0;j<row_num();j++){
            if(i<=j){
                U.set_init(i,j,get(i,j)-L.mul_sub(U,i,j,i-1));
                if(i == j)
                    L.set_init(i,i,1);

            }else{
                L.set_init(i,j,(get(i,j)-L.mul_sub(U,i,j,j-1))/U.get(j,j));
            }
            
        }
    }
    
    
}
void compressed_matrix::solve_L(const std::vector<double> &b, std::vector<double> &y){

}
void compressed_matrix::solve_U(const std::vector<double> &y, std::vector<double> &x){

}


compressed_matrix::~compressed_matrix(){

}
