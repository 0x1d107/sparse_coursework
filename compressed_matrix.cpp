#include "compressed_matrix.hpp"
#include <limits.h>
#include <utility>
#include <fstream>
#include <iostream>
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

compressed_matrix::~compressed_matrix(){
}
