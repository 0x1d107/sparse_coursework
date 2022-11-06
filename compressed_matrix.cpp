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
    std::cout<<"[";
	for(int i = 0;i<row_num();i++){
		auto begin = rows[i];
		auto end = rows[i+1];
        std::cout<<"[";
		for(int j = 0; j< row_num();j++){	
			double a = 0;
			if(begin<end&&begin<row_num()*row_num()&&j == cols[begin]){
				a = value[begin];
				begin++;
			}


			std::cout<<a<<", ";
		}
		std::cout<<"],"<<std::endl;
	}
    std::cout<<"]";
}
void compressed_matrix::set_init(int i, int j, double v){
    rows[i] = std::min(rows[i],(int)value.size());
    if(cols.size()&&j == cols.back()){
        value.back() = v;
        return;
    }
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
double compressed_matrix::get(int i,int j) const{
    if(rows[i]>=cols.size())
        return 0.0;
    auto it = std::lower_bound(cols.begin()+rows[i],cols.begin()+std::min(rows[i+1],(int)cols.size()),j);
    if(it!=cols.end()&&*it == j){
        return value[it - cols.begin()];
    }
    return 0.0;
}
double compressed_matrix::mul_sub(const compressed_matrix& other,int i,int j,int k){
    double sum = 0.0;
    /*
    for(int ri=rows[i];ri<rows[i+1];ri++){
        std::cout<<"ri = "<<ri<<std::endl;
        if(ri>=cols.size())
            continue;
        int rk = cols[ri];
        std::cout<<"rk = "<<rk<<std::endl;
        if(rk>=k)
            break;
        double rv = value[rk];
        double ov = other.get(rk,j);
        sum+=rv*ov;
        std::cout<<"("<<ri<<","<<rk<<") "<<rv<<" * "<<ov<<std::endl;


    }*/
    for(int t=0;t<k;t++){
        auto a = get(i,t);
        auto b = other.get(t,j);
        std::cout<<"("<<a<<"*"<<b<<":"<<t<<") +" ;
        sum+=a*b;
    };
    std::cout<<std::endl;
    return sum;
}
void compressed_matrix::LU_decomposition(compressed_matrix &L, compressed_matrix &U){
    L.set_size(row_num());
    U.set_size(row_num());
    for(int i=0;i<row_num();i++){
        for(int j=0;j<row_num();j++){
            if(i<=j){
                //if(i==j)
                //    L.set_init(i,j,1);
                double t = (get(i,j)-L.mul_sub(U,i,j,i));
                std::cout<< "U["<<i<<","<<j<<"]::"<<t<<std::endl;
                U.set_init(i,j,t);
                //U.print_matrix();
                //if(i==j)
                //    L.set_init(i,j,1);
                
            }else{
                std::cout<<"M[i,j] = "<<get(i,j)<<std::endl;
                double t = (get(i,j)-L.mul_sub(U,i,j,j))/U.get(j,j);
                std::cout<< "L["<<i<<","<<j<<"]::"<<t<<"| U[j,j] ="<<U.get(j,j)<<std::endl;
                L.set_init(i,j,t);
                //L.print_matrix();
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

