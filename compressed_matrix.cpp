#include "compressed_matrix.hpp"
#include <limits.h>
#include <utility>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstdlib>
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
int compressed_matrix::elem_num()const{
    return value.size();
}
int compressed_matrix::row_num()const{
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
    auto it = std::lower_bound(cols.begin()+rows[i],
            cols.begin()+std::min(rows[i+1],(int)cols.size()),j);
    if(it!=cols.end()&&*it == j){
        return value[it - cols.begin()];
    }
    return 0.0;
}
double compressed_matrix::mul_sub(const compressed_matrix& other,int i,int j,int k){
    double sum = 0.0;
    for(int t=0;t<k;t++){
        auto a = get(i,t);
        auto b = other.get(t,j);
        sum+=a*b;
    };
    return sum;
}
void compressed_matrix::LU_decomposition(compressed_matrix &L, compressed_matrix &U){
    L.set_size(row_num());
    U.set_size(row_num());
    for(int i=0;i<row_num();i++){
        for(int j=0;j<row_num();j++){
            if(i<=j){
                double t = (get(i,j)-L.mul_sub(U,i,j,i));
                U.set_init(i,j,t);
                if(i==j)
                    L.set_init(i,j,1);
                
            }else{
                double t = (get(i,j)-L.mul_sub(U,i,j,j))/U.get(j,j);
                L.set_init(i,j,t);
            }
        }
    }
    
    
}
void compressed_matrix::ILU_decomposition(compressed_matrix &L, compressed_matrix &U){
    L.set_size(row_num());
    U.set_size(row_num());
    for(int i=0;i<row_num();i++){
        bool diag = false;
        for(int r=rows[i];r<std::min(rows[i+1],(int)cols.size());r++){
            int j = cols[r];
            double v= value[r];
            if(i<=j){
                double t = (v-L.mul_sub(U,i,j,i));
                U.set_init(i,j,t);
                if(!diag){
                    L.set_init(i,i,1);
                    diag = true;
                }
                
            }else{
                double t = (v-L.mul_sub(U,i,j,j))/U.get(j,j);
                L.set_init(i,j,t);
            }
        }
        if(!diag){
            L.set_init(i,i,1);
            diag = true;
        }
    }
}
void compressed_matrix::solve_L(const std::vector<double> &b, std::vector<double> &y){
    y = std::vector<double>(b);
    for(int i=0;i<row_num();i++){
        y[i]/=get(i,i);
        for(int j=i+1;j<row_num();j++){
            

            y[j]-=get(j,i)*y[i];
        }
    }
}
void compressed_matrix::solve_U(const std::vector<double> &y, std::vector<double> &x){
    x = std::vector<double>(y);
    for(int i=row_num()-1;i>=0;i--){
        x[i]/=get(i,i);
        for(int j=i-1;j>=0;j--){
            x[j]-=get(j,i)*x[i];
        }
    }

}
std::vector<double> compressed_matrix::operator*(const std::vector<double>& vec)const {
    std::vector<double> res(row_num(),0);
    for(int i=0;i<res.size();i++){
        for(int j=rows[i];j<std::min(rows[i+1],(int)cols.size());j++){
            res[i]+= value[j] * vec[cols[j]];
        } 
    }
    return res;

}
std::vector<double> operator+(const std::vector<double> &a, const std::vector<double> &b){
    std::vector<double> c(a.size());
    for(int i=0;i<a.size();i++){
        c[i]=a[i]+b[i];
    }
    return c;
}
std::vector<double> operator*(double a, const std::vector<double> &b){
    std::vector<double> res(b);
    for(double &bi : res)
        bi*=a;
    return res;
}
std::vector<double> operator-(const std::vector<double> &a, const std::vector<double> &b){
    std::vector<double> c(a.size());
    for(int i=0;i<a.size();i++){
        c[i]=a[i]-b[i];
    }
    return c;
}

double operator*(const std::vector<double> &a, const std::vector<double> &b){
    double sum=0.0;
    for(int i=0;i<a.size();i++){
        sum+=a[i]*b[i];
    }
    return sum;
}
std::vector<double> operator/(std::vector<double> &a, double b){
    std::vector<double> c(a.size());
    for(int i=0;i<a.size();i++){
        c[i]=a[i]/b;
    }
    return c;
}

double norm(const std::vector<double> &x){
    double sum = 0;
    for(double s:x)
        sum+=s*s;
    return sqrt(sum);
}

std::vector<double> compressed_matrix::T_prod(const std::vector<double> &vec)const {
    std::vector<double> res(vec.size(),0);
    for(int i=0;i<vec.size();i++){
        for(int k =0;k<row_num();k++){
            res[i]+=get(k,i)*vec[k];
            }
        }
    return res;
}

double compressed_matrix::BiCGStab_solve(const std::vector<double> &b, std::vector<double>& x,int m){
    const compressed_matrix &A = *this;
    std::vector<double> r0 = b - A*x;
    std::vector<double> r0_ = r0;
    std::vector<double> p0 = r0;
    std::vector<double> p0_ = r0;
    for(int j=0;j<m;j++){
        double alpha = (r0*r0_)/((A*p0) * r0_);
        std::vector<double> s = r0 - alpha * (A*p0);
        double omega = (A * s) * s/((A*s)*(A*s));
        x = x + alpha*p0 +omega * s;
        std::vector<double> r1 = s - omega * (A *s);
        double beta = (r1*r0_)/(r0*r0_)*alpha/omega;
        p0 = r1 + beta * (p0 - omega*(A*p0));
        r0 = r1;
        double rn = norm(r0);
        if(rn<1e-5){
            return rn;
        }
    }
    return norm(r0);
}
double rnd(){
    return rand()*2.0/RAND_MAX - 1.0;
}
void compressed_matrix::generate(int seed,double chance){
    srand(seed);
    for(int i=0;i<row_num();i++){
        double sum = fabs(rnd())*10;
        for(int j=0;j<row_num();j++){
            if(i==j)
                set_init(i,j,sum*2);
            else if(j<i){
                set_init(i,j,get(j,i));
            }else{
                if(sum>0&&fabs(rnd())<chance){
                    double val = sum*rnd();
                    sum-=fabs(val);
                    set_init(i,j,val);
                }
            }
        }
    }
    
}

bool compressed_matrix::check(std::vector<double> b, std::vector<double> x){
    std::vector<double> r = (*this)*x - b;
    for(int i=0;i<r.size();i++){
        if(std::abs(r[i])>=1e-5)
            return false;
    }
    return true;
}


compressed_matrix::~compressed_matrix(){
}

