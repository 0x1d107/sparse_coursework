#pragma once
#include <vector>
#include <string>
class compressed_matrix{
    protected:
        std::vector<int> rows;
        std::vector<int> cols;
        std::vector<double> value;
        void set_init(int i,int j, double v);
		void set_size(int n);
        virtual double mul_sub(const compressed_matrix& other,int i,int j,int k);
    public:
        compressed_matrix(int n);
		void read( const std::string & filename);
        void generate(int seed,double chance=0.5);
        virtual ~compressed_matrix();
        int row_num() const;
        int elem_num() const;
        void print_matrix();
        double get(int i,int j) const ;
		void LU_decomposition(compressed_matrix &L,compressed_matrix &U);
        void ILU_decomposition(compressed_matrix &L,compressed_matrix &U);
		void solve_L(const std::vector<double> &b, std::vector<double>& y);
		void solve_U(const std::vector<double> &y, std::vector<double> &x);
        std::vector<double> operator*(const std::vector<double>& vec)const ;
        
        std::vector<double> T_prod(const std::vector<double>& vec)const ;
		double BiCGStab_solve(const std::vector<double> &b, std::vector<double>& x,int n = 1000);
	bool check(std::vector<double> b, std::vector<double> x);
};
