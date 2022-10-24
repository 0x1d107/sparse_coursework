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
    public:
        compressed_matrix(int n);
		void read( const std::string & filename);
        virtual ~compressed_matrix();
        int row_num();
        int elem_num();
        void print_matrix();
		void LU_decomposition(compressed_matrix &L,compressed_matrix &U);
		void solve_L(const std::vector<double> &b, std::vector<double>& y);
		void solve_U(const std::vector<double> &y, std::vector<double> &x);
};
