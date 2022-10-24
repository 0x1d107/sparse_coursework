#pragma once
#include <vector>
class compressed_matrix{
    protected:
        std::vector<int> rows;
        std::vector<int> cols;
        std::vector<double> value;
    public:
        compressed_matrix(int n);
        virtual ~compressed_matrix();
        int row_num();
        int elem_num();
        void print_matrix();
        void set(int i,int j, double v);
};
