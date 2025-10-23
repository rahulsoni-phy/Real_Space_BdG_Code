#include <vector>
#include <complex>
using namespace std;
#ifndef tensors
#define tensors

typedef vector< complex<double> >  Mat_1_Complex_doub;
typedef vector<Mat_1_Complex_doub> Mat_2_Complex_doub;
typedef vector<Mat_2_Complex_doub> Mat_3_Complex_doub;
typedef vector<Mat_3_Complex_doub> Mat_4_Complex_doub;
typedef vector<Mat_4_Complex_doub> Mat_5_Complex_doub;

typedef vector<double> Mat_1_doub;
typedef vector<Mat_1_doub> Mat_2_doub;
typedef vector<Mat_2_doub> Mat_3_doub;
typedef vector<Mat_3_doub> Mat_4_doub;
typedef vector<Mat_4_doub> Mat_5_doub;

typedef vector<int> Mat_1_int;
typedef vector<Mat_1_int> Mat_2_int;
typedef vector<Mat_2_int> Mat_3_int;
typedef vector<Mat_3_int> Mat_4_int;

typedef vector<string> Mat_1_string;
typedef vector<Mat_1_string> Mat_2_string;
typedef vector<Mat_2_string> Mat_3_string;
typedef vector<Mat_3_string> Mat_4_string;

complex<double> const One_Complex(1.0,0.0);
complex<double> const Iota_Complex(0.0,1.0);
complex<double> const Zero_Complex(0.0,0.0);

struct pair_int{
    int first;
    int second;

    // Constructor accepting two integers
    pair_int(int a, int b) : first(a), second(b) {}

    // Overload the < operator for sorting
    bool operator<(const pair_int &other) const {
        if (first != other.first)
            return first < other.first;
        return second < other.second;
    }

    // Overload the == operator for uniqueness
    bool operator==(const pair_int &other) const {
        return (first == other.first) && (second == other.second);
    }
};

typedef vector<pair_int> Mat_1_intpair;
typedef vector<Mat_1_intpair> Mat_2_intpair;

typedef pair<double, double> pair_doub;
typedef vector<pair_doub> Mat_1_doubpair;

#endif