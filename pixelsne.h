#ifndef PIXELSNE_H
#define PIXELSNE_H

#include "ptree.h"

static inline double sign(double x) { return (x == .0 ? .0 : (x < .0 ? -1.0 : 1.0)); }

class PixelSNE
{
private:
    PTree* tree;
public:
    PixelSNE();
    ~PixelSNE();
    void run(double* X, int N, int D, double* Y, int no_dims, double perplexity, double theta,
        unsigned int bins, int p_method, int rand_seed, bool skip_random_init, int max_iter=1000, int stop_lying_iter=250, 
        int mom_switch_iter=250);
    bool load_data(double** data, int* n, int* d, int* no_dims, double* theta, double* perplexity, unsigned int* bins, int* p_method, int* rand_seed);
    void save_data(double* data, int* landmarks, double* costs, int n, int d);
    void symmetrizeMatrix(unsigned long long** row_P, unsigned long long** col_P, double** val_P, int N); // should be static!

private:
    void computeGradient(double* P, unsigned long long* inp_row_P, unsigned long long* inp_col_P, double* inp_val_P, double* Y, int N, int D, double* dC, double theta, double beta, unsigned int bins, int iter_cnt);
    void computeExactGradient(double* P, double* Y, int N, int D, double* dC);
    double evaluateError(double* P, double* Y, int N, int D);
    double evaluateError(unsigned long long* row_P, unsigned long long* col_P, double* val_P, double* Y, int N, int D, double theta, double beta, unsigned int bins, int iter_cnt);
    void zeroMean(double* X, int N, int D);
    double minmax(double* X, int N, int D, double beta, unsigned int bins, int iter_cnt);
    void computeGaussianPerplexity(double* X, int N, int D, double* P, double perplexity);
    void computeGaussianPerplexity(double* X, int N, int D, unsigned long long** _row_P, unsigned long long** _col_P, double** _val_P, double perplexity, int K);
    void computeSquaredEuclideanDistance(double* X, int N, int D, double* DD);
    double randn();
};

#endif
