#ifndef PTREE_H
#define PTREE_H

using namespace std;

class Cell {

    unsigned int dimension; // number of dimension
    int* corner;         // position of corner
    int* width;          // width. height is same with width.
    
public:
    Cell(unsigned int inp_dimension);
    Cell(unsigned int inp_dimension, int* inp_corner, int* inp_width);
    ~Cell();
    
    int getCorner(unsigned int d);
    int getWidth(unsigned int d);

    void setCorner(unsigned int d, int val);
    void setWidth(unsigned int d, int val);

    bool containsPoint(double point[]);
};

class PTree
{

 public:   

	// A buffer we use when doing force computations
    double* buff;
    
    // Properties of this node in the tree
    PTree* parent;
    unsigned int dimension;
    
    unsigned int size;
    unsigned int cum_size;
        
    // Axis-aligned bounding box stored as a center with half-dimensions to represent the boundaries of this quad tree
    Cell* boundary;
    
    // Indices in this space-partitioning tree node, corresponding center-of-mass, and list of all children
    unsigned int data_size;
    double* data;

    PTree** children;
    unsigned int no_children;

    unsigned int pixel_width;

    bool is_leaf;
    int level;
    int iter_count;

public:
    PTree(unsigned int D, double* inp_data, unsigned int N, unsigned int bins, int lv, int iter_cnt);
    PTree(unsigned int D, double* inp_data, unsigned int N, int* inp_corner, int* inp_width, unsigned int pixel_width, int lv, int iter_cnt);
    PTree(unsigned int D, double* inp_data, int* inp_corner, int* inp_width, unsigned int pixel_width, int lv, int iter_cnt);
    PTree(PTree* inp_parent, unsigned int D, double* inp_data, int* inp_corner, int* inp_width, unsigned int pixel_width, int lv, int iter_cnt);
    PTree(PTree* inp_parent, unsigned int D, double* inp_data, unsigned int N, int* inp_corner, int* inp_width, unsigned int pixel_width, int lv, int iter_cnt);
    ~PTree();

    void clean(int iter_cnt);
    void setData(double* inp_data);
    PTree* getParent();
    void construct(Cell boundary);
    bool insert(unsigned int new_index, int iter_cnt);
    void subdivide();

    unsigned int getDepth();
    void computeNonEdgeForces(unsigned int point_index, double theta, double neg_f[], double* sum_Q, double beta, int iter_cnt);
    void computeEdgeForces(unsigned long long* row_P, unsigned long long* col_P, double* val_P, int N, double* pos_f, double beta);
    void print(int level);
    void fill(unsigned int N, int iter_cnt);

    
private:
    void init(PTree* inp_parent, unsigned int D, double* inp_data, int* inp_corner, int* inp_width, unsigned int pix_width, int lv, int iter_cnt);
    bool isChild(unsigned int test_index, unsigned int start, unsigned int end);
};



#endif