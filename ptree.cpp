#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include "ptree.h"

extern long long num_insert;

// Constructs cell
Cell::Cell(unsigned int inp_dimension) {
    dimension = inp_dimension;
    corner = (int*) malloc(dimension * sizeof(int));
    width  = (int*) malloc(dimension * sizeof(int));
}

Cell::Cell(unsigned int inp_dimension, int* inp_corner, int* inp_width) {
    dimension = inp_dimension;
    corner = (int*) malloc(dimension * sizeof(int));
    width  = (int*) malloc(dimension * sizeof(int));
    for(int d = 0; d < dimension; d++) setCorner(d, inp_corner[d]);
    for(int d = 0; d < dimension; d++) setWidth(d, inp_width[d]);
}

// Destructs cell
Cell::~Cell() {
    free(corner);
    free(width);
}

int Cell::getCorner(unsigned int d) {
    return corner[d];
}

int Cell::getWidth(unsigned int d) {
    return width[d];
}

void Cell::setCorner(unsigned int d, int val) {
    corner[d] = val;
}

void Cell::setWidth(unsigned int d, int val) {
    width[d] = val;
}

// Checks whether a point lies in a cell
bool Cell::containsPoint(double point[])
{
    for(int d = 0; d < dimension; d++) {
        if(point[d] <= corner[d] - width[d]) return false;
        if(corner[d] + width[d] < point[d]) return false;
    }
    return true;
}

PTree::PTree(unsigned int D, double* inp_data, unsigned int N, unsigned int bins, int lv, int iter_cnt)
{

    num_insert = 0;

    // #ifdef USE_BITWISE_OP
    //     printf("ptree.cpp USE_BITWISE_OP\n");
    // #else
    //     printf("ptree.cpp not USE_BITWISE_OP\n");
    // #endif

    // Compute mean, width, and height of current map (boundaries of SPTree)
    int nD = 0;
    int* mean_Y = (int*) calloc(D,  sizeof(int));

    #ifdef USE_BITWISE_OP
        for (int d = 0; d < D; d++) mean_Y[d] = (bins >> 1) - 1;/*Op*/
    #else
        for(int d = 0; d < D; d++) mean_Y[d] = (bins/2) - 1;/*Op*/
    #endif

    // Construct PTree
    int* width = (int*) malloc(D * sizeof(int));

    #ifdef USE_BITWISE_OP
        for (int d = 0; d < D; d++) width[d] = (bins >> 1);/*Op*/
    #else
        for (int d = 0; d < D; d++) width[d] = (bins /2);/*Op*/
    #endif

    init(NULL, D, inp_data, mean_Y, width, bins, lv, iter_cnt);
    fill(N, iter_cnt);    // fill every data.
    
    // Clean up memory
    free(mean_Y);
    free(width);
}

// Constructor for PTree with particular size and parent -- build the tree, too!
PTree::PTree(unsigned int D, double* inp_data, unsigned int N, int* inp_corner, int* inp_width, unsigned int pixel_width, int lv, int iter_cnt)
{
    init(NULL, D, inp_data, inp_corner, inp_width, pixel_width, lv, iter_cnt);
    fill(N, iter_cnt);
}


// Constructor for PTree with particular size (do not fill the tree)
PTree::PTree(unsigned int D, double* inp_data, int* inp_corner, int* inp_width, unsigned int pixel_width, int lv, int iter_cnt)
{
    init(NULL, D, inp_data, inp_corner, inp_width, pixel_width, lv, iter_cnt);
}


// Constructor for PTree with particular size and parent (do not fill tree)
PTree::PTree(PTree* inp_parent, unsigned int D, double* inp_data, int* inp_corner, int* inp_width, unsigned int pixel_width, int lv, int iter_cnt) {
    init(inp_parent, D, inp_data, inp_corner, inp_width, pixel_width, lv, iter_cnt);
}


// Constructor for PTree with particular size and parent -- build the tree, too!
PTree::PTree(PTree* inp_parent, unsigned int D, double* inp_data, unsigned int N, int* inp_corner, int* inp_width, unsigned int pixel_width, int lv, int iter_cnt)
{
    init(inp_parent, D, inp_data, inp_corner, inp_width, pixel_width, lv, iter_cnt);
    fill(N, iter_cnt);
}

// Main initialization function
void PTree::init(PTree* inp_parent, unsigned int D, double* inp_data, int* inp_corner, int* inp_width, unsigned int pix_width, int lv, int iter_cnt)
{
    iter_count = iter_cnt;
    level = lv;
    parent = inp_parent;
    dimension = D;
    no_children=(1<<D); 
    data = inp_data;
    size = 0;
    cum_size = 0;

    pixel_width = pix_width; 

    //printf("tree level: %d, pixel_width: %d\n", level, pixel_width);

    if (pixel_width == 1) {
        //printf("pixel_width = 1, level: %d\n", level);
        is_leaf = true;
    } else {
        is_leaf = false;
    }
    
    boundary = new Cell(dimension);
    for(unsigned int d = 0; d < D; d++) boundary->setCorner(d, inp_corner[d]);
    for(unsigned int d = 0; d < D; d++) boundary->setWidth(d, inp_width[d]);
    
    children = (PTree**) malloc(no_children * sizeof(PTree*));
    for(unsigned int i = 0; i < no_children; i++) children[i] = NULL;
    
    buff = (double*) malloc(D * sizeof(double));
    for(unsigned int i = 0; i < D; i++) buff[i] = .0;
}

void PTree::clean(int iter_cnt) {

    size = 0;
    cum_size = 0;
    iter_count = iter_cnt;
}

// Destructor for PTree
PTree::~PTree()
{
    for(unsigned int i = 0; i < no_children; i++) {
        if(children[i] != NULL) delete children[i];
    }

    free(children);
    free(buff);
    delete boundary;
}


// Update the data underlying this tree
void PTree::setData(double* inp_data)
{
    data = inp_data;
}


// Get the parent of the current tree
PTree* PTree::getParent()
{
    return parent;
}


// Insert a point into the PTree
bool PTree::insert(unsigned int new_index, int iter_cnt)
{

    if (is_leaf) {
        size++;
        num_insert++;
        return true;
    } 

    // Ignore objects which do not belong in this quad tree
    double* point = data + new_index * dimension;
 
    if(!boundary->containsPoint(point)) {
        return false;
    }

    num_insert++;
    cum_size++;
    
    bool hasPoint = false;
    for (int i = 0 ; i < no_children ; i++) {       
        if (children[i] != NULL){
            hasPoint = true;
            break;
        }
    }
    if (hasPoint == false){
        subdivide();
    } 

    bool isOldTree = false;
    for (int i = 0 ; i < no_children ; i++) {       
        if (children[i]->iter_count < iter_count) {
            isOldTree = true;
            continue;
        }
    }
    if (isOldTree == true) {
        for(unsigned int i = 0; i < no_children; i++) {
            children[i]->clean(iter_cnt);
        }
    }
    // Find out where the point can be inserted
    for(unsigned int i = 0; i < no_children; i++) {
        if(children[i]->insert(new_index, iter_cnt)) return true;
    }

    // Otherwise, the point cannot be inserted (this should never happen)
    return false;
}

    
// Create four children which fully divide this cell into four quads of equal area
void PTree::subdivide() {
    
    // Create new children
    int* new_corner = (int*) malloc(dimension * sizeof(int));
    int* new_width  = (int*) malloc(dimension * sizeof(int));

    #ifdef USE_BITWISE_OP
    	unsigned int new_pixel_width = (pixel_width >> 1); /*Op*/
    #else
        unsigned int new_pixel_width = (pixel_width / 2); /*Op*/
    #endif

    // printf("subdivide, level: %d\n", level);

    for(unsigned int i = 0; i < no_children; i++) {
        unsigned int div = 1;
        for(unsigned int d = 0; d < dimension; d++) {
            new_width[d] = .5 * boundary->getWidth(d);

            #ifdef USE_BITWISE_OP
    			if (((i / div) & 1) == 1) new_corner[d] = boundary->getCorner(d) - .5 * boundary->getWidth(d); /*Op*/
            #else
                if (((i / div) % 2) == 1) new_corner[d] = boundary->getCorner(d) - .5 * boundary->getWidth(d); /*Op*/
            #endif
			
            else
                new_corner[d] = boundary->getCorner(d) + .5 * boundary->getWidth(d);

            #ifdef USE_BITWISE_OP
    			div = (div << 1);/*Op*/
            #else
                div = (div * 2);/*Op*/
            #endif
		}

        children[i] = new PTree(this, dimension, data, new_corner, new_width, new_pixel_width, level+1, iter_count);
    }
    free(new_corner);
    free(new_width);
    
    // Empty parent node
    size = 0;
    is_leaf = false;
}

// Build PTree on dataset
void PTree::fill(unsigned int N, int iter_cnt)
{
    for(unsigned int i = 0; i < N; i++) insert(i, iter_cnt);
}

unsigned int PTree::getDepth() {
    if(is_leaf) return 1;
    int depth = 0;
    for(unsigned int i = 0; i < no_children; i++) depth = fmax(depth, children[i]->getDepth());
    return 1 + depth;
}

// Compute non-edge forces using Barnes-Hut algorithm
void PTree::computeNonEdgeForces(unsigned int point_index, double theta, double neg_f[], double* sum_Q, double beta, int iter_cnt)
{
    
    // Make sure that we spend no time on empty nodes or self-interactions
    if(cum_size == 0 || iter_count < iter_cnt) 
        return;

    if (is_leaf == true) {
        double* point = data + point_index * dimension;
        if (boundary->containsPoint(point))
            return;
    }
    
    // Compute distance between point and center-of-mass
    double D = .0;
    unsigned int ind = point_index * dimension;
    for(unsigned int d = 0; d < dimension; d++) {buff[d] = data[ind + d] - boundary->getCorner(d); }
    for(unsigned int d = 0; d < dimension; d++) D += buff[d] * buff[d];

    int max_width = 0;
    int cur_width;
    for(unsigned int d = 0; d < dimension; d++) {
        cur_width = boundary->getWidth(d);
        max_width = (max_width > cur_width) ? max_width : cur_width;
    }

    #ifdef USE_BITWISE_OP
    	if (is_leaf || (D!=0 && (max_width*max_width / D < theta*theta))) {/*Op*/
    #else
        if (is_leaf || (D!=0 && (max_width / sqrt(D) < theta))) {/*Op*/
    #endif
        
        // Compute and add t-SNE force between point and current node
        D = beta / (beta + D);
        double mult = cum_size * D;
        *sum_Q += mult;
        mult *= D;
        for(unsigned int d = 0; d < dimension; d++) neg_f[d] += mult * buff[d];
    }
    else {
        // Recursively apply Barnes-Hut to children
        for(unsigned int i = 0; i < no_children; i++) {
            children[i]->computeNonEdgeForces(point_index, theta, neg_f, sum_Q, beta, iter_cnt);
        }
    }
}


// Computes edge forces
void PTree::computeEdgeForces(unsigned long long* row_P, unsigned long long* col_P, double* val_P, int N, double* pos_f, double beta)
{
    // Loop over all edges in the graph
    unsigned long long ind1 = 0;
    unsigned long long ind2 = 0;
    double D;
    for(unsigned int n = 0; n < N; n++) {
        for(unsigned long long i = row_P[n]; i < row_P[n + 1]; i++) {
        
            // Compute pairwise distance and Q-value
            D = beta;
            ind2 = col_P[i] * dimension;
            for(unsigned int d = 0; d < dimension; d++) buff[d] = data[ind1 + d] - data[ind2 + d];
            for(unsigned int d = 0; d < dimension; d++) D += buff[d] * buff[d];
            D = val_P[i] * beta / D;

            // Sum positive force
            for(unsigned int d = 0; d < dimension; d++) pos_f[ind1 + d] += D * buff[d];
        }
        ind1 += dimension;
    }
}

// Print out tree

void PTree::print() 
{
    // if(cum_size == 0) {
    //     printf("[%d]Empty node\n", level);
    //     return;
    // }

    //printf("level: %d, is_leaf: %d\n", level, is_leaf);

    if (level == 0)
        printf("total data size: %d\n", cum_size);

    if(is_leaf) {
        printf("Leaf node: %d, size: %d\n", level, size);
        // printf("[%d]Leaf node; data = [", level);
        // for(int i = 0; i < size; i++) {
        //     if(i < size - 1) printf("%.3f, ", data[i]);
        //     else printf("]\n");
        // }
    }
    else {
        // for(int d = 0; d < dimension; d++) {
        // }
        for(int i = 0; i < no_children; i++) 
            if (children[i] != NULL)
                children[i]->print();
    }
}























