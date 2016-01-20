#include <vector>
#include <glog/logging.h>
#include "Sample.h"
using namespace std;

#ifndef DETREE_H
#define DETREE_H

class DETreeNode;

class DETree{

public:
    DETree();
    DETree(vector<Sample> &sample_set, vector<double> *sample_low, vector<double> *sample_high);
    ~DETree();

    void create_tree(const vector<Sample> &sample_set, vector<double> *sample_low, vector<double> *sample_high);
    DETree create_conditional_tree(Sample given);
    double density_value(Sample sample, double rho);

    DETreeNode* get_root();

    vector<DETreeNode*>* depth_first();
    string depth_first_str();

private:
    DETreeNode *root;

    void depth_first(vector<DETreeNode*> *& nodes, DETreeNode* current_node);
};

class DETreeNode{

public:
    DETreeNode();
    DETreeNode(vector<Sample>sub_sample, int level, char node_type, vector<double> low_bounds, vector<double> high_bounds);
    ~DETreeNode();

    string str();

    const double min_diff_interval = 0.05;

    vector<Sample> samples;

    bool leaf_node = false;
    int level = 0;
    int node_size = 0;
    char node_type = 'R';
    double node_sigma = 0.0;
    double node_f_hat = 0.0;
    double node_v_size = 1.0;

    vector<double> node_lower_bounds;
    vector<double> node_higher_bounds;

    double node_longest_interval = 0.0;
    int node_split_dimension = 0;
    double node_mid_point = 0.0;

    DETreeNode *left_child;
    DETreeNode *right_child;
    DETreeNode *parent;
};

#endif
