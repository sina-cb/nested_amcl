#include "DETree.h"
#include <iostream>
#include <limits>
#include <cmath>
#include <algorithm>
#include <string>
#include <sstream>
using namespace std;

DETree::DETree(){
    LOG(INFO) << "The DETree is not created, you should use create_tree explicitly to initialize the DETree";
}

DETree::DETree(vector<Sample> &sample_set, vector<double> *sample_low, vector<double> *sample_high){
    create_tree(sample_set, sample_low, sample_high);
}

DETree::~DETree(){
    vector<DETreeNode*>* nodes = depth_first();
    for (size_t i = 0; i < nodes->size(); i++){
        delete (*nodes)[i];
    }

//    LOG(INFO) << "DETREE DECONSTRUCTOR!";
}

double DETree::density_value(Sample sample, double rho){

    DETreeNode *node = this->get_root();

    bool cond = !node->leaf_node;
    while(cond){
        int dimension = node->node_split_dimension;

        cond = !node->leaf_node;

        if (cond && sample.values[dimension] < node->node_mid_point){
            node = node->left_child;
        }else if (cond){
            node = node->right_child;
        }
    }

    cond = !(node->node_type == 'R');
    double p = 0.0;
    double temp_rho = 1.0;
    while (cond){
        if (node->node_type != 'R'){
            p = p + (1 - rho) * temp_rho * node->node_f_hat;
        } else {
            p = p + temp_rho * node->node_f_hat;
        }
        temp_rho = temp_rho * rho;

        cond = !(node->node_type == 'R');
        node = node->parent;
    }

    return p;

}

void DETree::create_tree(const vector<Sample> &sample_set, vector<double> *sample_low, vector<double> *sample_high){
    root = new DETreeNode(sample_set, 0, 'R', *sample_low, *sample_high);
}

string DETree::depth_first_str(){
    stringbuf buf;
    ostream os (&buf);
    vector<DETreeNode*> * nodes = depth_first();

    double f_hat_sum = 0.0;
    for (size_t i = 0; i < nodes->size(); i++){
        os << "Level: " << (*nodes)[i]->level << "\t";
        os << "Side: " << (*nodes)[i]->node_type << "\t";
        os << "Size: " << (*nodes)[i]->node_size << "\t";
        os << "Sigma: " << (*nodes)[i]->node_sigma << "\t";
        os << "F Hat: " << (*nodes)[i]->node_f_hat << "\t";
        os << "V Size: " << (*nodes)[i]->node_v_size << "\t";
        os << "\n";

        f_hat_sum += (*nodes)[i]->node_f_hat;
    }

    os << "\n" << "F Hat Sum: " << f_hat_sum << endl;

    return buf.str();
}

vector<DETreeNode*>* DETree::depth_first(){
    vector<DETreeNode*>* results = new vector<DETreeNode*>();
    depth_first(results, root);
    return results;
}

void DETree::depth_first(vector<DETreeNode *> *& nodes, DETreeNode *current_node){
    if (current_node->leaf_node){
        nodes->push_back(current_node);
        return;
    }

    nodes->push_back(current_node);

    if (current_node->left_child != NULL){
        depth_first(nodes, current_node->left_child);
    }

    if (current_node->right_child != NULL){
        depth_first(nodes, current_node->right_child);
    }
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

DETreeNode* DETree::get_root(){
    return root;
}

DETreeNode::DETreeNode(){

}

DETreeNode::DETreeNode(vector<Sample> sub_sample, int level, char node_type, vector<double> low_bounds, vector<double> high_bounds){

    // Set the current node's sample set
    for (size_t i = 0; i < sub_sample.size(); i++){
        samples.push_back(sub_sample[i]);
        node_sigma = node_sigma + sub_sample[i].p;
    }

    // Set the variables in the Node
    node_size = samples.size();
    this->node_type = node_type;
    this->level = level;
    node_lower_bounds = low_bounds;
    node_higher_bounds = high_bounds;

    low_bounds.clear();
    high_bounds.clear();

    // Compute the node's region size
    for (size_t i = 0; i < node_lower_bounds.size(); i++){
        node_v_size *= (node_higher_bounds[i] - node_lower_bounds[i]);
    }

    node_f_hat = node_sigma / node_v_size;

    // Check if we have a leaf node or not
    if (samples.size() <= 1){
        leaf_node = true;
        return;
    }

    // Find the longest dimension
    node_longest_interval = node_higher_bounds[0] - node_lower_bounds[0];
    for (size_t i = 1; i < node_lower_bounds.size(); i++){
        if (node_higher_bounds[i] - node_lower_bounds[i] > node_longest_interval){
            node_longest_interval = node_higher_bounds[i] - node_lower_bounds[i];
            node_split_dimension = i;
        }
    }

    // Compute the new bounds
    vector<double> left_lower_intervals;
    vector<double> left_higher_intervals;

    vector<double> right_lower_intervals;
    vector<double> right_higher_intervals;

    node_mid_point = (node_longest_interval / 2.0) + node_lower_bounds[node_split_dimension];

    for (int i = 0; i < (int)node_lower_bounds.size(); i++){
        if (i == node_split_dimension){
            left_lower_intervals.push_back(node_lower_bounds[i]);
            left_higher_intervals.push_back(node_mid_point);
            right_lower_intervals.push_back(node_mid_point);
            right_higher_intervals.push_back(node_higher_bounds[i]);
        }else{
            left_lower_intervals.push_back(node_lower_bounds[i]);
            left_higher_intervals.push_back(node_higher_bounds[i]);
            right_lower_intervals.push_back(node_lower_bounds[i]);
            right_higher_intervals.push_back(node_higher_bounds[i]);
        }
    }

    // Sort the samples based on the value in the dimension we want to split
    std::sort(samples.begin(), samples.end(),
              [&, this](const Sample& a, const Sample& b) {
                    return a.values[node_split_dimension] < b.values[node_split_dimension];
                }
    );

    // Split the samples into two sub sample sets
    int cut_index = 0;
    for (size_t i = 0; i < samples.size(); i++){
        if (samples[i].values[node_split_dimension] < node_mid_point){
            cut_index = i;
        }
    }

    if (cut_index == (int)samples.size() - 1){
        cut_index--;
    }

    vector<Sample> sub_sample_left (
                samples.begin(),
                samples.begin() + cut_index + 1
                );
    vector<Sample> sub_sample_right(
                samples.begin() + cut_index + 1,
                samples.end()
                );

    this->left_child = new DETreeNode(sub_sample_left, level + 1, 'l', left_lower_intervals, left_higher_intervals);
    this->right_child = new DETreeNode(sub_sample_right, level + 1, 'r', right_lower_intervals, right_higher_intervals);

    this->left_child->parent = this;
    this->right_child->parent = this;
}

DETreeNode::~DETreeNode(){
    samples.clear();
    node_lower_bounds.clear();
    node_higher_bounds.clear();

//    LOG(INFO) << "DETREENODE DECONSTRUCTOR!!!";
}

string DETreeNode::str(){
    stringbuf buf;
    ostream os(&buf);
    for (size_t i = 0; i < samples.size(); i++){
        for (size_t j = 0; j < samples[0].size(); j++){
            os << samples[i].values[j] << "\t";
        }
        os << "P: " << samples[i].p << endl;
    }

    return buf.str();
}
