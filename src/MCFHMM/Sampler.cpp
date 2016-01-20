#include "Sampler.h"
#include <tuple>
#include <random>
#include <chrono>
using namespace std;

Sample Sampler::likelihood_weighted_sampler(vector<Sample> &sample_set){
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine gen(seed);
    uniform_real_distribution<double> dist(0, 1.0);

    vector<double> low_p;
    vector<double> high_p;

    low_p.push_back(0.0);
    high_p.push_back(sample_set[0].p);
    for (size_t i = 1; i < sample_set.size(); i++){
        low_p.push_back(high_p[i - 1]);
        high_p.push_back(high_p[i - 1] + sample_set[i].p);
    }

    double temp_p = dist(gen);
    int index = 0;
    for (size_t j = 0; j < low_p.size(); j++){
        if (temp_p >= low_p[j] && temp_p <= high_p[j]){
            index = j;
            break;
        }
    }

    return sample_set[index];
}

vector<Sample> Sampler::uniform_sampling(vector<double> *sample_low, vector<double> *sample_high, size_t N){
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine gen(seed);

    vector<Sample> results;

    for (size_t i = 0; i < N; i++){
        Sample sample;
        for (size_t j = 0; j < sample_low->size(); j++){
            uniform_real_distribution<double> dist((*sample_low)[j], (*sample_high)[j]);
            sample.values.push_back(dist(gen));
        }
        results.push_back(sample);
    }

    return results;
}

vector<Sample> Sampler::likelihood_weighted_resampler(vector<Sample> &sample_set, int size){
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine gen(seed);
    uniform_real_distribution<double> dist(0, 1.0);

    if (size == -1){
        size = sample_set.size();
    }

    vector<double> low_p;
    vector<double> high_p;

    low_p.push_back(0.0);
    high_p.push_back(sample_set[0].p);
    for (size_t i = 1; i < sample_set.size(); i++){
        low_p.push_back(high_p[i - 1]);
        high_p.push_back(high_p[i - 1] + sample_set[i].p);
    }

    vector<Sample> temp;
    for (int i = 0; i < size; i++){
        double temp_p = dist(gen);
        int index = 0;
        for (size_t j = 0; j < low_p.size(); j++){
            if (temp_p >= low_p[j] && temp_p <= high_p[j]){
                index = j;
                break;
            }
        }

        if (index >= (int)sample_set.size()) index = sample_set.size() - 1;
        temp.push_back(sample_set[index]);
    }

    return temp;
}

Sample Sampler::sample(DETree *tree){
    Sample sample;
    sample.init_rand(&tree->get_root()->node_lower_bounds, &tree->get_root()->node_higher_bounds);
    sample.p = 1.0;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine gen(seed);
    DETreeNode *node = tree->get_root();

    bool cond = !node->leaf_node;
    while(cond){
        int dimension = node->node_split_dimension;
        double min_value = node->node_lower_bounds[dimension];
        double max_value = node->node_higher_bounds[dimension];

        if (node->leaf_node){
            uniform_real_distribution<double> dist(min_value, max_value);
            double temp = dist(gen);
            sample.values[dimension] = temp;

            cond = !node->leaf_node;
        }else{
            uniform_real_distribution<double> child_decision_gen(0, node->left_child->node_sigma + node->right_child->node_sigma);
            double decision = child_decision_gen(gen);

            if (decision < node->left_child->node_sigma){
                min_value = node->left_child->node_lower_bounds[dimension];
                max_value = node->left_child->node_higher_bounds[dimension];
            }else{
                min_value = node->right_child->node_lower_bounds[dimension];
                max_value = node->right_child->node_higher_bounds[dimension];
            }

            uniform_real_distribution<double> dist(min_value, max_value);
            double temp = dist(gen);
            sample.values[dimension] = temp;

            cond = !node->leaf_node;

            if (sample.values[dimension] < node->node_mid_point){
                node = node->left_child;
            }else{
                node = node->right_child;
            }
        }
    }

    return sample;
}

Sample Sampler::sample_given(DETree *tree, Sample &given){
    Sample sample;
    sample.init_rand(&tree->get_root()->node_lower_bounds, &tree->get_root()->node_higher_bounds);
    sample.p = 1.0;

    int offset = sample.size() - given.size();
    for (size_t i = offset; i < sample.size(); i++){
        sample.values[i] = given.values[i - offset];
    }

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine gen(seed);
    DETreeNode *node = tree->get_root();

    bool cond = !node->leaf_node;
    while(cond){
        int dimension = node->node_split_dimension;
        double min_value = node->node_lower_bounds[dimension];
        double max_value = node->node_higher_bounds[dimension];

        // If we need to actually generate values
        if (dimension < offset){
            if (node->leaf_node){
                uniform_real_distribution<double> dist(min_value, max_value);
                double temp = dist(gen);
                sample.values[dimension] = temp;

                cond = !node->leaf_node;
            }else{
                uniform_real_distribution<double> child_decision_gen(0, node->left_child->node_f_hat + node->right_child->node_f_hat);
                double decision = child_decision_gen(gen);  // This can be moved out of the while loop to make the process faster!

                if (decision < node->left_child->node_f_hat){
                    min_value = node->left_child->node_lower_bounds[dimension];
                    max_value = node->left_child->node_higher_bounds[dimension];
                }else{
                    min_value = node->right_child->node_lower_bounds[dimension];
                    max_value = node->right_child->node_higher_bounds[dimension];
                }

                uniform_real_distribution<double> dist(min_value, max_value);
                double temp = dist(gen);
                sample.values[dimension] = temp;

                cond = !node->leaf_node;

                if (sample.values[dimension] < node->node_mid_point){
                    node = node->left_child;
                }else{
                    node = node->right_child;
                }
            }
        }else{ // If we have the values given
            if (!node->leaf_node && sample.values[dimension] < node->node_mid_point){
                node = node->left_child;
            }else if (!node->leaf_node){
                node = node->right_child;
            }
        }
    }

    return sample;
}

vector<Sample> Sampler::resample_from(DETree *tree, size_t sample_set_size){
    vector<Sample> results;
    for (size_t i = 0; i < sample_set_size; i++){
        Sample smp = sample(tree);
        smp.p = 1.0 / sample_set_size;
        results.push_back(smp);
    }
    return results;
}
