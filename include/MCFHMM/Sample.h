#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <glog/logging.h>
using namespace std;

#ifndef SAMPLE_H
#define SAMPLE_H

class Sample{

public:
    ~Sample();

    vector<double> values;
    double p;

    void init_rand(vector<double> *low_limit, vector<double> *high_limit);
    Sample combine(vector<double> second);
    size_t size();
    string str();
};

#endif
